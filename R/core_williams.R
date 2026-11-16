# Numerical core of Williams' test (internal) ---------------------------------
#
# Williams DA (1971) Biometrics 27:103-117; Williams DA (1972) Biometrics
# 28:519-531. Two versions of the step-down procedure:
#   * "all" (default): the statistic of the step that tests dose i uses the
#     maximum-likelihood (isotonic) estimate of mu_i computed from all the
#     doses and the variance pooled over all groups; its critical value is
#     that of Williams' distribution for i doses with the degrees of freedom
#     of that variance (Williams 1971, 1972).
#   * "reduced": each step re-analyses the design without the higher doses,
#     with the isotonic estimate and the pooled variance of the control and
#     doses 1..i only (the original population_noec() script of the package
#     authors).
# The null distribution is simulated with the algorithm of that script.

# Maximum-likelihood mean of the highest dose under the order restriction
# mu_1 <= ... <= mu_i (control excluded): the largest weighted mean of any
# set of consecutive treated groups that ends at the highest dose.
williams_last_iso <- function(mu, n) {
  max(cumsum(rev(mu * n)) / cumsum(rev(n)))
}

# Maximum-likelihood means of all treated groups under mu_1 <= ... <= mu_k
# (weighted isotonic regression), with the max-min formula of Williams (1971):
# M_i = max_{u <= i} min_{v >= i} sum_{j=u}^{v} n_j y_j / sum_{j=u}^{v} n_j.
williams_iso_all <- function(mu, n) {
  k <- length(mu)
  vapply(seq_len(k), function(i) {
    max(vapply(seq_len(i), function(u) {
      min(vapply(i:k, function(v) sum(n[u:v] * mu[u:v]) / sum(n[u:v]), numeric(1)))
    }, numeric(1)))
  }, numeric(1))
}

# Simulated null distribution of the statistic for the replication `n`
# (control first, at least two doses) and residual degrees of freedom `df`,
# sorted. Normal errors with a common variance: dose means N(0, 1/n_j),
# residual SD sqrt(chi2_df / df). The order of the random draws is that of
# the original script, so that the same seed gives the same distribution.
williams_null <- function(n, n_sim, df = sum(n - 1)) {
  k <- length(n)
  z0 <- stats::rnorm(n_sim, sd = 1 / sqrt(n[1L]))
  weighted_sum <- numeric(n_sim)
  n_pool <- 0
  ymax <- rep(-Inf, n_sim)
  for (i in seq.int(k, 2L)) {
    weighted_sum <- weighted_sum + stats::rnorm(n_sim, sd = sqrt(n[i]))
    n_pool <- n_pool + n[i]
    ymax <- pmax(ymax, weighted_sum / n_pool)
  }
  residual_sd <- sqrt(stats::rchisq(n_sim, df = df) / df)
  sort((ymax - z0) / (residual_sd * sqrt(1 / n[k] + 1 / n[1L])))
}

# One null distribution per replication pattern and degrees of freedom. With
# a seed, every pattern is simulated from that seed, so a result does not
# depend on which other taxa or dates are analysed, nor in which order.
williams_null_cached <- function(n, df, n_sim, seed, cache) {
  key <- paste(c(n, "df", df), collapse = "_")
  out <- cache[[key]]
  if (is.null(out)) {
    out <- if (is.null(seed)) {
      williams_null(n, n_sim, df)
    } else {
      withr::with_seed(seed, williams_null(n, n_sim, df))
    }
    assign(key, out, envir = cache)
  }
  out
}

# Pooled SD, degrees of freedom and standard error of the step that tests
# dose i (group i + 1; control = group 1).
williams_scale <- function(n, ss, i, step_variance) {
  ix <- if (step_variance == "reduced") seq_len(i) else seq_along(n)
  df <- sum(n[ix] - 1)
  s <- sqrt(sum(ss[ix]) / df)
  if (!is.finite(s) || s <= 0) {
    return(list(df = df, s = NA_real_, se = NA_real_))
  }
  list(df = df, s = s, se = s * sqrt(1 / n[i] + 1 / n[1L]))
}

# One-sided p-value of the statistic of a step with groups `n` (control and
# the doses up to the one tested).
williams_p <- function(stat, n, df, n_sim, seed, cache) {
  if (length(n) == 2L) {
    return(stats::pt(stat, df, lower.tail = FALSE))
  }
  null <- williams_null_cached(n, df, n_sim, seed, cache)
  n_below <- findInterval(stat, null, left.open = TRUE)
  (1 + length(null) - n_below) / (length(null) + 1)
}

# One-sided critical value at level `a`.
williams_critical <- function(n, df, a, n_sim, seed, cache) {
  if (length(n) == 2L) {
    return(stats::qt(1 - a, df))
  }
  stats::quantile(williams_null_cached(n, df, n_sim, seed, cache), 1 - a, type = 1,
                  names = FALSE)
}

# Step-down procedure for one direction: from the highest dose down, stopping
# at the first non-rejection. `sign` is 1 for increases and -1 for decreases.
# Returns one row per dose (lowest first).
williams_path <- function(mu, n, ss, sign, n_dir, alpha, n_sim, seed, cache,
                          step_variance) {
  k <- length(n)
  tested <- rep(FALSE, k - 1L)
  stat <- rep(NA_real_, k - 1L)
  p <- rep(NA_real_, k - 1L)
  p_adj <- rep(NA_real_, k - 1L)
  rejected <- rep(NA, k - 1L)
  iso_all <- if (step_variance == "all") williams_iso_all(sign * mu[-1L], n[-1L]) else NULL
  for (i in seq.int(k, 2L)) {
    j <- i - 1L
    sc <- williams_scale(n, ss, i, step_variance)
    tested[j] <- TRUE
    if (!is.na(sc$s)) {
      m_i <- if (step_variance == "reduced") {
        williams_last_iso(sign * mu[2:i], n[2:i])
      } else {
        iso_all[j]
      }
      stat[j] <- (m_i - sign * mu[1L]) / sc$se
      p[j] <- williams_p(stat[j], n[seq_len(i)], sc$df, n_sim, seed, cache)
      p_adj[j] <- min(1, n_dir * p[j])
      rejected[j] <- p_adj[j] < alpha
    }
    if (!isTRUE(rejected[j])) break
  }
  list(tested = tested, statistic = stat, p_value = p, p_adjusted = p_adj,
       significant = rejected)
}
