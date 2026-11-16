#' Compare every treatment with the control
#'
#' `test_control()` compares, for every taxon (or variable) and sampling time,
#' each treatment with the control on the scale of the analyses (see
#' [transform_response()]), with Williams' step-down test or with Dunnett's
#' test. Its result feeds [derive_noec()] and [calc_mdd()].
#'
#' @param x A [meso_data()] object.
#' @param method `"williams"` or `"dunnett"`; no default. Williams' test
#'   assumes that the response is monotonic in the concentration and needs
#'   concentrations; Dunnett's test does not assume a monotonic response.
#' @param direction `"both"` (default), `"decrease"` or `"increase"`: the
#'   effects looked for.
#' @param alpha Significance level for the procedure of each taxon and time
#'   (default 0.05, the level used by Amador et al. 2024).
#' @param taxa `NULL` (default, all) or a character vector with the taxa to
#'   analyse, for example those marked as eligible by [taxon_eligibility()].
#' @param na_rm If `FALSE` (default), a taxon and time with an unknown value
#'   (`NA`) is not tested (`status = "missing_values"`). If `TRUE`, the
#'   unknown values are left out and the test uses the remaining replicates
#'   (`n_missing` says how many were left out).
#' @param step_variance Williams' test only: `"all"` (default) uses, at every
#'   step of the procedure, the isotonic means estimated from all the doses
#'   and the variance pooled over all groups (Williams 1971, 1972);
#'   `"reduced"` re-analyses at each step the design without the higher
#'   doses, with the isotonic mean and pooled variance of the control and the
#'   doses up to the one tested (the original analysis script of the package
#'   authors). The two can give different NOECs: `"all"` has more degrees of
#'   freedom; `"reduced"` is not affected by the variance of the higher doses.
#' @param n_sim Number of simulated draws of the null distribution of
#'   Williams' statistic (default 99999). Not used by Dunnett's test.
#' @param seed Seed of the simulation of Williams' null distribution
#'   (default 123). With a number, results are reproducible and the global
#'   random number stream is restored afterwards; every replication pattern
#'   is simulated from this seed, so the result of a taxon does not depend on
#'   the other taxa analysed. `NULL` uses the current stream.
#'
#' @return An object of class `meso_test`, a list with
#'   * `tests`: one row per taxon, time, treatment and direction tested, in
#'     the order taxon (C locale), time, direction and concentration (or
#'     treatment), with columns `taxon`, `time`, `treatment`,
#'     `concentration`, `direction` (`"decrease"`, `"increase"`, or `"both"`
#'     for the two-sided Dunnett test), `n_control`, `n_treatment`,
#'     `mean_control_value` and `mean_treatment_value` (arithmetic means of
#'     the values, in the unit of the data), `mean_control` and
#'     `mean_treatment` (means on the scale of the analyses), `difference`
#'     (treatment minus control, scale of the analyses), `se` and `df` (its
#'     standard error and degrees of freedom), `critical_value` (of the test
#'     statistic at the level used), `lower` and `upper` (confidence interval
#'     of the difference, level `1 - alpha`), `ci_type` (`"simultaneous"` for
#'     Dunnett, `"marginal"` for Williams), `tested`, `statistic` (oriented so
#'     that large values are evidence for the direction of the row; signed for
#'     `"both"`), `p_value` (Williams: one-sided p-value of the step;
#'     Dunnett: single-step adjusted p-value), `p_adjusted` (Williams:
#'     `p_value` times the number of directions, at most 1; Dunnett: equal to
#'     `p_value`) and `significant` (`p_adjusted < alpha`; `NA` if not
#'     tested);
#'   * `cells`: one row per taxon and time with `status` (`"ok"`,
#'     `"absent"` (all values zero), `"no_data"`, `"missing_values"`,
#'     `"insufficient_replicates"` (a group with fewer than two values) or
#'     `"zero_residual_variance"`), `reason`, `n_values` and `n_missing`;
#'   * `settings`: the method and options used, the transformation and the
#'     order of the treatments.
#'
#' @details
#' **Williams.** The step that tests dose \eqn{i} uses
#' \deqn{\bar t_i = \frac{\hat M_i - \bar y_0}{s \sqrt{1/n_i + 1/n_0}}}{t_i =
#'   (M_i - ybar_0) / (s sqrt(1/n_i + 1/n_0))}
#' where \eqn{\hat M_i} is the maximum-likelihood estimate of the mean of dose
#' \eqn{i} under a monotonic response. With `step_variance = "all"`,
#' \eqn{\hat M_i = \max_{u \le i} \min_{v \ge i} \sum_{j=u}^{v} n_j \bar y_j /
#' \sum_{j=u}^{v} n_j}{M_i = max over u <= i of min over v >= i of the weighted
#' mean of doses u to v} over all doses and \eqn{s^2} is the variance pooled
#' over all groups (Williams 1971); with `"reduced"`, \eqn{\hat M_i} is the
#' largest weighted mean of doses \eqn{u, \dots, i} and \eqn{s^2} the variance
#' pooled over the control and doses \eqn{1, \dots, i}. For decreases the
#' means change sign. The test starts at the highest dose and stops at the
#' first non-significant one. The null distribution of \eqn{\bar t_i}
#' (Williams' distribution for \eqn{i} doses with the degrees of freedom of
#' \eqn{s}; normal errors, common variance) is simulated with `n_sim` draws
#' and the actual replication; with a single dose it is Student's t. With
#' `direction = "both"` the two one-sided procedures are run and their
#' p-values doubled (Bonferroni).
#'
#' **Dunnett.** \eqn{t_i = (\bar y_i - \bar y_0) / (s \sqrt{1/n_i + 1/n_0})}
#' with the variance pooled over all groups; adjusted p-values and critical
#' values come from the multivariate t distribution of the \eqn{t_i}
#' (two-sided: their maximum absolute value), computed by numerical
#' integration. Confidence intervals are simultaneous.
#'
#' @section Missing values and zeros:
#' Zeros are observed values. Unknown values stop the test of that taxon and
#' time unless `na_rm = TRUE`. A taxon whose values are all zero at a time is
#' `"absent"` and not tested.
#'
#' @section Assumptions:
#' Values are independent between mesocosms and, on the scale of the
#' analyses, normal with a common variance in all treatments. Williams' test
#' also assumes that the mean response is monotonic in the concentration.
#'
#' @section Reproducibility:
#' Williams' p-values have a Monte Carlo error of about
#' \eqn{\sqrt{p (1 - p) / n_{sim}}}{sqrt(p (1 - p) / n_sim)}. The same seed
#' can give different draws in other versions of R.
#'
#' @references
#' Williams DA (1971). A test for differences between treatment means when
#' several dose levels are compared with a zero dose control. *Biometrics*
#' 27, 103–117. \doi{10.2307/2528930}
#'
#' Williams DA (1972). The comparison of several dose levels with a zero dose
#' control. *Biometrics* 28, 519–531. \doi{10.2307/2556164}
#'
#' Dunnett CW (1955). A multiple comparison procedure for comparing several
#' treatments with a control. *Journal of the American Statistical
#' Association* 50, 1096–1121. \doi{10.1080/01621459.1955.10501294}
#'
#' @seealso [derive_noec()], [calc_mdd()].
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" &
#'                   habitat == "eutrophic_lake" & taxon %in% c("Chironomini", "Cloeon"))
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' d <- transform_response(d, A = 2)
#' tw <- test_control(d, method = "williams", n_sim = 9999)
#' tw
#' td <- test_control(d, method = "dunnett")
#' head(as.data.frame(td))
test_control <- function(x, method, direction = "both", alpha = 0.05, taxa = NULL,
                         na_rm = FALSE, step_variance = "all", n_sim = 99999, seed = 123) {
  check_meso_data(x)
  if (missing(method)) {
    abort_input(paste(
      "`method` is missing and has no default: use \"williams\" (assumes a",
      "monotonic response) or \"dunnett\" (does not)."
    ))
  }
  check_choice(method, c("williams", "dunnett"))
  check_choice(direction, c("both", "decrease", "increase"))
  check_number(alpha, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  check_flag(na_rm)
  check_choice(step_variance, c("all", "reduced"))
  check_number(n_sim, lower = 999, upper = 1e7)
  if (n_sim != round(n_sim)) abort_input("`n_sim` must be a whole number.")
  check_number(seed, allow_null = TRUE)
  d <- select_taxa(x, taxa)
  levels_trt <- x$settings$treatments
  conc <- x$design$concentration[match(levels_trt, x$design$treatment)]
  if (method == "williams") {
    check_dose_order(x, levels_trt, conc)
  }
  directions <- if (direction == "both") c("decrease", "increase") else direction
  res <- run_tests(d, levels_trt, conc, method = method, direction = direction,
                   directions = directions, alpha = alpha, na_rm = na_rm,
                   step_variance = step_variance, n_sim = n_sim, seed = seed)
  settings <- list(
    method = method, direction = direction, alpha = alpha, na_rm = na_rm,
    n_sim = if (method == "williams") n_sim else NA_real_,
    seed = if (method == "williams") seed else NULL,
    step_variance = if (method == "williams") step_variance else NA_character_,
    transform = x$settings$transform, unit = x$settings$unit,
    concentration_unit = x$settings$concentration_unit,
    control = x$settings$control, treatments = levels_trt, concentrations = conc,
    applications = x$settings$applications, taxa = sort(unique(d$taxon), method = "radix")
  )
  structure(list(tests = res$tests, cells = res$cells, settings = settings),
            class = "meso_test")
}

# Rows of the taxa requested.
select_taxa <- function(x, taxa, call = sys.call(-1)) {
  d <- x$data
  if (is.null(taxa)) return(d)
  if (!is.character(taxa) || length(taxa) == 0L || anyNA(taxa)) {
    abort_input("`taxa` must be NULL or a character vector of taxa.", call = call)
  }
  unknown <- setdiff(taxa, d$taxon)
  if (length(unknown) > 0L) {
    abort_input(sprintf("Taxa not in the data: %s.", format_items(quote_items(unknown))),
                call = call)
  }
  d[d$taxon %in% taxa, , drop = FALSE]
}

check_dose_order <- function(x, levels_trt, conc, call = sys.call(-1)) {
  if (!isTRUE(x$settings$has_concentration)) {
    abort_input(paste(
      "Williams' test needs the concentration of each treatment: give",
      "`concentration_col` in meso_data(), or use method = \"dunnett\"."
    ), call = call)
  }
  if (anyDuplicated(conc)) {
    abort_input(sprintf(
      paste(
        "Williams' test needs one treatment per concentration; these share a",
        "concentration: %s. Use method = \"dunnett\" or analyse them separately."
      ),
      format_items(levels_trt[conc %in% conc[duplicated(conc)]])
    ), call = call)
  }
  invisible(TRUE)
}

# Group statistics of one taxon and time. `g` is the index of the treatment
# level of every value (control = 1).
prepare_cell <- function(y, value, g, n_groups, na_rm) {
  n_missing <- sum(is.na(y))
  if (n_missing > 0L && !na_rm) {
    return(list(status = "missing_values",
                reason = sprintf("%d unknown value(s)", n_missing),
                n_missing = n_missing, n_values = length(y)))
  }
  keep <- !is.na(y)
  y <- y[keep]
  value <- value[keep]
  g <- factor(g[keep], levels = seq_len(n_groups))
  yy <- split(y, g)
  vv <- split(value, g)
  n <- unname(lengths(yy))
  # mean() (not sum / n) so that equal values give a residual of exactly zero
  mu <- unname(vapply(yy, function(v) if (length(v) > 0L) mean(v) else NA_real_, numeric(1)))
  mean_value <- unname(vapply(vv, function(v) if (length(v) > 0L) mean(v) else NA_real_,
                              numeric(1)))
  ss <- unname(vapply(yy, function(v) if (length(v) > 0L) sum((v - mean(v))^2) else 0,
                      numeric(1)))
  status <- "ok"
  reason <- NA_character_
  if (length(y) == 0L) {
    status <- "no_data"
    reason <- "no values"
  } else if (all(value == 0)) {
    status <- "absent"
    reason <- "all values are zero"
  } else if (any(n < 2L)) {
    status <- "insufficient_replicates"
    reason <- "a treatment has fewer than 2 values"
  } else if (sum(ss) == 0) {
    status <- "zero_residual_variance"
    reason <- "all replicates of every treatment are equal"
  }
  list(status = status, reason = reason, n = n, mu = mu, ss = ss,
       mean_value = mean_value, n_missing = n_missing, n_values = length(y))
}

run_tests <- function(d, levels_trt, conc, method, direction, directions, alpha,
                      na_rm, step_variance, n_sim, seed) {
  n_groups <- length(levels_trt)
  cells <- unique(d[c("taxon", "time")])
  cells <- cells[order(cells$taxon, cells$time, method = "radix"), , drop = FALSE]
  rownames(cells) <- NULL
  cell_key <- paste(d$taxon, format_time_key(d$time), sep = "\u001f")
  rows_of <- split(seq_len(nrow(d)), factor(cell_key, levels = unique(
    paste(cells$taxon, format_time_key(cells$time), sep = "\u001f"))))
  cache <- new.env(parent = emptyenv())
  out_tests <- vector("list", nrow(cells))
  out_cells <- vector("list", nrow(cells))
  for (ci in seq_len(nrow(cells))) {
    i <- rows_of[[ci]]
    g <- match(d$treatment[i], levels_trt)
    prep <- prepare_cell(d$response[i], d$value[i], g, n_groups, na_rm)
    out_cells[[ci]] <- data.frame(taxon = cells$taxon[ci], time = cells$time[ci],
                                  status = prep$status, reason = prep$reason,
                                  n_values = prep$n_values, n_missing = prep$n_missing,
                                  stringsAsFactors = FALSE)
    tab <- if (method == "williams") {
      williams_cell(prep, n_groups, directions, alpha, step_variance, n_sim, seed, cache)
    } else {
      dunnett_cell(prep, n_groups, direction, alpha, cache)
    }
    tab <- cbind(data.frame(taxon = rep(cells$taxon[ci], nrow(tab)),
                            time = rep(cells$time[ci], nrow(tab)),
                            stringsAsFactors = FALSE), tab)
    out_tests[[ci]] <- tab
  }
  tests <- do.call(rbind, out_tests)
  tests$treatment <- levels_trt[tests$group]
  tests$concentration <- conc[tests$group]
  tests <- tests[c("taxon", "time", "treatment", "concentration", "direction",
                   setdiff(names(tests), c("taxon", "time", "treatment", "concentration",
                                           "direction", "group")))]
  rownames(tests) <- NULL
  cells_df <- do.call(rbind, out_cells)
  rownames(cells_df) <- NULL
  list(tests = tests, cells = cells_df)
}

# Empty rows (cell not testable) with the columns of a test table.
blank_rows <- function(groups, directions, prep) {
  n_rows <- length(groups) * length(directions)
  na <- rep(NA_real_, n_rows)
  grp <- rep(groups, times = length(directions))
  dirs <- rep(directions, each = length(groups))
  has_n <- !is.null(prep$n)
  data.frame(
    group = grp, direction = dirs,
    n_control = if (has_n) rep(prep$n[1L], n_rows) else rep(NA_integer_, n_rows),
    n_treatment = if (has_n) prep$n[grp] else rep(NA_integer_, n_rows),
    mean_control_value = if (has_n) rep(prep$mean_value[1L], n_rows) else na,
    mean_treatment_value = if (has_n) prep$mean_value[grp] else na,
    mean_control = if (has_n) rep(prep$mu[1L], n_rows) else na,
    mean_treatment = if (has_n) prep$mu[grp] else na,
    difference = if (has_n) prep$mu[grp] - prep$mu[1L] else na,
    se = na, df = na, critical_value = na, lower = na, upper = na,
    ci_type = rep(NA_character_, n_rows), tested = rep(FALSE, n_rows),
    statistic = na, p_value = na, p_adjusted = na, significant = rep(NA, n_rows),
    stringsAsFactors = FALSE
  )
}

williams_cell <- function(prep, n_groups, directions, alpha, step_variance, n_sim, seed,
                          cache) {
  groups <- seq.int(2L, n_groups)
  tab <- blank_rows(groups, directions, prep)
  if (!identical(prep$status, "ok")) return(tab)
  n_dir <- length(directions)
  a <- alpha / n_dir
  # Scale of the step that tests each dose (whether or not it is reached).
  for (i in groups) {
    sc <- williams_scale(prep$n, prep$ss, i, step_variance)
    rows <- tab$group == i
    tab$se[rows] <- sc$se
    tab$df[rows] <- sc$df
    if (!is.na(sc$s)) {
      tab$critical_value[rows] <- williams_critical(prep$n[seq_len(i)], sc$df, a, n_sim,
                                                    seed, cache)
      half <- stats::qt(1 - alpha / 2, sc$df) * sc$se
      tab$lower[rows] <- tab$difference[rows] - half
      tab$upper[rows] <- tab$difference[rows] + half
    }
    tab$ci_type[rows] <- "marginal"
  }
  for (dir in directions) {
    sign <- if (dir == "increase") 1 else -1
    path <- williams_path(prep$mu, prep$n, prep$ss, sign, n_dir, alpha, n_sim, seed, cache,
                          step_variance)
    rows <- which(tab$direction == dir)
    tab$tested[rows] <- path$tested
    tab$statistic[rows] <- path$statistic
    tab$p_value[rows] <- path$p_value
    tab$p_adjusted[rows] <- path$p_adjusted
    tab$significant[rows] <- path$significant
  }
  tab
}

dunnett_cell <- function(prep, n_groups, direction, alpha, cache) {
  groups <- seq.int(2L, n_groups)
  tab <- blank_rows(groups, direction, prep)
  if (!identical(prep$status, "ok")) return(tab)
  n <- prep$n
  df <- sum(n - 1)
  s <- sqrt(sum(prep$ss) / df)
  se <- s * sqrt(1 / n[-1L] + 1 / n[1L])
  t <- (prep$mu[-1L] - prep$mu[1L]) / se
  lambda <- sqrt(n[-1L] / (n[-1L] + n[1L]))
  two_sided <- direction == "both"
  crit <- dunnett_critical(lambda, df, alpha, two_sided, cache)
  stat <- switch(direction, both = t, increase = t, decrease = -t)
  p <- dunnett_p(stat, lambda, df, two_sided)
  diff <- tab$difference
  tab$se <- se
  tab$df <- rep(df, length(se))
  tab$critical_value <- rep(crit, length(se))
  tab$lower <- switch(direction, both = diff - crit * se, increase = diff - crit * se,
                      decrease = rep(-Inf, length(se)))
  tab$upper <- switch(direction, both = diff + crit * se, increase = rep(Inf, length(se)),
                      decrease = diff + crit * se)
  tab$ci_type <- rep("simultaneous", length(se))
  tab$tested <- rep(TRUE, length(se))
  tab$statistic <- stat
  tab$p_value <- p
  tab$p_adjusted <- p
  tab$significant <- p < alpha
  tab
}

#' @export
print.meso_test <- function(x, ...) {
  s <- x$settings
  cells <- x$cells
  cat(sprintf("<meso_test> %s, direction: %s, alpha = %s\n",
              if (s$method == "williams") "Williams' test" else "Dunnett's test", s$direction,
              format(s$alpha)))
  cat(sprintf("Scale: %s; %d taxa x times\n", describe_transform(s$transform), nrow(cells)))
  if (s$method == "williams") {
    cat(sprintf("Steps: %s; null distribution: %s simulated draws, seed %s\n",
                if (s$step_variance == "all") {
                  "variance and isotonic means of all groups"
                } else {
                  "reduced design (groups up to the dose tested)"
                },
                format(s$n_sim, big.mark = ","), if (is.null(s$seed)) "none" else format(s$seed)))
  }
  st <- table(factor(cells$status, levels = unique(c("ok", sort(unique(cells$status))))))
  cat("Status of the taxa x times: ",
      paste(sprintf("%s %d", names(st), as.integer(st)), collapse = ", "), "\n", sep = "")
  sig <- x$tests[x$tests$significant %in% TRUE, , drop = FALSE]
  n_cells <- nrow(unique(sig[c("taxon", "time")]))
  cat(sprintf("Taxa x times with a significant difference from the control: %d\n", n_cells))
  invisible(x)
}

#' @export
as.data.frame.meso_test <- function(x, ...) {
  x$tests
}
