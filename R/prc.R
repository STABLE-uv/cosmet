#' Principal response curves
#'
#' `fit_prc()` fits the principal response curves (PRC) of van den Brink and
#' ter Braak (1999) with [vegan::prc()]: a partial redundancy analysis that
#' removes the changes over time shared by all mesocosms and summarises, on
#' its first axis, how each treatment departs from the control at each time.
#'
#' @param x A [meso_data()] object, usually after [transform_response()].
#' @param taxa `NULL` (default, all) or the taxa to include.
#' @param scaling Scaling of the coefficients and taxon weights, as in
#'   [vegan::scores.rda()]: `"symmetric"` (default), `"sites"`, `"species"`
#'   or `"none"`. Only the relative size of coefficients and weights changes.
#' @param orientation `"highest_down"` (default) changes the sign of the axis,
#'   which is arbitrary, so that the coefficient of largest absolute value of
#'   the highest treatment is negative; `"none"` keeps the sign of vegan.
#'
#' @return An object of class `meso_prc`, a list with
#'   * `curves`: one row per treatment and time (control included, with
#'     `cdt = 0`), columns `treatment`, `concentration`, `time` and `cdt`, the
#'     PRC coefficient;
#'   * `weights`: one row per taxon (C-locale order), columns `taxon` and
#'     `bk`, the taxon weight. The product \eqn{b_k c_{dt}} approximates the
#'     difference between treatment and control on the scale of the analyses;
#'     positive weights follow the curves, negative ones go the opposite way;
#'   * `variance`: a one-row data frame with the total variance, the
#'     variance explained by time and by treatment (main effect and
#'     interaction with time), as values and as % of the total, the
#'     eigenvalue of the first axis and its % of the treatment variance and
#'     of the total;
#'   * `model`: the [vegan::prc()] object; `samples`: the samples in the order
#'     of the rows of the model; `settings`.
#'
#' @details
#' The model is \eqn{Y \sim treatment + treatment:time + Condition(time)},
#' with the control as reference level (vegan uses treatment contrasts). The
#' coefficients \eqn{c_{dt}} and weights \eqn{b_k} are those of
#' `summary(model, scaling = scaling)`.
#'
#' @section Missing values and zeros:
#' The analysis needs a complete table: an unknown value (`NA`) in any sample
#' is an error. Use `zero_fill = TRUE` in [meso_data()] if missing taxa are
#' absences, or remove the samples or taxa with unknown values.
#'
#' @references
#' van den Brink PJ, ter Braak CJF (1999). Principal response curves:
#' analysis of time-dependent multivariate responses of biological community
#' to stress. *Environmental Toxicology and Chemistry* 18, 138–148.
#' \doi{10.1002/etc.5620180207}
#' @seealso [test_prc()], [plot_prc()], [community_noec()].
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' prc <- fit_prc(transform_response(d, A = 2))
#' prc
#' head(prc$weights[order(-abs(prc$weights$bk)), ])
fit_prc <- function(x, taxa = NULL, scaling = "symmetric", orientation = "highest_down") {
  check_meso_data(x)
  check_choice(scaling, c("symmetric", "sites", "species", "none"))
  check_choice(orientation, c("highest_down", "none"))
  m <- community_matrix(x, taxa)
  smp <- m$samples
  trt <- factor(smp$treatment, levels = x$settings$treatments)
  if (nlevels(droplevels(trt)) < 2L) {
    abort_input("The PRC needs the control and at least one other treatment.")
  }
  times <- sort(unique(smp$time))
  if (length(times) < 2L) {
    abort_input("The PRC needs at least two sampling times; use an RDA for a single time.")
  }
  time <- factor(smp$time, levels = times, labels = format(times, trim = TRUE))
  model <- vegan::prc(response = m$y, treatment = trt, time = time)
  s <- summary(model, axis = 1, scaling = scaling)
  cdt <- s$coefficients
  bk <- drop(s$sp)
  if (is.null(cdt) || length(cdt) == 0L || is.null(bk) || length(bk) == 0L) {
    abort_internal("PRC coefficients or taxon weights could not be extracted from vegan::prc().")
  }
  highest <- nrow(cdt)
  dominant <- cdt[highest, which.max(abs(cdt[highest, ]))]
  reversed <- orientation == "highest_down" && isTRUE(dominant > 0)
  if (reversed) {
    cdt <- -cdt
    bk <- -bk
  }
  treated <- rownames(cdt)
  conc <- x$design$concentration[match(x$settings$treatments, x$design$treatment)]
  curves <- data.frame(
    treatment = rep(c(x$settings$control, treated), times = length(times)),
    time = rep(times, each = length(treated) + 1L),
    cdt = as.vector(rbind(0, cdt)),
    stringsAsFactors = FALSE
  )
  curves$concentration <- conc[match(curves$treatment, x$settings$treatments)]
  curves <- curves[c("treatment", "concentration", "time", "cdt")]
  weights <- data.frame(taxon = colnames(m$y), bk = unname(bk[colnames(m$y)]),
                        stringsAsFactors = FALSE)
  total <- model$tot.chi
  var_time <- if (is.null(model$pCCA)) 0 else model$pCCA$tot.chi
  var_trt <- model$CCA$tot.chi
  eig1 <- unname(model$CCA$eig[1L])
  variance <- data.frame(
    total = total, time = var_time, treatment = var_trt,
    residual = model$CA$tot.chi,
    time_pct = 100 * var_time / total, treatment_pct = 100 * var_trt / total,
    axis1 = eig1, axis1_pct_treatment = 100 * eig1 / var_trt,
    axis1_pct_total = 100 * eig1 / total
  )
  settings <- list(scaling = scaling, orientation = orientation, reversed = reversed,
                   transform = x$settings$transform, control = x$settings$control,
                   treatments = x$settings$treatments, concentrations = conc,
                   concentration_unit = x$settings$concentration_unit,
                   applications = x$settings$applications, taxa = colnames(m$y))
  structure(list(curves = curves, weights = weights, variance = variance, model = model,
                 samples = smp, response = m$y, settings = settings),
            class = "meso_prc")
}

# Samples x taxa matrix on the scale of the analyses, rows ordered by
# mesocosm and time, columns by taxon (C locale). No unknown values allowed.
community_matrix <- function(x, taxa, call = sys.call(-1)) {
  d <- select_taxa(x, taxa, call = call)
  bad <- unique(d$taxon[is.na(d$response)])
  if (length(bad) > 0L) {
    abort_input(sprintf(
      paste(
        "Multivariate analyses need a complete table, but %d taxa have unknown",
        "values (%s). If a missing taxon is an absence, use zero_fill = TRUE in",
        "meso_data(); otherwise leave those taxa out with `taxa`."
      ),
      length(bad), format_items(bad)
    ), call = call)
  }
  smp <- unique(d[c("mesocosm", "time", "treatment", "concentration")])
  smp <- smp[order(smp$mesocosm, smp$time, method = "radix"), , drop = FALSE]
  rownames(smp) <- NULL
  taxa_all <- sort(unique(d$taxon), method = "radix")
  y <- matrix(NA_real_, nrow(smp), length(taxa_all),
              dimnames = list(NULL, taxa_all))
  r <- match(paste(d$mesocosm, format_time_key(d$time), sep = "\u001f"),
             paste(smp$mesocosm, format_time_key(smp$time), sep = "\u001f"))
  y[cbind(r, match(d$taxon, taxa_all))] <- d$response
  list(y = y, samples = smp)
}

#' @export
print.meso_prc <- function(x, ...) {
  v <- x$variance
  cat("<meso_prc> principal response curves\n")
  cat(sprintf("%d samples, %d taxa; scale: %s; scaling: %s%s\n", nrow(x$samples),
              nrow(x$weights), describe_transform(x$settings$transform), x$settings$scaling,
              if (x$settings$reversed) " (axis reversed)" else ""))
  cat(sprintf("Variance explained by time: %.1f%%, by treatment: %.1f%%\n",
              v$time_pct, v$treatment_pct))
  cat(sprintf("First axis: %.1f%% of the treatment variance (%.1f%% of the total)\n",
              v$axis1_pct_treatment, v$axis1_pct_total))
  w <- x$weights[order(-abs(x$weights$bk)), , drop = FALSE]
  cat("Taxa with the largest weights:\n")
  print(utils::head(w, 8L), row.names = FALSE, digits = 3)
  invisible(x)
}

#' Test the principal response curves
#'
#' `test_prc()` tests whether the first PRC axis explains more of the
#' treatment effect than expected by chance, by permuting the complete time
#' series of the mesocosms among treatments, and, for every sampling time,
#' whether the treatments differ (redundancy analysis of that time).
#'
#' @param x A `meso_prc` object from [fit_prc()].
#' @param n_perm Number of permutations (default 999).
#' @param by_time If `TRUE` (default), also test every sampling time.
#' @param seed Seed of the permutations (default 123). With a number, the
#'   result is reproducible and the global random number stream is restored
#'   afterwards; with `NULL`, the current stream is used.
#'
#' @return An object of class `meso_prc_test`, a list with `global` (one row:
#'   `df`, `variance` of the first axis, `F`, `p_value`, `n_perm`),
#'   `by_time` (one row per time: `time`, `n`, `r_squared`,
#'   `adj_r_squared`, `F`, `p_value`, `n_perm`; `NULL` with `by_time =
#'   FALSE`) and `settings`.
#'
#' @details
#' Global test: F statistic of the first constrained axis of the PRC model
#' ([vegan::anova.cca()] with `first = TRUE`), with the mesocosms as whole
#' plots permuted freely and no permutation within them, so that each time
#' series stays intact. Tests by time: RDA of the samples of that time on the
#' treatment, with free permutations.
#'
#' @section Assumptions:
#' The global test needs every mesocosm to be sampled at every time; an
#' unbalanced design is an error.
#'
#' @section Reproducibility:
#' p-values have a Monte Carlo error of about
#' \eqn{\sqrt{p (1 - p) / n_{perm}}}{sqrt(p (1 - p) / n_perm)}. The same seed
#' can give different permutations in other versions of R or 'permute'.
#'
#' @references
#' van den Brink PJ, ter Braak CJF (1999). Principal response curves:
#' analysis of time-dependent multivariate responses of biological community
#' to stress. *Environmental Toxicology and Chemistry* 18, 138–148.
#' \doi{10.1002/etc.5620180207}
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' prc <- fit_prc(transform_response(d, A = 2))
#' test_prc(prc, n_perm = 199)
test_prc <- function(x, n_perm = 999, by_time = TRUE, seed = 123) {
  if (!inherits(x, "meso_prc")) {
    abort_input("`x` must be the result of fit_prc().")
  }
  check_number(n_perm, lower = 1, upper = 1e6)
  if (n_perm != round(n_perm)) abort_input("`n_perm` must be a whole number.")
  check_flag(by_time)
  check_number(seed, allow_null = TRUE)
  smp <- x$samples
  times_per <- table(smp$mesocosm)
  if (length(unique(as.integer(times_per))) != 1L ||
      nrow(smp) != length(unique(smp$time)) * length(times_per)) {
    abort_input(paste(
      "The permutation of whole time series needs every mesocosm sampled at",
      "every time; check_design() lists the missing samples."
    ))
  }
  run <- function(expr) if (is.null(seed)) expr else withr::with_seed(seed, expr)
  plots <- factor(smp$mesocosm, levels = sort(unique(smp$mesocosm), method = "radix"))
  ctrl <- permute::how(plots = permute::Plots(strata = plots, type = "free"),
                       within = permute::Within(type = "none"), nperm = n_perm)
  a <- run(vegan::anova.cca(x$model, permutations = ctrl, first = TRUE))
  global <- data.frame(df = a$Df[1L], variance = a$Variance[1L], F = a$F[1L],
                       p_value = a[["Pr(>F)"]][1L], n_perm = n_perm)
  by <- NULL
  if (by_time) {
    trt_all <- factor(smp$treatment, levels = x$settings$treatments)
    rows <- lapply(sort(unique(smp$time)), function(t) {
      idx <- smp$time == t
      trt <- droplevels(trt_all[idx])
      if (nlevels(trt) < 2L) return(NULL)
      y <- x$response[idx, , drop = FALSE]
      m_t <- vegan::rda(y ~ trt)
      a_t <- run(vegan::anova.cca(m_t, permutations = n_perm))
      r2 <- vegan::RsquareAdj(m_t)
      data.frame(time = t, n = sum(idx), r_squared = r2$r.squared,
                 adj_r_squared = r2$adj.r.squared, F = a_t$F[1L],
                 p_value = a_t[["Pr(>F)"]][1L], n_perm = n_perm)
    })
    by <- do.call(rbind, rows)
    rownames(by) <- NULL
  }
  structure(list(global = global, by_time = by,
                 settings = list(n_perm = n_perm, seed = seed,
                                 permutations = "whole time series of mesocosms, free; none within")),
            class = "meso_prc_test")
}

#' @export
print.meso_prc_test <- function(x, ...) {
  g <- x$global
  cat("<meso_prc_test>\n")
  cat(sprintf("First PRC axis: F = %.3f, p = %.4g (%d permutations of whole time series)\n",
              g$F, g$p_value, g$n_perm))
  if (!is.null(x$by_time)) {
    cat("Tests by time (RDA on treatment):\n")
    print(x$by_time, row.names = FALSE, digits = 3)
  }
  invisible(x)
}

#' Community NOECs
#'
#' `community_noec()` derives one NOEC per sampling time for the whole
#' community: for every time, a principal component analysis of the samples
#' of that time, and Williams' test (or Dunnett's) on the scores of its
#' first axis, as in Amador et al. (2024).
#'
#' @param x A [meso_data()] object, usually after [transform_response()].
#' @param method `"williams"` (default, as Amador et al. 2024) or
#'   `"dunnett"`.
#' @param pre_times `NULL` or the sampling times before the application (see
#'   [derive_noec()]).
#' @param taxa `NULL` (default, all) or the taxa to include.
#' @param alpha,step_variance,n_sim,seed As in [test_control()].
#'
#' @return An object of class `meso_noec` (see [derive_noec()]) with the
#'   taxon `"PC1"`, plus `pca`, a data frame with one row per time: `time`,
#'   `n`, `pc1_pct` (% of the variance of that time on the first axis) and
#'   `status`. Both directions are tested, because the sign of a principal
#'   component is arbitrary.
#'
#' @section Missing values and zeros:
#' As in [fit_prc()], the table must be complete.
#'
#' @references
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @seealso [derive_noec()], [fit_prc()].
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' cn <- community_noec(transform_response(d, A = 2), pre_times = -7, n_sim = 9999)
#' noec_table(cn)
#' cn$pca
community_noec <- function(x, method = "williams", pre_times = NULL, taxa = NULL,
                           alpha = 0.05, step_variance = "all", n_sim = 99999, seed = 123) {
  check_meso_data(x)
  check_choice(method, c("williams", "dunnett"))
  m <- community_matrix(x, taxa)
  smp <- m$samples
  rows <- list()
  pca <- list()
  for (t in sort(unique(smp$time))) {
    idx <- which(smp$time == t)
    y <- m$y[idx, , drop = FALSE]
    keep <- apply(y, 2L, stats::var) > 0
    if (!any(keep)) {
      pca[[length(pca) + 1L]] <- data.frame(time = t, n = length(idx), pc1_pct = NA_real_,
                                            status = "zero_variance")
      score <- rep(NA_real_, length(idx))
    } else {
      p <- vegan::rda(y[, keep, drop = FALSE])
      score <- as.vector(vegan::scores(p, display = "sites", choices = 1, scaling = "sites"))
      pca[[length(pca) + 1L]] <- data.frame(time = t, n = length(idx),
                                            pc1_pct = unname(100 * p$CA$eig[1L] / p$tot.chi),
                                            status = "ok")
    }
    rows[[length(rows) + 1L]] <- data.frame(mesocosm = smp$mesocosm[idx], time = t,
                                            taxon = "PC1", value = score,
                                            stringsAsFactors = FALSE)
  }
  scores <- do.call(rbind, rows)
  scores <- merge_design(scores, x)
  d <- new_meso_data(scores, x$design, x$samples,
                     utils::modifyList(x$settings, list(unit = NULL,
                       transform = list(method = "none", A = NA_real_))),
                     list())
  tests <- test_control(d, method = method, direction = "both", alpha = alpha,
                        step_variance = step_variance, n_sim = n_sim, seed = seed)
  out <- derive_noec(tests, pre_times = pre_times)
  # The sign of a principal component is arbitrary: an "increase" of the
  # scores is not an increase of the community.
  dir <- out$noec$direction
  out$noec$direction[dir %in% "increase"] <- "higher_scores"
  out$noec$direction[dir %in% "decrease"] <- "lower_scores"
  out$pca <- do.call(rbind, pca)
  rownames(out$pca) <- NULL
  out$settings$community <- "first principal component of each time"
  out
}

# Treatment and concentration of each mesocosm, keeping rows and order.
merge_design <- function(d, x) {
  i <- match(d$mesocosm, x$design$mesocosm)
  data.frame(mesocosm = d$mesocosm, time = d$time, treatment = x$design$treatment[i],
             concentration = x$design$concentration[i], taxon = d$taxon, value = d$value,
             response = d$value, value_status = "recorded", stringsAsFactors = FALSE)
}
