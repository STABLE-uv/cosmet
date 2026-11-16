#' Minimum detectable differences
#'
#' `calc_mdd()` computes, for every taxon, time and treatment, the minimum
#' detectable difference (MDD) of the test used in [test_control()]: the
#' smallest difference from the control that the test would have declared
#' significant, with the same variance, degrees of freedom and critical value.
#' It is returned with the estimated effect and its confidence interval,
#' because a non-significant result with a small MDD is not proof that there
#' is no effect.
#'
#' @param x A `meso_test` object from [test_control()].
#'
#' @return A data frame with one row per taxon, time and treatment (order of
#'   `x$tests`), with columns `taxon`, `time`, `treatment`, `concentration`,
#'   `status` (of the taxon and time in [test_control()]), `n_control`,
#'   `n_treatment`, `control_value` (control mean on the original scale: the
#'   back-transformed mean \eqn{(e^{\bar y_0} - 1)/A} with
#'   [transform_response()], the arithmetic mean otherwise),
#'   `difference`, `lower`, `upper` and `ci_type` (on the scale of the
#'   analyses, from [test_control()]), `effect_pct` (difference between the
#'   treatment and the control on the original scale, % of the control),
#'   `se`, `df`, `critical_value`, `mdd` (scale of the analyses), `mdd_value`
#'   (original scale), `mdd_pct` (% of the control) and `mdd_class` (`"0"`,
#'   `"I"`, `"II"`, `"III"` or `"IV"`). Values are `NA` where the taxon and
#'   time could not be tested or the control mean is zero.
#'
#' @details
#' Following Brock et al. (2015, p. 1162, eq. 2):
#' \deqn{MDD = t_{1-\alpha, df, k} \, s \sqrt{1/n_0 + 1/n}}{MDD = t * s *
#'   sqrt(1/n_0 + 1/n)}
#' where \eqn{t} is the critical value of the test (Student's t with a single
#' treatment, Williams' critical value of the step, or Dunnett's critical
#' value; with `direction = "both"`, at level `alpha / 2` for Williams'
#' test), and \eqn{s} the pooled standard deviation of the same test.
#' Without transformation, \eqn{MDD\% = 100 \, MDD / \bar x_0} (eq. 3). With
#' \eqn{y = \ln(A x + 1)} (eq. 5 and 6):
#' \deqn{MDD_{abu} = \frac{e^{\bar y_0} - e^{\bar y_0 - MDD}}{A}, \quad
#'   MDD\% = 100 \, \frac{MDD_{abu}}{(e^{\bar y_0} - 1)/A}}{MDD_abu =
#'   (exp(y0) - exp(y0 - MDD)) / A; MDD\% = 100 MDD_abu / ((exp(y0) - 1) / A)}
#' with \eqn{\bar y_0} the mean of the transformed control values. These
#' express the MDD as a decrease; the classes follow Table 1 of Brock et al.
#' (2015): 0 above 100 % (no effect can be determined), I from 90 to 100 %,
#' II from 70 to 90 %, III from 50 to 70 % and IV below 50 % (small effects
#' can be determined). A value on a boundary goes to the class above it
#' (for example 70 % is class II), except 100 %, which is class I.
#'
#' This MDD has no power term: a true difference equal to it is detected
#' about half of the time. Mair et al. (2020) recommend judging
#' non-significant results with the confidence interval of the effect, which
#' is why it is returned alongside.
#'
#' @references
#' Brock TCM, Hammers-Wirtz M, Hommen U, Preuss TG, Ratte H-T, Roessink I,
#' Strauss T, Van den Brink PJ (2015). The minimum detectable difference (MDD)
#' and the interpretation of treatment-related effects of pesticides in
#' experimental ecosystems. *Environmental Science and Pollution Research* 22,
#' 1160–1174. \doi{10.1007/s11356-014-3398-2}
#'
#' Mair MM, Kattwinkel M, Jakoby O, Hartig F (2020). The minimum detectable
#' difference (MDD) concept for establishing trust in nonsignificant results:
#' a critical review. *Environmental Toxicology and Chemistry* 39, 2109–2123.
#' \doi{10.1002/etc.4847}
#' @seealso [test_control()], [plot_mdd()].
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" &
#'                   habitat == "eutrophic_lake" & taxon %in% c("Chironomini", "Cloeon"))
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' d <- transform_response(d, A = 2)
#' mdd <- calc_mdd(test_control(d, method = "dunnett"))
#' mdd[mdd$concentration == 200, c("taxon", "time", "mdd_pct", "mdd_class")]
calc_mdd <- function(x) {
  check_meso_test(x)
  tests <- x$tests
  # Williams with two directions: both rows of a dose share the scale.
  first_dir <- tests$direction == tests$direction[1L]
  tab <- tests[first_dir, , drop = FALSE]
  cells <- x$cells
  status <- cells$status[match(paste(tab$taxon, format_time_key(tab$time), sep = "\u001f"),
                               paste(cells$taxon, format_time_key(cells$time), sep = "\u001f"))]
  tr <- x$settings$transform
  mdd <- tab$critical_value * tab$se
  if (identical(tr$method, "log")) {
    back <- function(y) (exp(y) - 1) / tr$A
    control_value <- back(tab$mean_control)
    mdd_value <- (exp(tab$mean_control) - exp(tab$mean_control - mdd)) / tr$A
    effect <- back(tab$mean_treatment) - control_value
  } else {
    control_value <- tab$mean_control_value
    mdd_value <- mdd
    effect <- tab$mean_treatment_value - control_value
  }
  ok <- status == "ok" & !is.na(control_value) & control_value > 0
  pct <- function(v) ifelse(ok, 100 * v / control_value, NA_real_)
  out <- data.frame(
    taxon = tab$taxon, time = tab$time, treatment = tab$treatment,
    concentration = tab$concentration, status = status, n_control = tab$n_control,
    n_treatment = tab$n_treatment, control_value = ifelse(status == "ok", control_value, NA_real_),
    difference = tab$difference, lower = tab$lower, upper = tab$upper, ci_type = tab$ci_type,
    effect_pct = pct(effect), se = tab$se, df = tab$df, critical_value = tab$critical_value,
    mdd = mdd, mdd_value = ifelse(status == "ok", mdd_value, NA_real_), mdd_pct = pct(mdd_value),
    stringsAsFactors = FALSE
  )
  out$mdd_class <- mdd_class(out$mdd_pct)
  rownames(out) <- NULL
  out
}

# Classes of Table 1 of Brock et al. (2015).
mdd_class <- function(pct) {
  out <- rep(NA_character_, length(pct))
  out[!is.na(pct) & pct < 50] <- "IV"
  out[!is.na(pct) & pct >= 50 & pct < 70] <- "III"
  out[!is.na(pct) & pct >= 70 & pct < 90] <- "II"
  out[!is.na(pct) & pct >= 90 & pct <= 100] <- "I"
  out[!is.na(pct) & pct > 100] <- "0"
  out
}
