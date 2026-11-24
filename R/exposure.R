#' Measured exposure in a mesocosm experiment
#'
#' `meso_exposure()` validates a table of concentrations of the test substance
#' measured in the mesocosms, with one row per mesocosm and sampling occasion,
#' and assigns every sample to the application that preceded it. Its result is
#' the input of [fit_dissipation()], [calc_twa()], [compare_nominal()] and
#' [plot_exposure()].
#'
#' @param data A data frame with one row per mesocosm and sampling occasion.
#'   Tibbles are accepted; grouped data frames are not.
#' @param control The level of the treatment column that is the control
#'   (character; numbers are compared as text, so `"0"` for a column of
#'   concentrations). Controls are kept but left out of fits, time-weighted
#'   averages and comparisons with the nominal concentration.
#' @param applications The applications of the substance (no default). Either
#'   a vector with their times, on the scale of `time_col`, when all the
#'   treatments were applied at the same times (a single number for a single
#'   application); or, when each treatment has its own schedule or doses, a
#'   data frame with one row per treatment and application and the columns
#'   `treatment` (a level of the treatment column), `time` and, optionally,
#'   `nominal` (the nominal concentration that the application adds, in
#'   `unit`). Without `nominal`, each application adds the nominal
#'   concentration of its treatment (`concentration_col`). A treatment without
#'   rows in the data frame (usually the control) has no applications.
#' @param mesocosm_col,time_col,treatment_col,value_col Names of the columns
#'   with the mesocosm, the sampling time (numeric, for example days relative
#'   to the first application), the treatment and the measured concentration
#'   (numeric, not negative; `NA` when the sample has no value).
#' @param concentration_col Name of the column with the nominal concentration
#'   of each mesocosm (numeric, zero in the control), or `NULL` (default). It
#'   orders the treatments, labels the plots and is the reference of the
#'   time-weighted averages as a percentage of the nominal concentration.
#' @param timing_col Name of a column that says whether the samples taken at
#'   an application time were taken `"before"` or `"after"` the application
#'   (`NA` for samples taken at other times), or `NULL` (default). Without it,
#'   samples taken at an application time count as taken after it, and a
#'   mesocosm cannot have two samples at the same time.
#' @param unit Unit of the measured and nominal concentrations (for example
#'   `"µg/L"`), or `NULL`.
#' @param time_unit Unit of time, used to label rate constants and
#'   half-lives (default `"d"`).
#'
#' @return An object of class `meso_exposure`: a list with
#'   * `data`: one row per sample, ordered by mesocosm, time and timing
#'     (before, then after), with columns `mesocosm` (character), `treatment`
#'     (character), `concentration` (nominal; `NA` without
#'     `concentration_col`), `time`, `timing` (`"before"`, `"after"` or `NA`),
#'     `period` (integer: number of the application of its treatment that the
#'     sample follows, 0 before the first one), `application` (time of that
#'     application, `NA` in period 0), `elapsed` (time since that
#'     application), `value` (in `unit`) and `value_status` (`"measured"` for
#'     a positive value, `"zero"` or `"missing"`);
#'   * `applications`: one row per treatment and application, with
#'     `treatment`, `application` (1, 2, ...), `time` and `nominal` (the
#'     nominal concentration that it adds; `NA` when unknown);
#'   * `design`: one row per mesocosm with its treatment and concentration;
#'   * `settings`: the control, the order of the treatments, whether all the
#'     treatments share the schedule (`shared_schedule`), the units, whether
#'     the timing at the applications was assumed (`timing_assumed`) and the
#'     original column names;
#'   * `diagnostics`: counts of measured, zero and missing values, and of
#'     samples of treated mesocosms taken before their first application.
#'
#' @section Missing values and zeros:
#' A missing value (`NA`) is a sample without a usable measurement: not
#' analysed, lost, or reported below a limit without a number. It is kept,
#' marked `"missing"` and left out of [fit_dissipation()], [calc_twa()] and
#' [compare_nominal()]; it is never replaced. A zero is kept as given. Values
#' below a limit of quantification are not substituted by the package: if a
#' number is wanted for them (for example half the limit), put it in `data`.
#'
#' @section Assumptions:
#' The mesocosm is the experimental unit, with a single treatment for the
#' whole experiment, and every application of a treatment reaches all its
#' mesocosms at the same time.
#'
#' @seealso [fit_dissipation()], [calc_twa()], [compare_nominal()],
#'   [plot_exposure()], [azoxy_exposure].
#' @export
#' @examples
#' # Two applications, on days 0 and 14, in all the treatments
#' expo <- subset(azoxy_exposure, habitat == "eutrophic_lake")
#' e <- meso_exposure(expo, control = "0", applications = c(0, 14),
#'                    mesocosm_col = "mesocosm", time_col = "day",
#'                    treatment_col = "concentration",
#'                    concentration_col = "concentration",
#'                    value_col = "azoxystrobin", timing_col = "timing",
#'                    unit = "µg/L")
#' e
#' head(e$data)
#'
#' # A schedule for each treatment: one pulse of 20 µg/L against four weekly
#' # pulses of 5 µg/L
#' schedule <- data.frame(treatment = c("single", rep("weekly", 4)),
#'                        time = c(0, 0, 7, 14, 21), nominal = c(20, 5, 5, 5, 5))
#' pulses <- data.frame(
#'   mesocosm = rep(c("C1", "S1", "S2", "W1", "W2"), each = 4),
#'   treatment = rep(c("control", "single", "single", "weekly", "weekly"), each = 4),
#'   time = rep(c(0, 3, 10, 24), 5),
#'   value = c(rep(NA, 4), 19.1, 14.2, 7.6, 1.9, 20.6, 15.1, 7.1, 1.7,
#'             4.8, 3.6, 5.7, 7.1, 5.1, 3.8, 5.3, 6.7)
#' )
#' p <- meso_exposure(pulses, control = "control", applications = schedule,
#'                    unit = "µg/L")
#' p
#' p$applications
meso_exposure <- function(data, control, applications, mesocosm_col = "mesocosm",
                          time_col = "time", treatment_col = "treatment",
                          value_col = "value", concentration_col = NULL,
                          timing_col = NULL, unit = NULL, time_unit = "d") {
  if (missing(control)) {
    abort_input(paste(
      "`control` is missing: give the level of the treatment column that is",
      "the control, for example control = \"0\" or control = \"C\"."
    ))
  }
  if (missing(applications)) {
    abort_input(paste(
      "`applications` is missing: give the times of the applications of the",
      "substance, for example applications = c(0, 14), or a data frame with a",
      "schedule for each treatment."
    ))
  }
  check_data_frame(data, allow_empty = FALSE)
  check_string(mesocosm_col)
  check_string(time_col)
  check_string(treatment_col)
  check_string(value_col)
  check_string(concentration_col, allow_null = TRUE)
  check_string(timing_col, allow_null = TRUE)
  check_string(unit, allow_null = TRUE)
  check_string(time_unit)
  cols <- c(mesocosm_col = mesocosm_col, time_col = time_col,
            treatment_col = treatment_col, value_col = value_col,
            concentration_col = concentration_col, timing_col = timing_col)
  check_columns(data, cols)

  mesocosm <- key_text(data[[mesocosm_col]], mesocosm_col)
  treatment <- key_text(data[[treatment_col]], treatment_col)
  if (inherits(data[[time_col]], c("Date", "POSIXt"))) {
    abort_input(sprintf(
      paste(
        "Column '%s' holds dates. Give the sampling time as a number, for",
        "example the days relative to the first application:",
        "as.numeric(date - first_application)."
      ),
      time_col
    ))
  }
  check_numeric_column(data, time_col)
  check_numeric_column(data, value_col, lower = 0, allow_na = TRUE)
  concentration <- rep(NA_real_, nrow(data))
  if (!is.null(concentration_col)) {
    check_numeric_column(data, concentration_col, lower = 0)
    concentration <- as.numeric(data[[concentration_col]])
  }
  if (is.numeric(control) && length(control) == 1L && is.finite(control)) {
    control <- as.character(control)
  }
  check_string(control)
  time <- as.numeric(data[[time_col]])
  obs <- data.frame(mesocosm = mesocosm, treatment = treatment,
                    concentration = concentration, time = time,
                    value = as.numeric(data[[value_col]]), stringsAsFactors = FALSE)
  design <- build_design(obs, control, has_concentration = !is.null(concentration_col),
                         cols = cols)
  treatments <- design_levels(design, control)
  schedule <- exposure_schedule(applications, design, treatments, control, cols)
  shared <- !is.data.frame(applications)
  app_times <- split(schedule$time, factor(schedule$treatment, levels = treatments))

  # Samples taken at an application time of their treatment: before or after it.
  at_app <- vapply(seq_along(time), function(i) time[i] %in% app_times[[treatment[i]]],
                   logical(1))
  treated <- treatment != control
  timing <- rep(NA_character_, nrow(data))
  timing_assumed <- FALSE
  if (!is.null(timing_col)) {
    tm <- data[[timing_col]]
    if (is.factor(tm)) tm <- as.character(tm)
    if (!is.character(tm) && !all(is.na(tm))) {
      abort_input(sprintf("Column '%s' must hold \"before\", \"after\" or NA, not %s.",
                          timing_col, class(tm)[1L]))
    }
    tm <- as.character(tm)
    bad <- which(!is.na(tm) & !tm %in% c("before", "after"))
    if (length(bad) > 0L) {
      abort_input(sprintf(
        "Column '%s' must hold \"before\", \"after\" or NA; it has %s in %s.",
        timing_col, format_items(quote_items(unique(tm[bad]))), format_rows(bad)
      ))
    }
    no_timing <- which(at_app & is.na(tm))
    if (length(no_timing) > 0L) {
      abort_input(sprintf(
        paste(
          "Samples taken at an application time need \"before\" or \"after\" in",
          "column '%s', which is NA in %s."
        ),
        timing_col, format_rows(no_timing)
      ))
    }
    any_app <- time %in% schedule$time
    off_app <- which(treated & !any_app & !is.na(tm))
    if (length(off_app) > 0L) {
      abort_input(sprintf(
        paste(
          "Column '%s' says \"before\" or \"after\" for samples of treated",
          "mesocosms that were not taken at an application time (%s; times %s).",
          "Check `applications`, or set the column to NA for those samples."
        ),
        timing_col, format_rows(off_app),
        format_items(sort(unique(time[off_app])))
      ))
    }
    # Samples taken at the application of another treatment (or in the
    # controls) keep their labels: they only tell apart samples taken at the
    # same time.
    keep <- at_app | any_app | !treated
    timing[keep] <- tm[keep]
  } else if (any(at_app)) {
    timing[at_app] <- "after"
    timing_assumed <- TRUE
  }
  obs$timing <- timing

  key <- paste(obs$mesocosm, format_time_key(obs$time), obs$timing, sep = "\u001f")
  dup <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
  if (length(dup) > 0L) {
    hint <- if (is.null(timing_col) && any(at_app[dup])) {
      paste(" Samples taken before and after an application at the same time",
            "need a column that says which is which (timing_col).")
    } else {
      ""
    }
    abort_input(sprintf(
      "Each mesocosm can have one sample per time (and timing); these rows repeat one: %s.%s",
      format_rows(dup), hint
    ))
  }

  # Application period: number of applications of the treatment at or before
  # the time; a sample taken right before an application belongs to the
  # previous one.
  period <- integer(nrow(obs))
  application <- rep(NA_real_, nrow(obs))
  for (trt in unique(obs$treatment)) {
    i <- which(obs$treatment == trt)
    tt <- app_times[[trt]]
    if (length(tt) == 0L) next
    p <- findInterval(obs$time[i], tt)
    before <- obs$timing[i] %in% "before" & at_app[i]
    p[before] <- p[before] - 1L
    period[i] <- p
    application[i[p > 0L]] <- tt[p[p > 0L]]
  }
  status <- ifelse(is.na(obs$value), "missing", ifelse(obs$value == 0, "zero", "measured"))
  out <- data.frame(
    mesocosm = obs$mesocosm, treatment = obs$treatment,
    concentration = obs$concentration, time = obs$time, timing = obs$timing,
    period = as.integer(period), application = application,
    elapsed = obs$time - application, value = obs$value, value_status = status,
    stringsAsFactors = FALSE
  )
  ord <- order(out$mesocosm, out$time, match(out$timing, c("before", "after"), nomatch = 0L),
               method = "radix")
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  settings <- list(
    control = control, treatments = treatments,
    has_concentration = !is.null(concentration_col), shared_schedule = shared,
    unit = unit, time_unit = time_unit, timing_assumed = timing_assumed,
    columns = as.list(cols)
  )
  treated_out <- out$treatment != control
  diagnostics <- list(
    n_measured = sum(out$value_status == "measured"),
    n_zero = sum(out$value_status == "zero"),
    n_missing = sum(out$value_status == "missing"),
    n_before_first = sum(treated_out & out$period == 0L)
  )
  structure(list(data = out, applications = schedule, design = design,
                 settings = settings, diagnostics = diagnostics),
            class = "meso_exposure")
}

# Schedule of applications: one row per treatment and application, with its
# time and the nominal concentration that it adds.
exposure_schedule <- function(applications, design, treatments, control, cols,
                              call = sys.call(-1)) {
  conc <- design$concentration[match(treatments, design$treatment)]
  if (!is.data.frame(applications)) {
    if (!is.numeric(applications) || length(applications) == 0L ||
        anyNA(applications) || any(!is.finite(applications))) {
      abort_input(paste(
        "`applications` must be a vector of finite numbers (the times of the",
        "applications) or a data frame with columns treatment, time and,",
        "optionally, nominal."
      ), call = call)
    }
    tt <- sort(unique(as.numeric(applications)))
    out <- data.frame(treatment = rep(treatments, each = length(tt)),
                      application = rep(seq_along(tt), length(treatments)),
                      time = rep(tt, length(treatments)),
                      nominal = rep(conc, each = length(tt)), stringsAsFactors = FALSE)
    return(out)
  }
  check_data_frame(applications, allow_empty = FALSE, call = call)
  extra <- setdiff(names(applications), c("treatment", "time", "nominal"))
  if (length(extra) > 0L) {
    abort_input(sprintf(
      paste(
        "`applications` has columns that are not used: %s. Its columns are",
        "'treatment', 'time' and, optionally, 'nominal'."
      ),
      format_items(quote_items(extra))
    ), call = call)
  }
  check_columns(applications, c("treatment", "time"), arg = "applications", call = call)
  trt <- key_text(applications$treatment, "treatment", call = call)
  check_numeric_column(applications, "time", arg = "applications", call = call)
  unknown <- setdiff(trt, treatments)
  if (length(unknown) > 0L) {
    abort_input(sprintf(
      "`applications` has treatments that are not in column '%s': %s. Its levels are %s.",
      cols[["treatment_col"]], format_items(quote_items(unknown)),
      format_items(quote_items(treatments), max = 20L)
    ), call = call)
  }
  nominal <- conc[match(trt, treatments)]
  if ("nominal" %in% names(applications)) {
    check_numeric_column(applications, "nominal", lower = 0, allow_na = TRUE,
                         arg = "applications", call = call)
    nominal <- as.numeric(applications$nominal)
  }
  out <- data.frame(treatment = trt, time = as.numeric(applications$time),
                    nominal = nominal, stringsAsFactors = FALSE)
  key <- paste(out$treatment, format_time_key(out$time), sep = "\u001f")
  dup <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
  if (length(dup) > 0L) {
    abort_input(sprintf(
      "`applications` repeats an application (same treatment and time) in %s.",
      format_rows(dup)
    ), call = call)
  }
  if (any(out$treatment == control & !is.na(out$nominal) & out$nominal > 0)) {
    abort_input(sprintf(
      paste(
        "The control (%s) cannot receive the substance: its nominal",
        "concentration in `applications` must be 0."
      ),
      quote_items(control)
    ), call = call)
  }
  none <- setdiff(treatments[-1L], out$treatment)
  if (length(none) > 0L) {
    abort_input(sprintf(
      paste(
        "Treatments without applications in `applications`: %s. Give at least",
        "one for each treatment."
      ),
      format_items(quote_items(none))
    ), call = call)
  }
  out <- out[order(match(out$treatment, treatments), out$time, method = "radix"), ,
             drop = FALSE]
  out$application <- as.integer(stats::ave(out$time, out$treatment, FUN = seq_along))
  out <- out[c("treatment", "application", "time", "nominal")]
  rownames(out) <- NULL
  out
}

# Application times of a treatment, in order.
schedule_times <- function(schedule, treatment) {
  schedule$time[schedule$treatment == treatment]
}

check_meso_exposure <- function(x, arg = arg_name(substitute(x)), call = sys.call(-1)) {
  if (!inherits(x, "meso_exposure")) {
    abort_input(sprintf("`%s` must be the result of meso_exposure(), not %s.",
                        arg, class(x)[1L]), call = call)
  }
  invisible(x)
}

#' @export
print.meso_exposure <- function(x, ...) {
  s <- x$settings
  d <- x$data
  g <- x$diagnostics
  sch <- x$applications
  cat("<meso_exposure>\n")
  cat(sprintf("%d mesocosms, %d samples: %d measured, %d zero, %d missing\n",
              nrow(x$design), nrow(d), g$n_measured, g$n_zero, g$n_missing))
  conc <- x$design$concentration[match(s$treatments, x$design$treatment)]
  trt <- s$treatments
  if (s$has_concentration && any(trt != format_conc(conc))) {
    trt <- paste0(trt, " (", format_conc(conc), ")")
  }
  trt[1L] <- paste(trt[1L], "[control]")
  cat("Treatments: ", paste(trt, collapse = ", "),
      if (is.null(s$unit)) "" else paste0("; unit ", s$unit), "\n", sep = "")
  times <- sort(unique(d$time))
  both <- vapply(times, function(t) {
    all(c("before", "after") %in% d$timing[d$time == t])
  }, logical(1))
  lab <- format(times, trim = TRUE, drop0trailing = TRUE)
  lab[both] <- paste(lab[both], "(before and after)")
  cat("Sampling times: ", paste(lab, collapse = ", "), "\n", sep = "")
  if (isTRUE(s$shared_schedule)) {
    tt <- unique(sch$time)
    cat("Applications at time ", paste(format(tt, trim = TRUE), collapse = ", "),
        if (s$has_concentration) ", each at the nominal concentration of the treatment" else "",
        "\n", sep = "")
  } else {
    cat("Applications per treatment, time (nominal):\n")
    for (t in unique(sch$treatment)) {
      z <- sch[sch$treatment == t, , drop = FALSE]
      cat("  ", t, ": ", paste0(format(z$time, trim = TRUE),
                                ifelse(is.na(z$nominal), "",
                                       paste0(" (", format_conc(z$nominal), ")")),
                                collapse = ", "), "\n", sep = "")
    }
  }
  if (isTRUE(s$timing_assumed)) {
    cat("Samples taken at an application time count as taken after it (no timing column)\n")
  }
  if (g$n_before_first > 0L) {
    cat(sprintf("%d samples of treated mesocosms before their first application (not used)\n",
                g$n_before_first))
  }
  invisible(x)
}

#' @export
as.data.frame.meso_exposure <- function(x, ...) {
  x$data
}

#' Dissipation of the test substance after each application
#'
#' `fit_dissipation()` fits first-order (SFO, single first-order) dissipation
#' to the concentrations measured after each application in a
#' [meso_exposure()] object, and derives the rate constant, the half-life
#' (DT50) and the DT90.
#'
#' @param x A `meso_exposure` object.
#' @param method `"log_linear"` (default): least squares regression of ln C on
#'   the time since the application; `"nonlinear"`: nonlinear least squares
#'   of C itself ([stats::nls()]), started from the log-linear fit.
#' @param by `"treatment"` (default): one fit per treatment, with the values of
#'   all its mesocosms; `"mesocosm"`: one fit per mesocosm.
#' @param k_by `"period"` (default): a rate constant for each application,
#'   from the values of its period alone; `"series"`: one rate constant for all
#'   the applications, with an initial concentration for each.
#' @param periods The application periods to fit (whole numbers: 1 for the
#'   period after the first application of each treatment, 2 after the second,
#'   ...), or `NULL` (default) for all of them. A treatment with fewer
#'   applications has no rows for the periods it lacks.
#' @param level Confidence level of the intervals (default 0.95).
#'
#' @details
#' After application \eqn{p}, made at time \eqn{t_p}, the concentration
#' follows
#' \deqn{C(t) = C_{0,p}\, e^{-k (t - t_p)}, \qquad t_p \le t < t_{p+1},}{
#'   C(t) = C0p exp(-k (t - tp)), tp <= t < tp+1,}
#' where \eqn{C_{0,p}} is the concentration right after the application and
#' \eqn{k} the first-order rate constant; then \eqn{DT_{50} = \ln 2 / k}{DT50
#' = ln 2 / k} and \eqn{DT_{90} = \ln 10 / k}{DT90 = ln 10 / k} (FOCUS 2014,
#' section 4.2.1). The fits use the treated mesocosms and the samples with a
#' positive value; controls, zeros and missing values are left out (see
#' `diagnostics$excluded`).
#'
#' The log-linear regression fits \eqn{\ln C = a_p - k (t - t_p)}{ln C = ap -
#' k (t - tp)} by least squares, so \eqn{C_{0,p} = e^{a_p}}{C0p = exp(ap)} is
#' on the scale of the geometric mean. The logarithm gives the same relative
#' weight to all concentrations; FOCUS (2014, section 6.1.3) recommends fits
#' of the untransformed concentrations (`method = "nonlinear"`), because the
#' transformation changes the weight of each value and can lengthen the DT50
#' when the decline is fast at first and slower later.
#'
#' With `k_by = "period"` each period is fitted separately, with its own
#' residual variance; with `k_by = "series"` all the periods of a treatment
#' (or mesocosm) are fitted together. Each treatment follows its own schedule
#' of applications (see [meso_exposure()]).
#'
#' `c0_pct_nominal` compares \eqn{C_{0,p}}{C0p} with the nominal concentration
#' that the application adds: after the first application it measures the
#' dosing; after later ones it also includes what remained of the earlier
#' applications, so it can exceed 100 %.
#'
#' @section Uncertainty:
#' The intervals of \eqn{k} are Wald confidence intervals with the
#' \eqn{t}{t} distribution and the residual degrees of freedom of the fit
#' (`df`); those of the DT50 invert them. They include only the uncertainty
#' of the parameters, given the model. When the lower limit of \eqn{k} is not
#' positive, the data are compatible with no dissipation and the upper limit
#' of the DT50 is `Inf`.
#'
#' @section Estimates not available:
#' Estimates that cannot be computed are `NA`, with the cause in `status` and
#' `reason`: `"not_estimable"` (fewer than two sampling times with a positive
#' value, or a nonlinear fit that did not converge), `"no_decline"` (the
#' fitted \eqn{k} is not positive, so there is no DT50) and
#' `"no_uncertainty"` (no residual degrees of freedom, so no intervals).
#' Otherwise `status` is `"ok"`.
#'
#' @return An object of class `meso_dissipation`: a list with
#'   * `estimates`: one row per treatment (or mesocosm) and period, with
#'     `treatment`, `concentration`, `mesocosm` (with `by = "mesocosm"`),
#'     `period`, `application` (its time), `n` (values used), `n_times`
#'     (distinct sampling times), `df` (residual degrees of freedom of the
#'     fit), `c0` (fitted concentration right after the application, in the
#'     unit of `x`), `nominal` (the nominal concentration that the application
#'     adds, from `x$applications`), `c0_pct_nominal` (100 `c0` / `nominal`;
#'     `NA` when the nominal concentration is unknown or zero), `k`, `k_se`,
#'     `k_lower` and `k_upper` (per `time_unit`),
#'     `dt50`, `dt50_lower`, `dt50_upper` and `dt90` (in `time_unit`),
#'     `status` and `reason`. With `k_by = "series"`, `k`, its intervals,
#'     the DT50 and the DT90 repeat in the rows of each treatment (or
#'     mesocosm);
#'   * `fitted`: the values used, with their fitted values;
#'   * `settings`: the options and what [predict()], [calc_twa()] and
#'     [plot_exposure()] need from `x`;
#'   * `diagnostics`: the samples left out (`excluded`) and their number.
#'
#' @references
#' FOCUS (2014). *Generic guidance for estimating persistence and degradation
#' kinetics from environmental fate studies on pesticides in EU registration*,
#' version 1.1 (18 December 2014). Update of the report of the FOCUS Work
#' Group on Degradation Kinetics, EC document reference Sanco/10058/2005
#' version 2.0 (2006).
#'
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @seealso [meso_exposure()], [predict.meso_dissipation()], [calc_twa()],
#'   [plot_exposure()].
#' @export
#' @examples
#' expo <- subset(azoxy_exposure, habitat == "eutrophic_lake")
#' e <- meso_exposure(expo, control = "0", applications = c(0, 14),
#'                    mesocosm_col = "mesocosm", time_col = "day",
#'                    treatment_col = "concentration",
#'                    concentration_col = "concentration",
#'                    value_col = "azoxystrobin", timing_col = "timing",
#'                    unit = "µg/L")
#' dis <- fit_dissipation(e)
#' dis
#' dis$estimates[c("treatment", "period", "k", "dt50", "dt50_lower", "dt50_upper")]
#'
#' # Untransformed concentrations, as recommended by FOCUS (2014)
#' fit_dissipation(e, method = "nonlinear")
fit_dissipation <- function(x, method = "log_linear", by = "treatment", k_by = "period",
                            periods = NULL, level = 0.95) {
  check_meso_exposure(x)
  check_choice(method, c("log_linear", "nonlinear"))
  check_choice(by, c("treatment", "mesocosm"))
  check_choice(k_by, c("period", "series"))
  check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  s <- x$settings
  sch <- x$applications
  m <- max(table(sch$treatment[sch$treatment != s$control]))
  if (is.null(periods)) {
    periods <- seq_len(m)
  } else {
    if (!is.numeric(periods) || length(periods) == 0L || anyNA(periods) ||
        any(!is.finite(periods)) || any(periods != round(periods)) ||
        any(periods < 1 | periods > m)) {
      abort_input(sprintf(
        paste(
          "`periods` must be NULL or whole numbers from 1 to %d (the most",
          "applications of a treatment), not %s."
        ),
        m, describe_value(periods)
      ))
    }
    periods <- sort(unique(as.integer(periods)))
  }

  d <- x$data
  d <- d[d$treatment != s$control & d$period %in% periods, , drop = FALSE]
  used <- d$value_status == "measured"
  excluded <- d[!used, c("treatment", "mesocosm", "time", "timing", "period", "value_status"),
                drop = FALSE]
  rownames(excluded) <- NULL
  d <- d[used, , drop = FALSE]

  groups <- exposure_groups(x, by)
  gkey <- if (by == "treatment") "treatment" else "mesocosm"
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    dg <- d[d[[gkey]] == groups[[gkey]][i], , drop = FALSE]
    apps <- sch[sch$treatment == groups$treatment[i], , drop = FALSE]
    p_g <- periods[periods <= nrow(apps)]
    dg <- dg[dg$period %in% p_g, , drop = FALSE]
    if (k_by == "period") {
      for (p in p_g) {
        dp <- dg[dg$period == p, , drop = FALSE]
        f <- sfo_fit(dp$value, dp$elapsed, rep(1L, nrow(dp)), method, level)
        rows[[length(rows) + 1L]] <- sfo_rows(groups[i, , drop = FALSE], p, apps$time[p],
                                              apps$nominal[p], f, 1L)
      }
    } else {
      present <- p_g[p_g %in% dg$period]
      f <- sfo_fit(dg$value, dg$elapsed, match(dg$period, present), method, level)
      for (p in p_g) {
        j <- match(p, present)
        if (is.na(j)) {
          f_p <- sfo_fit(numeric(), numeric(), integer(), method, level)
          j <- 1L
        } else {
          f_p <- f
        }
        rows[[length(rows) + 1L]] <- sfo_rows(groups[i, , drop = FALSE], p, apps$time[p],
                                              apps$nominal[p], f_p, j)
      }
    }
  }
  est <- do.call(rbind, rows)
  rownames(est) <- NULL

  # Values used, with the fitted curve of their treatment (or mesocosm) and period
  fitted <- d[c("treatment", "concentration", "mesocosm", "period", "time", "timing",
                "elapsed", "value")]
  pos <- match(paste(fitted[[gkey]], fitted$period, sep = "\u001f"),
               paste(est[[gkey]], est$period, sep = "\u001f"))
  fitted$fitted <- est$c0[pos] * exp(-est$k[pos] * fitted$elapsed)
  rownames(fitted) <- NULL

  settings <- list(
    method = method, by = by, k_by = k_by, periods = periods, level = level,
    interval = "confidence (Wald, t distribution)", uncertainty = "parameter",
    control = s$control, treatments = s$treatments, design = x$design,
    has_concentration = s$has_concentration, applications = sch,
    unit = s$unit, time_unit = s$time_unit
  )
  diagnostics <- list(
    excluded = excluded,
    n_zero = sum(excluded$value_status == "zero"),
    n_missing = sum(excluded$value_status == "missing")
  )
  structure(list(estimates = est, fitted = fitted, settings = settings,
                 diagnostics = diagnostics),
            class = "meso_dissipation")
}

# Treated groups of a meso_exposure object: treatments, or mesocosms, in the
# order of the treatments.
exposure_groups <- function(x, by) {
  s <- x$settings
  des <- x$design[x$design$treatment != s$control, , drop = FALSE]
  des <- des[order(match(des$treatment, s$treatments), des$mesocosm, method = "radix"), ,
             drop = FALSE]
  if (by == "treatment") {
    out <- unique(des[c("treatment", "concentration")])
  } else {
    out <- des[c("treatment", "concentration", "mesocosm")]
  }
  rownames(out) <- NULL
  out
}

# First-order fit of positive values `y` at times `e` since the application,
# with an initial concentration for each period index in `g` (1, 2, ...) and
# a common rate constant. Returns the estimates and their status.
sfo_fit <- function(y, e, g, method, level) {
  n_g <- if (length(g) == 0L) 1L else max(g)
  n_times <- vapply(seq_len(n_g), function(j) length(unique(e[g == j])), integer(1))
  out <- list(c0 = rep(NA_real_, n_g), k = NA_real_, k_se = NA_real_, k_lower = NA_real_,
              k_upper = NA_real_, df = NA_integer_, n = tabulate(g, nbins = n_g),
              n_times = n_times, status = "not_estimable", reason = NA_character_)
  if (length(y) == 0L) {
    out$reason <- "no positive values in this period"
    return(out)
  }
  if (all(n_times < 2L)) {
    out$reason <- "fewer than two sampling times with a positive value"
    return(out)
  }
  npar <- n_g + 1L
  df <- length(y) - npar
  X <- cbind((outer(g, seq_len(n_g), "==")) + 0, -e)
  ll <- stats::lm.fit(X, log(y))
  if (ll$rank < npar) {
    out$reason <- "the sampling times do not allow the rate constant to be estimated"
    return(out)
  }
  cf <- ll$coefficients
  c0 <- exp(unname(cf[seq_len(n_g)]))
  k <- unname(cf[npar])
  se <- NA_real_
  if (df > 0L) {
    piv <- ll$qr$pivot[seq_len(npar)]
    V <- matrix(NA_real_, npar, npar)
    V[piv, piv] <- chol2inv(ll$qr$qr[seq_len(npar), seq_len(npar), drop = FALSE])
    se <- sqrt(sum(ll$residuals^2) / df * V[npar, npar])
  }
  if (method == "nonlinear") {
    sc <- max(y)
    dat <- list(y = y / sc, e = e, g = g)
    fml <- if (n_g == 1L) y ~ c0 * exp(-k * e) else y ~ c0[g] * exp(-k * e)
    nl <- tryCatch(
      stats::nls(fml, data = dat, start = list(c0 = c0 / sc, k = k),
                 control = stats::nls.control(maxiter = 200L, tol = 1e-8, scaleOffset = 1)),
      error = function(cnd) cnd
    )
    if (inherits(nl, "error")) {
      out$reason <- paste0("the nonlinear fit did not converge (", conditionMessage(nl), ")")
      return(out)
    }
    cf <- stats::coef(nl)
    c0 <- unname(cf[seq_len(n_g)]) * sc
    k <- unname(cf[npar])
    se <- if (df > 0L) sqrt(stats::vcov(nl)[npar, npar]) else NA_real_
  }
  out$c0 <- c0
  out$k <- k
  out$df <- as.integer(df)
  if (df > 0L) {
    q <- stats::qt((1 + level) / 2, df)
    out$k_se <- se
    out$k_lower <- k - q * se
    out$k_upper <- k + q * se
  }
  no_df <- "no residual degrees of freedom, so no confidence interval"
  if (k <= 0) {
    out$status <- "no_decline"
    out$reason <- paste0("the fitted rate constant is not positive: no decline",
                         if (df == 0L) paste0("; ", no_df) else "")
  } else if (df == 0L) {
    out$status <- "no_uncertainty"
    out$reason <- no_df
  } else {
    out$status <- "ok"
  }
  out
}

# One row of fit_dissipation()$estimates: period `p` (index `j` of the fit).
sfo_rows <- function(group, p, application, nominal, f, j) {
  k <- f$k
  positive <- is.finite(k) && k > 0
  dt50_lower <- if (is.finite(f$k_upper) && f$k_upper > 0) log(2) / f$k_upper else NA_real_
  dt50_upper <- if (is.finite(f$k_lower) && is.finite(f$k_upper) && f$k_upper > 0) {
    if (f$k_lower > 0) log(2) / f$k_lower else Inf
  } else {
    NA_real_
  }
  data.frame(
    group, period = as.integer(p), application = application, n = f$n[j],
    n_times = f$n_times[j], df = f$df, c0 = f$c0[j], nominal = nominal,
    c0_pct_nominal = pct_of(f$c0[j], nominal), k = k, k_se = f$k_se,
    k_lower = f$k_lower, k_upper = f$k_upper,
    dt50 = if (positive) log(2) / k else NA_real_, dt50_lower = dt50_lower,
    dt50_upper = dt50_upper, dt90 = if (positive) log(10) / k else NA_real_,
    status = f$status, reason = f$reason, stringsAsFactors = FALSE, row.names = NULL
  )
}

# Labels of the treatments in tables and plots of the exposure: the nominal
# concentration, preceded by the name of the treatment when the name is not
# that concentration (for example "weekly (5 µg/L)" against "single (20
# µg/L)"); the names alone when the concentrations are unknown.
exposure_labels <- function(levels_trt, conc, control, unit, labels = NULL,
                            call = sys.call(-1)) {
  lab <- treatment_labels(levels_trt, conc, control, unit, labels, call = call)
  if (!is.null(labels)) {
    if (anyDuplicated(lab)) {
      abort_input("`labels` must be different for each treatment.", call = call)
    }
    return(lab)
  }
  if (!all(is.na(conc))) {
    named <- levels_trt != control & levels_trt != format_conc(conc)
    lab[named] <- paste0(levels_trt[named], " (", lab[named], ")")
  }
  lab
}

# 100 x / reference, NA when the reference is unknown or zero.
pct_of <- function(x, reference) {
  ifelse(is.finite(reference) & reference > 0, 100 * x / reference, NA_real_)
}

check_meso_dissipation <- function(x, arg = arg_name(substitute(x)), call = sys.call(-1)) {
  if (!inherits(x, "meso_dissipation")) {
    abort_input(sprintf("`%s` must be the result of fit_dissipation(), not %s.",
                        arg, class(x)[1L]), call = call)
  }
  invisible(x)
}

# The fit must come from the same experiment as `x`.
check_same_exposure <- function(fit, x, call = sys.call(-1)) {
  same <- identical(fit$settings$applications, x$applications) &&
    identical(fit$settings$treatments, x$settings$treatments) &&
    identical(fit$settings$design, x$design)
  if (!same) {
    abort_input(paste(
      "`fit` does not come from `x`: its mesocosms, treatments or applications",
      "differ. Fit it with fit_dissipation(x)."
    ), call = call)
  }
  invisible(fit)
}

#' @export
print.meso_dissipation <- function(x, digits = 3, ...) {
  s <- x$settings
  est <- x$estimates
  cat("<meso_dissipation>\n")
  cat(sprintf("First-order dissipation: %s; %s; %s\n",
              if (s$method == "log_linear") "log-linear regression of ln C" else
                "nonlinear least squares of C",
              if (s$k_by == "period") "one rate constant per application" else
                "one rate constant for all the applications",
              if (s$by == "treatment") "one fit per treatment" else "one fit per mesocosm"))
  cat(sprintf(paste0("%s%% confidence intervals (Wald, parameter uncertainty only; those of ",
                     "k in $estimates); k per %s, DT50 and DT90 in %s\n"),
              format(100 * s$level), s$time_unit, s$time_unit))
  fmt <- function(v) {
    ifelse(is.na(v), "NA", ifelse(is.infinite(v), "Inf",
                                  trimws(formatC(signif(v, digits), digits = digits,
                                                 format = "fg"))))
  }
  conc <- s$design$concentration[match(s$treatments, s$design$treatment)]
  lab <- exposure_labels(s$treatments, conc, s$control, s$unit)
  tab <- data.frame(treatment = unname(lab[est$treatment]), stringsAsFactors = FALSE)
  if (s$by == "mesocosm") tab$mesocosm <- est$mesocosm
  tab$period <- est$period
  tab$n <- est$n
  tab$c0 <- fmt(est$c0)
  if (any(!is.na(est$c0_pct_nominal))) {
    tab$c0 <- ifelse(is.na(est$c0_pct_nominal), tab$c0,
                     paste0(tab$c0, " (", round(est$c0_pct_nominal), "%)"))
    names(tab)[names(tab) == "c0"] <- "C0 (% nominal)"
  } else {
    names(tab)[names(tab) == "c0"] <- "C0"
  }
  tab$k <- fmt(est$k)
  tab$DT50 <- fmt(est$dt50)
  tab$dt50_ci <- ifelse(is.na(est$dt50_lower), "",
                        paste0("[", fmt(est$dt50_lower), ", ", fmt(est$dt50_upper), "]"))
  tab$DT90 <- fmt(est$dt90)
  tab$status <- est$status
  names(tab)[names(tab) == "dt50_ci"] <- "DT50 CI"
  print(tab, row.names = FALSE, right = FALSE)
  g <- x$diagnostics
  if (g$n_zero + g$n_missing > 0L) {
    cat(sprintf("Left out: %d zero and %d missing values of treated mesocosms\n",
                g$n_zero, g$n_missing))
  }
  invisible(x)
}

#' Fitted dissipation curves
#'
#' Concentrations predicted by a [fit_dissipation()] fit, on a grid of times
#' for each period or at the given times.
#'
#' @param object A `meso_dissipation` object.
#' @param times `NULL` (default), for `n` times from each application to the
#'   next one (to the last value used, in the last period); or a vector of
#'   times at which to predict. A time equal to an application time takes the
#'   curve of that application.
#' @param n Number of times per period when `times` is `NULL` (default 101).
#' @param ... Not used.
#'
#' @return A data frame with `treatment`, `concentration`, `mesocosm` (with
#'   `by = "mesocosm"`), `period`, `time`, `elapsed` (time since the
#'   application) and `fit` (predicted concentration). Periods whose fit is
#'   not available give no rows with `times = NULL` and `NA` otherwise, as do
#'   times before the first application of the treatment.
#' @seealso [fit_dissipation()].
#' @export
#' @examples
#' expo <- subset(azoxy_exposure, habitat == "eutrophic_lake")
#' e <- meso_exposure(expo, control = "0", applications = c(0, 14),
#'                    mesocosm_col = "mesocosm", time_col = "day",
#'                    treatment_col = "concentration",
#'                    concentration_col = "concentration",
#'                    value_col = "azoxystrobin", timing_col = "timing")
#' dis <- fit_dissipation(e)
#' predict(dis, times = c(0, 7, 14, 21, 56))
predict.meso_dissipation <- function(object, times = NULL, n = 101, ...) {
  check_meso_dissipation(object)
  check_number(n, lower = 2)
  if (n != round(n)) abort_input("`n` must be a whole number.")
  s <- object$settings
  est <- object$estimates
  gcols <- intersect(c("treatment", "concentration", "mesocosm"), names(est))
  gkey <- if (s$by == "treatment") "treatment" else "mesocosm"
  sch <- s$applications
  ok <- is.finite(est$c0) & is.finite(est$k)
  if (is.null(times)) {
    fit <- object$fitted
    rows <- lapply(which(ok), function(i) {
      p <- est$period[i]
      apps <- schedule_times(sch, est$treatment[i])
      if (p < length(apps)) {
        end <- apps[p + 1L]
      } else {
        in_p <- fit[[gkey]] == est[[gkey]][i] & fit$period == p
        end <- max(fit$time[in_p])
      }
      tt <- seq(apps[p], end, length.out = n)
      data.frame(est[i, gcols, drop = FALSE], period = p, time = tt, elapsed = tt - apps[p],
                 fit = est$c0[i] * exp(-est$k[i] * (tt - apps[p])),
                 stringsAsFactors = FALSE, row.names = NULL)
    })
  } else {
    if (!is.numeric(times) || length(times) == 0L || anyNA(times) || any(!is.finite(times))) {
      abort_input("`times` must be NULL or a vector of finite numbers.")
    }
    grp <- unique(est[gcols])
    rows <- lapply(seq_len(nrow(grp)), function(i) {
      apps <- schedule_times(sch, grp$treatment[i])
      p <- findInterval(times, apps)
      j <- match(paste(grp[[gkey]][i], p, sep = "\u001f"),
                 paste(est[[gkey]], est$period, sep = "\u001f"))
      j[!ok[j] %in% TRUE] <- NA_integer_
      a <- ifelse(p > 0L, apps[pmax(p, 1L)], NA_real_)
      data.frame(grp[rep(i, length(times)), , drop = FALSE], period = as.integer(p),
                 time = times, elapsed = times - a,
                 fit = est$c0[j] * exp(-est$k[j] * (times - a)),
                 stringsAsFactors = FALSE, row.names = NULL)
    })
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(est[0L, gcols, drop = FALSE], period = integer(), time = numeric(),
                      elapsed = numeric(), fit = numeric())
  }
  rownames(out) <- NULL
  out
}

#' Time-weighted average concentrations
#'
#' `calc_twa()` computes the time-weighted average (TWA) concentration of the
#' test substance over time windows that start at the first application of
#' each treatment (or at a given time), from the measured concentrations
#' (trapezoidal rule) or from the curves of [fit_dissipation()], also as a
#' percentage of the nominal concentration.
#'
#' @param x A `meso_exposure` object.
#' @param windows Lengths of the time windows, on the time scale of `x` (for
#'   example `c(7, 14, 21, 56)`; no default).
#' @param method `"trapezoid"` (default): linear interpolation between the
#'   measured concentrations; `"model"`: exact integral of the fitted curves of
#'   `fit`.
#' @param fit With `method = "model"`, the result of `fit_dissipation(x)`.
#' @param start Start of the windows: `NULL` (default) for the first
#'   application of each treatment, or a time, the same for all the
#'   treatments.
#' @param by `"treatment"` (TWA of the mean concentration of each treatment)
#'   or `"mesocosm"` (TWA of each mesocosm). Defaults to `"treatment"` with
#'   `method = "trapezoid"` and to the `by` of `fit` with `method = "model"`.
#'
#' @details
#' The TWA over the window \eqn{[a, b]}, with \eqn{b = a + w}, is
#' \deqn{TWA = \frac{1}{b - a} \int_a^b C(t)\, dt.}{
#'   TWA = (1 / (b - a)) integral from a to b of C(t) dt.}
#'
#' With `method = "trapezoid"`, \eqn{C(t)} is the linear interpolation between
#' consecutive samples of the same application period. With
#' `by = "treatment"`, each sample is the mean of the mesocosms of the
#' treatment with a value; a missing value is not replaced, and the
#' interpolation passes over a sample without any value. The profile jumps at
#' each application: the sample taken before it closes a period and the one
#' taken after it opens the next. A window that needs the concentration where
#' there are no samples to interpolate between (before the first sample,
#' after the last one, or right before or after an application without a
#' sample there) is not estimable.
#'
#' With `method = "model"`, \eqn{C(t) = C_{0,p} e^{-k (t - t_p)}}{C(t) = C0p
#' exp(-k (t - tp))} in each period, integrated exactly:
#' \deqn{\int_u^v C(t)\, dt = \frac{C_{0,p}}{k}
#'   \left(e^{-k (u - t_p)} - e^{-k (v - t_p)}\right).}{
#'   integral from u to v of C(t) dt = (C0p / k) (exp(-k (u - tp)) - exp(-k (v - tp))).}
#' Windows that end after the last measured value use the curve beyond the
#' data and are marked `"extrapolated"`.
#'
#' When the concentration declines exponentially, the straight line between
#' two samples lies above the curve, so the trapezoidal rule tends to give
#' larger TWAs than the model, more so the farther apart the samples are.
#'
#' `twa_pct_nominal` refers the TWA to the nominal concentration of the
#' treatment (`concentration_col` of [meso_exposure()]), whatever the number
#' of applications. Controls are left out.
#'
#' @return A data frame with one row per treatment (or mesocosm) and window:
#'   `treatment`, `concentration`, `mesocosm` (with `by = "mesocosm"`),
#'   `method`, `start`, `end`, `window`, `twa` (in `unit`), `twa_pct_nominal`
#'   (100 `twa` / `concentration`; `NA` when the nominal concentration is
#'   unknown or zero), `n_samples` (sampling occasions with a value in the
#'   window), `n_missing` (samples without a value in the window), `unit`,
#'   `status` (`"ok"`, `"extrapolated"` or `"not_estimable"`) and `reason`.
#' @references
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @seealso [meso_exposure()], [fit_dissipation()].
#' @export
#' @examples
#' expo <- subset(azoxy_exposure, habitat == "eutrophic_lake")
#' e <- meso_exposure(expo, control = "0", applications = c(0, 14),
#'                    mesocosm_col = "mesocosm", time_col = "day",
#'                    treatment_col = "concentration",
#'                    concentration_col = "concentration",
#'                    value_col = "azoxystrobin", timing_col = "timing",
#'                    unit = "µg/L")
#' calc_twa(e, windows = c(7, 14, 21, 56))
#' calc_twa(e, windows = c(7, 14, 21, 56), method = "model",
#'          fit = fit_dissipation(e))
calc_twa <- function(x, windows, method = "trapezoid", fit = NULL, start = NULL, by = NULL) {
  check_meso_exposure(x)
  if (missing(windows)) {
    abort_input(paste(
      "`windows` is missing: give the lengths of the time windows, for example",
      "windows = c(7, 14, 21)."
    ))
  }
  if (!is.numeric(windows) || length(windows) == 0L || anyNA(windows) ||
      any(!is.finite(windows)) || any(windows <= 0)) {
    abort_input("`windows` must be a vector of positive numbers (the lengths of the time windows).")
  }
  windows <- sort(unique(as.numeric(windows)))
  check_choice(method, c("trapezoid", "model"))
  s <- x$settings
  if (!is.null(start)) check_number(start)
  if (method == "model") {
    if (is.null(fit)) {
      abort_input(paste(
        "`fit` is missing: method = \"model\" integrates the curves of",
        "fit_dissipation(); give fit = fit_dissipation(x)."
      ))
    }
    check_meso_dissipation(fit)
    check_same_exposure(fit, x)
    if (is.null(by)) by <- fit$settings$by
    check_choice(by, c("treatment", "mesocosm"))
    if (by != fit$settings$by) {
      abort_input(sprintf(
        "`by` is \"%s\" but `fit` has one fit per %s; fit it with fit_dissipation(x, by = \"%s\").",
        by, fit$settings$by, by
      ))
    }
  } else {
    if (!is.null(fit)) abort_input("`fit` is only used with method = \"model\".")
    if (is.null(by)) by <- "treatment"
    check_choice(by, c("treatment", "mesocosm"))
  }

  groups <- exposure_groups(x, by)
  gkey <- if (by == "treatment") "treatment" else "mesocosm"
  d <- x$data[x$data$treatment != s$control, , drop = FALSE]
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    dg <- d[d[[gkey]] == groups[[gkey]][i], , drop = FALSE]
    apps <- schedule_times(x$applications, groups$treatment[i])
    occ <- twa_occasions(dg, apps)
    a <- if (is.null(start)) apps[1L] else start
    res <- lapply(windows, function(w) {
      b <- a + w
      r <- if (method == "trapezoid") {
        twa_trapezoid(occ, a, b)
      } else {
        est_g <- fit$estimates[fit$estimates[[gkey]] == groups[[gkey]][i], , drop = FALSE]
        twa_model(est_g, occ, apps, a, b)
      }
      inside <- occ$inside(a, b)
      data.frame(groups[i, , drop = FALSE], method = method, start = a, end = b, window = w,
                 twa = r$twa, twa_pct_nominal = pct_of(r$twa, groups$concentration[i]),
                 n_samples = sum(inside & !is.na(occ$profile$value)),
                 n_missing = sum(occ$profile$n_missing[inside]),
                 unit = if (is.null(s$unit)) NA_character_ else s$unit,
                 status = r$status, reason = r$reason,
                 stringsAsFactors = FALSE, row.names = NULL)
    })
    rows[[i]] <- do.call(rbind, res)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Sampling occasions of a treatment (or mesocosm): mean value of the
# mesocosms with a value and number of missing values, in time order.
twa_occasions <- function(dg, apps) {
  key <- paste(format_time_key(dg$time), dg$timing, sep = "\u001f")
  first <- !duplicated(key)
  prof <- dg[first, c("time", "timing", "period"), drop = FALSE]
  k <- key[first]
  prof$value <- vapply(k, function(z) {
    v <- dg$value[key == z]
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  }, numeric(1), USE.NAMES = FALSE)
  prof$n_missing <- vapply(k, function(z) sum(is.na(dg$value[key == z])), integer(1),
                           USE.NAMES = FALSE)
  prof <- prof[order(prof$time, match(prof$timing, c("before", "after"), nomatch = 0L),
                     method = "radix"), , drop = FALSE]
  rownames(prof) <- NULL
  m <- length(apps)
  lo <- ifelse(prof$period == 0L, -Inf, apps[pmax(prof$period, 1L)])
  hi <- ifelse(prof$period == m, Inf, apps[pmin(prof$period + 1L, m)])
  list(
    profile = prof,
    # Occasions in the window: in [a, b] and in a period that overlaps (a, b)
    inside = function(a, b) prof$time >= a & prof$time <= b & lo < b & hi > a
  )
}

twa_trapezoid <- function(occ, a, b) {
  prof <- occ$profile[!is.na(occ$profile$value), , drop = FALSE]
  n <- nrow(prof)
  area <- 0
  cover <- matrix(numeric(), ncol = 2L)
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      if (prof$period[i] != prof$period[i + 1L]) next
      t1 <- prof$time[i]
      t2 <- prof$time[i + 1L]
      u <- max(a, t1)
      v <- min(b, t2)
      if (v <= u) next
      c1 <- prof$value[i]
      c2 <- prof$value[i + 1L]
      cu <- c1 + (c2 - c1) * (u - t1) / (t2 - t1)
      cv <- c1 + (c2 - c1) * (v - t1) / (t2 - t1)
      area <- area + (cu + cv) / 2 * (v - u)
      cover <- rbind(cover, c(u, v))
    }
  }
  gap <- first_gap(cover, a, b)
  if (is.null(gap)) {
    return(list(twa = area / (b - a), status = "ok", reason = NA_character_))
  }
  reason <- if (n == 0L) {
    "no measured values"
  } else if (gap[2L] >= b && b > max(prof$time)) {
    sprintf("the window ends after the last measured value (time %s)", format(max(prof$time)))
  } else if (gap[1L] <= a && a < min(prof$time)) {
    sprintf("the window starts before the first measured value (time %s)", format(min(prof$time)))
  } else {
    sprintf(paste(
      "no measured values to interpolate between from time %s to %s (values are",
      "interpolated only between samples of the same application period)"
    ), format(gap[1L]), format(gap[2L]))
  }
  list(twa = NA_real_, status = "not_estimable", reason = reason)
}

# First part of [a, b] not covered by the intervals (rows of `cover`), or NULL.
first_gap <- function(cover, a, b) {
  eps <- 1e-9 * max(1, abs(a), abs(b))
  cover <- cover[order(cover[, 1L]), , drop = FALSE]
  cur <- a
  for (i in seq_len(nrow(cover))) {
    if (cover[i, 1L] > cur + eps) return(c(cur, cover[i, 1L]))
    cur <- max(cur, cover[i, 2L])
  }
  if (cur < b - eps) c(cur, b) else NULL
}

# Integral of c0 exp(-k s) from s1 to s2.
sfo_integral <- function(c0, k, s1, s2) {
  if (k == 0) return(c0 * (s2 - s1))
  c0 * exp(-k * s1) * -expm1(-k * (s2 - s1)) / k
}

twa_model <- function(est_g, occ, apps, a, b) {
  fail <- function(reason) list(twa = NA_real_, status = "not_estimable", reason = reason)
  if (a < apps[1L]) {
    return(fail(sprintf("the fitted curves start at the first application (time %s)",
                        format(apps[1L]))))
  }
  m <- length(apps)
  area <- 0
  for (p in seq_len(m)) {
    lo <- apps[p]
    hi <- if (p < m) apps[p + 1L] else Inf
    u <- max(a, lo)
    v <- min(b, hi)
    if (v <= u) next
    j <- which(est_g$period == p)
    if (length(j) == 0L) return(fail(sprintf("period %d was not fitted", p)))
    if (!is.finite(est_g$c0[j]) || !is.finite(est_g$k[j])) {
      return(fail(sprintf("the fit of period %d is not available (%s)", p, est_g$reason[j])))
    }
    area <- area + sfo_integral(est_g$c0[j], est_g$k[j], u - lo, v - lo)
  }
  measured <- occ$profile$time[!is.na(occ$profile$value) & occ$profile$period > 0L]
  if (length(measured) > 0L && b > max(measured)) {
    return(list(twa = area / (b - a), status = "extrapolated",
                reason = sprintf("the window ends after the last measured value (time %s)",
                                 format(max(measured)))))
  }
  list(twa = area / (b - a), status = "ok", reason = NA_character_)
}

#' Measured concentrations against the nominal ones
#'
#' `compare_nominal()` compares, for each application of each treatment, the
#' concentration measured right after it with the nominal concentration that
#' the application adds: a check of the dosing and of the build-up of the
#' substance with repeated applications.
#'
#' @param x A `meso_exposure` object.
#' @param by `"treatment"` (default; mean of the mesocosms of the treatment
#'   with a value) or `"mesocosm"`.
#' @param lag Longest time between an application and a sample for the sample
#'   to count as taken right after it (or right before it). The default, 0,
#'   takes only the samples taken at the application time; with times in
#'   days, `lag = 1/24` also takes those taken up to one hour later.
#'
#' @details
#' For application \eqn{p}, with nominal concentration \eqn{N_p}
#' (`x$applications$nominal`), mean concentration measured right after it
#' \eqn{A_p} and right before it \eqn{B_p}:
#' \deqn{pct\_nominal = 100 A_p / N_p, \qquad pct\_added = 100 (A_p - B_p) / N_p.}{
#'   pct_nominal = 100 Ap / Np,  pct_added = 100 (Ap - Bp) / Np.}
#' After the first application, `pct_nominal` measures the dosing. After later
#' ones, it includes what remained of the earlier applications, so it can
#' exceed 100 %; `pct_added` measures the dosing of that application alone,
#' and needs a sample right before it.
#'
#' Right after an application means the earliest sampling occasion within
#' `lag` after it, with at least one value; right before, the latest within
#' `lag` before it, among the samples of the previous period. Missing values
#' are left out of the means and never replaced.
#'
#' @return A data frame with one row per treatment (or mesocosm) and
#'   application: `treatment`, `concentration`, `mesocosm` (with
#'   `by = "mesocosm"`), `application` (1, 2, ...), `time` (of the
#'   application), `nominal`, `time_before`, `before`, `n_before` (values in
#'   the mean), `time_after`, `after`, `n_after`, `pct_nominal`, `pct_added`,
#'   `unit`, `status` and `reason`. `status` is `"ok"`, `"no_before"`
#'   (`pct_added` is `NA`: no value right before the application) or
#'   `"not_estimable"` (no value right after the application, or a nominal
#'   concentration that is unknown or zero).
#' @seealso [meso_exposure()], [calc_twa()] (`twa_pct_nominal`),
#'   [fit_dissipation()] (`c0_pct_nominal`).
#' @export
#' @examples
#' expo <- subset(azoxy_exposure, habitat == "eutrophic_lake")
#' e <- meso_exposure(expo, control = "0", applications = c(0, 14),
#'                    mesocosm_col = "mesocosm", time_col = "day",
#'                    treatment_col = "concentration",
#'                    concentration_col = "concentration",
#'                    value_col = "azoxystrobin", timing_col = "timing",
#'                    unit = "µg/L")
#' compare_nominal(e)[c("treatment", "application", "after", "pct_nominal",
#'                      "pct_added", "status")]
compare_nominal <- function(x, by = "treatment", lag = 0) {
  check_meso_exposure(x)
  check_choice(by, c("treatment", "mesocosm"))
  check_number(lag, lower = 0)
  s <- x$settings
  sch <- x$applications
  groups <- exposure_groups(x, by)
  gkey <- if (by == "treatment") "treatment" else "mesocosm"
  d <- x$data[x$data$treatment != s$control & !is.na(x$data$value), , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    dg <- d[d[[gkey]] == groups[[gkey]][i], , drop = FALSE]
    apps <- sch[sch$treatment == groups$treatment[i], , drop = FALSE]
    for (j in seq_len(nrow(apps))) {
      tp <- apps$time[j]
      eps <- 1e-8 * max(1, abs(tp), lag)
      aft <- dg[dg$period == j & dg$time <= tp + lag + eps, , drop = FALSE]
      bef <- dg[dg$period == j - 1L & dg$time >= tp - lag - eps & dg$time <= tp + eps, ,
                drop = FALSE]
      ta <- if (nrow(aft) > 0L) min(aft$time) else NA_real_
      tb <- if (nrow(bef) > 0L) max(bef$time) else NA_real_
      va <- aft$value[aft$time %in% ta]
      vb <- bef$value[bef$time %in% tb]
      after <- if (length(va) > 0L) mean(va) else NA_real_
      before <- if (length(vb) > 0L) mean(vb) else NA_real_
      nominal <- apps$nominal[j]
      status <- "ok"
      reason <- NA_character_
      if (is.na(nominal) || nominal == 0) {
        status <- "not_estimable"
        reason <- if (is.na(nominal)) {
          "the nominal concentration of the application is unknown"
        } else {
          "the nominal concentration of the application is zero"
        }
      } else if (is.na(after)) {
        status <- "not_estimable"
        reason <- "no value right after the application"
      } else if (is.na(before)) {
        status <- "no_before"
        reason <- "no value right before the application, so no pct_added"
      }
      ok <- status != "not_estimable"
      rows[[length(rows) + 1L]] <- data.frame(
        groups[i, , drop = FALSE], application = j, time = tp, nominal = nominal,
        time_before = tb, before = before, n_before = length(vb),
        time_after = ta, after = after, n_after = length(va),
        pct_nominal = if (ok) 100 * after / nominal else NA_real_,
        pct_added = if (ok) 100 * (after - before) / nominal else NA_real_,
        unit = if (is.null(s$unit)) NA_character_ else s$unit,
        status = status, reason = reason, stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot the measured exposure
#'
#' `plot_exposure()` draws the concentrations of a [meso_exposure()] object
#' over time, one panel per treatment, with the applications of the treatment
#' as dashed lines, its nominal concentration as a dotted line and,
#' optionally, the curves of [fit_dissipation()]. Needs the optional packages
#' 'ggplot2' and 'viridisLite'.
#'
#' @param x A `meso_exposure` object.
#' @param fit `NULL` (default) or the result of `fit_dissipation(x)`, whose
#'   curves are drawn. Without it, lines join the samples of each application
#'   period, the profile that `calc_twa()` integrates.
#' @param show `"mean_sd"` (default): mean and standard deviation of the
#'   mesocosms of each treatment; `"replicates"`: the value of each
#'   mesocosm; `"both"`.
#' @param scale `"linear"` (default; the y axis starts at zero) or `"log"`
#'   (logarithmic y axis; zeros are not drawn, and an error bar whose lower
#'   end would not be positive is drawn from the mean up).
#' @param show_control Draw the panel of the controls too (default `FALSE`);
#'   it is empty when the controls have no values.
#' @param show_nominal Draw the nominal concentration of each treatment
#'   (`concentration_col` of [meso_exposure()]) as a dotted line (default
#'   `TRUE`; nothing is drawn when it is unknown).
#' @param colours,labels `NULL` (defaults: the treatments on the viridis
#'   "mako" scale, from light to dark as the concentration increases, and the
#'   control in grey; labels with the concentrations) or one colour or label
#'   per treatment, control first, optionally named by treatment.
#' @param xlab Label of the time axis (default `"Time"`).
#'
#' @return A 'ggplot2' object.
#' @seealso [meso_exposure()], [fit_dissipation()], [calc_twa()].
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("viridisLite", quietly = TRUE)
#' expo <- subset(azoxy_exposure, habitat == "eutrophic_lake")
#' e <- meso_exposure(expo, control = "0", applications = c(0, 14),
#'                    mesocosm_col = "mesocosm", time_col = "day",
#'                    treatment_col = "concentration",
#'                    concentration_col = "concentration",
#'                    value_col = "azoxystrobin", timing_col = "timing",
#'                    unit = "µg/L")
#' plot_exposure(e, fit = fit_dissipation(e), xlab = "Day")
#' plot_exposure(e, show = "replicates", scale = "log", xlab = "Day")
plot_exposure <- function(x, fit = NULL, show = "mean_sd", scale = "linear",
                          show_control = FALSE, show_nominal = TRUE, colours = NULL,
                          labels = NULL, xlab = "Time") {
  check_meso_exposure(x)
  check_installed(c("ggplot2", "viridisLite"), "plot the exposure")
  check_choice(show, c("mean_sd", "replicates", "both"))
  check_choice(scale, c("linear", "log"))
  check_flag(show_control)
  check_flag(show_nominal)
  check_string(xlab)
  if (!is.null(fit)) {
    check_meso_dissipation(fit)
    check_same_exposure(fit, x)
  }
  s <- x$settings
  lev <- s$treatments
  conc <- x$design$concentration[match(lev, x$design$treatment)]
  col <- treatment_colours(lev, colours)
  lab <- exposure_labels(lev, conc, s$control, s$unit, labels)
  shown <- if (show_control) lev else lev[-1L]
  panel <- function(trt) factor(lab[trt], levels = unique(lab[shown]))

  d <- x$data[x$data$treatment %in% shown & !is.na(x$data$value), , drop = FALSE]
  if (scale == "log") d <- d[d$value > 0, , drop = FALSE]
  d$occasion <- paste(d$treatment, format_time_key(d$time), d$timing, sep = "\u001f")
  first <- !duplicated(d$occasion)
  means <- d[first, c("treatment", "time", "timing", "period"), drop = FALSE]
  means$mean <- as.numeric(tapply(d$value, d$occasion, mean)[d$occasion[first]])
  means$sd <- as.numeric(tapply(d$value, d$occasion, stats::sd)[d$occasion[first]])
  means$lower <- means$mean - means$sd
  if (scale == "log") means$lower <- ifelse(means$lower > 0, means$lower, means$mean)
  means$upper <- means$mean + means$sd
  means <- means[order(means$treatment, means$time,
                       match(means$timing, c("before", "after"), nomatch = 0L),
                       method = "radix"), , drop = FALSE]
  d <- d[order(d$mesocosm, d$time, match(d$timing, c("before", "after"), nomatch = 0L),
               method = "radix"), , drop = FALSE]
  means$panel <- panel(means$treatment)
  d$panel <- panel(d$treatment)
  means$treatment <- factor(means$treatment, levels = lev)
  d$treatment <- factor(d$treatment, levels = lev)
  means$line <- paste(means$treatment, means$period)
  d$line <- paste(d$mesocosm, d$period)

  apps <- x$applications[x$applications$treatment %in% shown, c("treatment", "time"),
                         drop = FALSE]
  apps$panel <- panel(apps$treatment)
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(data = apps, ggplot2::aes(xintercept = .data$time),
                        linetype = "dashed", colour = "#898781", linewidth = 0.4)
  nominal <- data.frame(treatment = shown, nominal = conc[match(shown, lev)],
                        stringsAsFactors = FALSE)
  nominal <- nominal[is.finite(nominal$nominal) & nominal$nominal > 0, , drop = FALSE]
  if (show_nominal && nrow(nominal) > 0L) {
    nominal$panel <- panel(nominal$treatment)
    nominal$x <- max(x$data$time)
    p <- p +
      ggplot2::geom_hline(data = nominal, ggplot2::aes(yintercept = .data$nominal),
                          linetype = "dotted", colour = "#898781", linewidth = 0.5) +
      ggplot2::geom_text(data = nominal, ggplot2::aes(x = .data$x, y = .data$nominal),
                         label = "nominal", hjust = 1, vjust = -0.5, size = 3,
                         colour = "#5f5e5a")
  }
  if (is.null(fit)) {
    if (show == "replicates") {
      p <- p + ggplot2::geom_line(data = d, ggplot2::aes(x = .data$time, y = .data$value,
                                                         colour = .data$treatment,
                                                         group = .data$line),
                                  linewidth = 0.4, alpha = 0.6)
    } else {
      p <- p + ggplot2::geom_line(data = means, ggplot2::aes(x = .data$time, y = .data$mean,
                                                             colour = .data$treatment,
                                                             group = .data$line),
                                  linewidth = 0.6)
    }
  } else {
    pred <- stats::predict(fit)
    pred <- pred[pred$treatment %in% shown, , drop = FALSE]
    pred$panel <- panel(pred$treatment)
    pred$treatment <- factor(pred$treatment, levels = lev)
    pred$line <- paste(if (is.null(pred$mesocosm)) pred$treatment else pred$mesocosm,
                       pred$period)
    p <- p + ggplot2::geom_line(data = pred, ggplot2::aes(x = .data$time, y = .data$fit,
                                                          colour = .data$treatment,
                                                          group = .data$line),
                                linewidth = 0.7)
  }
  if (show %in% c("replicates", "both")) {
    p <- p + ggplot2::geom_point(data = d, ggplot2::aes(x = .data$time, y = .data$value,
                                                        colour = .data$treatment),
                                 size = 1.4, alpha = if (show == "both") 0.45 else 0.8)
  }
  if (show %in% c("mean_sd", "both")) {
    bars <- means[!is.na(means$sd), , drop = FALSE]
    p <- p +
      ggplot2::geom_linerange(data = bars, ggplot2::aes(x = .data$time, ymin = .data$lower,
                                                        ymax = .data$upper,
                                                        colour = .data$treatment),
                              linewidth = 0.5) +
      ggplot2::geom_point(data = means, ggplot2::aes(x = .data$time, y = .data$mean,
                                                     colour = .data$treatment),
                          size = 2)
  }
  ylab <- if (is.null(s$unit)) "Concentration" else paste0("Concentration (", s$unit, ")")
  times <- sort(unique(x$data$time))
  p <- p +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_y", drop = FALSE) +
    ggplot2::scale_colour_manual(values = col, breaks = lev, guide = "none") +
    ggplot2::scale_x_continuous(breaks = if (length(times) <= 12L) times else
                                  ggplot2::waiver()) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(strip.background = ggplot2::element_blank())
  if (scale == "log") {
    p <- p + ggplot2::scale_y_log10()
  } else {
    p <- p + ggplot2::scale_y_continuous(limits = c(0, NA), oob = function(x, ...) x)
  }
  p
}
