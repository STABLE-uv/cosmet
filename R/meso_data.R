#' Mesocosm experiment data
#'
#' `meso_data()` validates a long table of a mesocosm (or microcosm) experiment,
#' with one row per mesocosm, sampling time and taxon (or any other measured
#' variable), and keeps it with the design and the decisions that the analyses
#' need. Every other function of the package takes its result.
#'
#' @param data A data frame with one row per mesocosm, sampling time and taxon.
#'   Tibbles are accepted; grouped data frames are not.
#' @param control The level of the treatment column that is the control
#'   (character; numbers are compared as text, so `"0"` for a column of
#'   concentrations).
#' @param zero_fill `TRUE` or `FALSE`, no default. With `TRUE`, a taxon without
#'   a row in an analysed sample was looked for and not found: it becomes an
#'   observed zero (`value_status = "zero_filled"`). With `FALSE`, it is
#'   unknown (`NA`, `value_status = "not_recorded"`). Use `TRUE` for community
#'   counts stored as presences only; use `FALSE` for variables other than
#'   abundances, or when absences are not known. With complete tables (every
#'   taxon in every sample) both give the same result.
#' @param mesocosm_col,time_col,treatment_col,taxon_col,value_col Names of the
#'   columns with the experimental unit, the sampling time (numeric, for
#'   example days relative to the first application), the treatment, the taxon
#'   or variable, and the measured value (numeric; `NA` for an unknown value).
#' @param concentration_col Name of the column with the nominal concentration
#'   assigned to each mesocosm (numeric, non-negative, zero in the control), or
#'   `NULL` (default) if there is none. Williams' test and NOECs need it.
#' @param samples `NULL` (default) or a data frame with the columns
#'   `mesocosm_col` and `time_col` listing every sample that was analysed,
#'   including those in which nothing was found. Without it, the analysed
#'   samples are those with at least one row in `data`.
#' @param unit Unit of `value` (for example `"individuals/L"`), or `NULL`.
#'   Used to label results and to convert abundance thresholds; see
#'   [s11_criteria()].
#' @param concentration_unit Unit of the concentrations (for example
#'   `"ug/L"`), or `NULL`. Used in labels.
#' @param applications Times of the applications of the chemical, on the scale
#'   of `time_col`, or `NULL`. Used in plots.
#'
#' @return An object of class `meso_data`: a list with
#'   * `data`: a data frame with one row per analysed sample and taxon, ordered
#'     by mesocosm, time and taxon (C-locale order for text), with columns
#'     `mesocosm` (character), `time` (numeric), `treatment` (character),
#'     `concentration` (numeric, `NA` without `concentration_col`), `taxon`
#'     (character), `value` (numeric, in `unit`), `response` (numeric, the
#'     scale of the analyses: `value` until [transform_response()] is applied)
#'     and `value_status` (`"recorded"`, `"missing"` for a recorded `NA`,
#'     `"zero_filled"` or `"not_recorded"`);
#'   * `design`: one row per mesocosm with its treatment and concentration;
#'   * `samples`: one row per analysed sample (mesocosm and time);
#'   * `settings`: the control, the order of the treatments, `zero_fill`,
#'     units, application times, the transformation and the original column
#'     names;
#'   * `diagnostics`: counts of zero-filled, missing and not recorded values.
#'
#' @section Missing values and zeros:
#' A recorded zero is an observed absence. A recorded `NA` stays unknown and is
#' never replaced by zero. A taxon without a row in an analysed sample is a
#' zero only with `zero_fill = TRUE`. A mesocosm that was not sampled at a
#' time is not filled at all: it is absent from `data` and listed by
#' [check_design()].
#'
#' @section Assumptions:
#' The mesocosm is the experimental unit: each one has a single treatment and
#' concentration for the whole experiment. Values of a taxon are comparable
#' between samples (same sampling effort and unit).
#'
#' @seealso [check_design()], [transform_response()], [test_control()].
#' @export
#' @examples
#' counts <- data.frame(
#'   mesocosm = rep(c("C1", "C2", "T1", "T2"), each = 2),
#'   time = 7,
#'   treatment = rep(c("control", "control", "high", "high"), each = 2),
#'   taxon = rep(c("Cloeon", "Chironomini"), 4),
#'   value = c(10, 4, 12, 6, 2, 5, 0, 7)
#' )
#' d <- meso_data(counts, control = "control", zero_fill = TRUE)
#' d
meso_data <- function(data, control, zero_fill, mesocosm_col = "mesocosm",
                      time_col = "time", treatment_col = "treatment",
                      taxon_col = "taxon", value_col = "value",
                      concentration_col = NULL, samples = NULL, unit = NULL,
                      concentration_unit = NULL, applications = NULL) {
  if (missing(control)) {
    abort_input(paste(
      "`control` is missing: give the level of the treatment column that is",
      "the control, for example control = \"0\" or control = \"C\"."
    ))
  }
  if (missing(zero_fill)) {
    abort_input(paste(
      "`zero_fill` is missing and has no default. Use zero_fill = TRUE if a",
      "taxon without a row in an analysed sample was looked for and not found",
      "(an observed zero), or zero_fill = FALSE if it is unknown."
    ))
  }
  check_data_frame(data, allow_empty = FALSE)
  check_flag(zero_fill)
  check_string(mesocosm_col)
  check_string(time_col)
  check_string(treatment_col)
  check_string(taxon_col)
  check_string(value_col)
  check_string(concentration_col, allow_null = TRUE)
  check_string(unit, allow_null = TRUE)
  check_string(concentration_unit, allow_null = TRUE)
  cols <- c(mesocosm_col = mesocosm_col, time_col = time_col,
            treatment_col = treatment_col, taxon_col = taxon_col,
            value_col = value_col, concentration_col = concentration_col)
  check_columns(data, cols)
  if (!is.null(applications)) {
    if (!is.numeric(applications) || length(applications) == 0L ||
        anyNA(applications) || any(!is.finite(applications))) {
      abort_input("`applications` must be NULL or a vector of finite numbers (times of application).")
    }
    applications <- sort(unique(as.numeric(applications)))
  }

  mesocosm <- key_text(data[[mesocosm_col]], mesocosm_col)
  treatment <- key_text(data[[treatment_col]], treatment_col)
  taxon <- key_text(data[[taxon_col]], taxon_col)
  time <- data[[time_col]]
  if (inherits(time, c("Date", "POSIXt"))) {
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
  check_numeric_column(data, value_col, allow_na = TRUE)
  concentration <- rep(NA_real_, nrow(data))
  if (!is.null(concentration_col)) {
    check_numeric_column(data, concentration_col, lower = 0)
    concentration <- as.numeric(data[[concentration_col]])
  }
  if (is.numeric(control) && length(control) == 1L && is.finite(control)) {
    control <- as.character(control)
  }
  check_string(control)

  obs <- data.frame(mesocosm = mesocosm, time = as.numeric(time),
                    treatment = treatment, concentration = concentration,
                    taxon = taxon, value = as.numeric(data[[value_col]]),
                    stringsAsFactors = FALSE)
  check_unique_key(obs, c("mesocosm", "time", "taxon"), arg = "data")

  design <- build_design(obs, control, has_concentration = !is.null(concentration_col),
                         cols = cols)
  treatments <- design_levels(design, control)

  # Analysed samples
  observed_samples <- unique(obs[c("mesocosm", "time")])
  if (is.null(samples)) {
    smp <- observed_samples
    samples_source <- "inferred from data"
  } else {
    check_data_frame(samples, allow_empty = FALSE)
    check_columns(samples, c(mesocosm_col = mesocosm_col, time_col = time_col),
                  arg = "samples")
    check_numeric_column(samples, time_col, arg = "samples")
    smp <- data.frame(mesocosm = key_text(samples[[mesocosm_col]], mesocosm_col),
                      time = as.numeric(samples[[time_col]]), stringsAsFactors = FALSE)
    smp <- unique(smp)
    unknown <- setdiff(smp$mesocosm, design$mesocosm)
    if (length(unknown) > 0L) {
      abort_input(sprintf(
        paste(
          "`samples` lists mesocosms without any row in `data`, so their",
          "treatment is unknown: %s. Add them to `data` or remove them."
        ),
        format_items(unknown)
      ))
    }
    k_obs <- paste(observed_samples$mesocosm, observed_samples$time, sep = "\u001f")
    k_smp <- paste(smp$mesocosm, smp$time, sep = "\u001f")
    not_listed <- which(!k_obs %in% k_smp)
    if (length(not_listed) > 0L) {
      abort_input(sprintf(
        "`data` has values for samples that are not in `samples`: %s.",
        format_items(paste0(observed_samples$mesocosm[not_listed], " at time ",
                            observed_samples$time[not_listed]))
      ))
    }
    samples_source <- "given"
  }
  smp <- smp[order(smp$mesocosm, smp$time, method = "radix"), , drop = FALSE]
  rownames(smp) <- NULL
  idx <- match(smp$mesocosm, design$mesocosm)
  smp$treatment <- design$treatment[idx]
  smp$concentration <- design$concentration[idx]

  # Complete grid: analysed samples x taxa
  taxa <- sort(unique(obs$taxon), method = "radix")
  grid <- data.frame(
    mesocosm = rep(smp$mesocosm, each = length(taxa)),
    time = rep(smp$time, each = length(taxa)),
    taxon = rep(taxa, times = nrow(smp)),
    stringsAsFactors = FALSE
  )
  key <- function(m, t, x) paste(m, format_time_key(t), x, sep = "\u001f")
  pos <- match(key(grid$mesocosm, grid$time, grid$taxon),
               key(obs$mesocosm, obs$time, obs$taxon))
  value <- obs$value[pos]
  status <- rep("recorded", nrow(grid))
  status[!is.na(pos) & is.na(value)] <- "missing"
  absent_row <- is.na(pos)
  if (zero_fill) {
    value[absent_row] <- 0
    status[absent_row] <- "zero_filled"
  } else {
    status[absent_row] <- "not_recorded"
  }
  idx <- match(grid$mesocosm, design$mesocosm)
  out <- data.frame(
    mesocosm = grid$mesocosm, time = grid$time,
    treatment = design$treatment[idx], concentration = design$concentration[idx],
    taxon = grid$taxon, value = value, response = value, value_status = status,
    stringsAsFactors = FALSE
  )

  settings <- list(
    control = control, treatments = treatments,
    has_concentration = !is.null(concentration_col), zero_fill = zero_fill,
    samples = samples_source, unit = unit, concentration_unit = concentration_unit,
    applications = applications,
    transform = list(method = "none", A = NA_real_),
    columns = as.list(cols)
  )
  diagnostics <- list(
    n_recorded = sum(status == "recorded"),
    n_missing = sum(status == "missing"),
    n_zero_filled = sum(status == "zero_filled"),
    n_not_recorded = sum(status == "not_recorded")
  )
  new_meso_data(out, design, smp, settings, diagnostics)
}

new_meso_data <- function(data, design, samples, settings, diagnostics) {
  structure(list(data = data, design = design, samples = samples,
                 settings = settings, diagnostics = diagnostics),
            class = "meso_data")
}

# Identifier columns as text without missing or empty values.
key_text <- function(x, column, call = sys.call(-1)) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x) && !is.numeric(x) && !is.logical(x)) {
    abort_input(sprintf("Column '%s' must hold text or numbers, not %s.",
                        column, class(x)[1L]), call = call)
  }
  x <- as.character(x)
  bad <- which(is.na(x) | !nzchar(trimws(x)))
  if (length(bad) > 0L) {
    abort_input(sprintf("Column '%s' has missing or empty values in %s.",
                        column, format_rows(bad)), call = call)
  }
  x
}

# Times as keys without floating-point formatting surprises.
format_time_key <- function(t) {
  formatC(t, digits = 15, format = "g")
}

# One row per mesocosm; checks that treatment and concentration are fixed.
build_design <- function(obs, control, has_concentration, cols, call = sys.call(-1)) {
  design <- unique(obs[c("mesocosm", "treatment", "concentration")])
  dup <- unique(design$mesocosm[duplicated(design$mesocosm)])
  if (length(dup) > 0L) {
    abort_input(sprintf(
      paste(
        "Each mesocosm must have a single treatment and concentration for the",
        "whole experiment; these have more than one: %s. Check columns '%s'%s."
      ),
      format_items(dup), cols[["treatment_col"]],
      if (has_concentration) sprintf(" and '%s'", cols[["concentration_col"]]) else ""
    ), call = call)
  }
  design <- design[order(design$mesocosm, method = "radix"), , drop = FALSE]
  rownames(design) <- NULL
  if (!control %in% design$treatment) {
    abort_input(sprintf(
      "`control` = %s is not a level of column '%s', whose levels are %s.",
      quote_items(control), cols[["treatment_col"]],
      format_items(quote_items(sort(unique(design$treatment), method = "radix")), max = 20L)
    ), call = call)
  }
  if (length(unique(design$treatment)) < 2L) {
    abort_input("The experiment needs at least one treatment besides the control.", call = call)
  }
  if (has_concentration) {
    per_trt <- tapply(design$concentration, design$treatment, function(x) length(unique(x)))
    bad <- names(per_trt)[per_trt > 1L]
    if (length(bad) > 0L) {
      abort_input(sprintf(
        "Treatments with more than one concentration: %s. Each treatment must have one nominal concentration.",
        format_items(bad)
      ), call = call)
    }
    c0 <- unique(design$concentration[design$treatment == control])
    if (!isTRUE(all(c0 == 0))) {
      abort_input(sprintf(
        "The control (%s) has concentration %s; the control concentration must be 0.",
        quote_items(control), format(c0)
      ), call = call)
    }
  }
  design
}

# Treatment levels: control first, then by concentration and name.
design_levels <- function(design, control) {
  trt <- unique(design[c("treatment", "concentration")])
  trt <- trt[trt$treatment != control, , drop = FALSE]
  trt <- trt[order(trt$concentration, trt$treatment, method = "radix"), , drop = FALSE]
  c(control, trt$treatment)
}

#' @export
print.meso_data <- function(x, ...) {
  s <- x$settings
  d <- x$data
  cat("<meso_data>\n")
  cat(sprintf("%d mesocosms, %d sampling time(s), %d taxa or variables, %d values\n",
              nrow(x$design), length(unique(d$time)), length(unique(d$taxon)), nrow(d)))
  trt <- s$treatments
  if (s$has_concentration) {
    conc <- format_conc(x$design$concentration[match(trt, x$design$treatment)])
    unit <- if (is.null(s$concentration_unit)) "" else paste0(" ", s$concentration_unit)
    trt <- ifelse(trt == conc, paste0(conc, unit), sprintf("%s (%s%s)", trt, conc, unit))
  }
  cat(sprintf("Treatments: %s; control: %s\n", paste(trt, collapse = ", "), s$control))
  cat(sprintf("Times: %s\n", paste(format(sort(unique(d$time)), trim = TRUE), collapse = ", ")))
  if (!is.null(s$applications)) {
    cat(sprintf("Applications at: %s\n", paste(format(s$applications, trim = TRUE), collapse = ", ")))
  }
  cat(sprintf("Unit: %s\n", if (is.null(s$unit)) "not given" else s$unit))
  cat(sprintf("Analysis scale: %s\n", describe_transform(s$transform)))
  dg <- x$diagnostics
  cat(sprintf("Values: %d recorded, %d zero-filled, %d missing, %d not recorded\n",
              dg$n_recorded, dg$n_zero_filled, dg$n_missing, dg$n_not_recorded))
  invisible(x)
}

#' @export
as.data.frame.meso_data <- function(x, ...) {
  x$data
}

describe_transform <- function(tr) {
  if (identical(tr$method, "log")) {
    sprintf("ln(%s x + 1)", format(tr$A, trim = TRUE, drop0trailing = TRUE))
  } else {
    "untransformed values"
  }
}

check_meso_data <- function(x, arg = arg_name(substitute(x)), call = sys.call(-1)) {
  force(arg)
  if (!inherits(x, "meso_data")) {
    abort_input(sprintf("`%s` must be a meso_data object created with meso_data(), not %s.",
                        arg, describe_value(x)), call = call)
  }
  invisible(x)
}

#' Check the design of a mesocosm experiment
#'
#' `check_design()` reports the replication of every treatment at every
#' sampling time, the samples that are missing (a mesocosm not sampled at a
#' time when others were), and the values that are missing or were filled with
#' zeros. It changes nothing.
#'
#' @param x A [meso_data()] object.
#'
#' @return An object of class `meso_design_check`, a list with `replicates`
#'   (a data frame with one row per treatment and columns `treatment`,
#'   `concentration` and one column per time with the number of sampled
#'   mesocosms), `missing_samples` (mesocosm, treatment and time of every
#'   sample not taken), `value_status` (counts per status), `balanced`
#'   (`TRUE` if every mesocosm was sampled at every time) and `issues`
#'   (character vector of the problems found, empty if none). It prints as a
#'   report.
#' @export
#' @examples
#' counts <- data.frame(
#'   mesocosm = c("C1", "C2", "T1", "T2", "C1", "T1"),
#'   time = c(7, 7, 7, 7, 14, 14),
#'   treatment = c("control", "control", "high", "high", "control", "high"),
#'   taxon = "Cloeon",
#'   value = c(10, 12, 2, 0, 9, 3)
#' )
#' check_design(meso_data(counts, control = "control", zero_fill = TRUE))
check_design <- function(x) {
  check_meso_data(x)
  smp <- x$samples
  trt <- x$settings$treatments
  times <- sort(unique(smp$time))
  design <- x$design
  counts <- vapply(times, function(t) {
    m <- smp$mesocosm[smp$time == t]
    vapply(trt, function(tr) sum(design$treatment[match(m, design$mesocosm)] == tr), integer(1))
  }, integer(length(trt)))
  counts <- matrix(counts, nrow = length(trt), dimnames = list(trt, format(times, trim = TRUE)))
  replicates <- data.frame(treatment = trt,
                           concentration = design$concentration[match(trt, design$treatment)],
                           counts, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(replicates) <- NULL
  full <- expand.grid(time = times, mesocosm = design$mesocosm, stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)
  k_full <- paste(full$mesocosm, format_time_key(full$time), sep = "\u001f")
  k_smp <- paste(smp$mesocosm, format_time_key(smp$time), sep = "\u001f")
  miss <- full[!k_full %in% k_smp, c("mesocosm", "time"), drop = FALSE]
  miss$treatment <- design$treatment[match(miss$mesocosm, design$mesocosm)]
  miss <- miss[order(miss$time, miss$mesocosm, method = "radix"), c("mesocosm", "treatment", "time"),
               drop = FALSE]
  rownames(miss) <- NULL
  st <- x$data$value_status
  status <- table(factor(st, levels = c("recorded", "zero_filled", "missing", "not_recorded")))
  issues <- character()
  if (nrow(miss) > 0L) {
    issues <- c(issues, sprintf("%d sample(s) not taken (a mesocosm not sampled at a time).", nrow(miss)))
  }
  if (any(counts < 2L)) {
    issues <- c(issues, "Some treatment has fewer than 2 replicates at some time: no test is possible there.")
  }
  if (status[["missing"]] > 0L) {
    issues <- c(issues, sprintf("%d recorded value(s) are NA.", status[["missing"]]))
  }
  if (status[["not_recorded"]] > 0L) {
    issues <- c(issues, sprintf(
      "%d taxon x sample combination(s) have no row and are unknown (zero_fill = FALSE).",
      status[["not_recorded"]]))
  }
  structure(list(replicates = replicates, missing_samples = miss,
                 value_status = as.list(status), balanced = nrow(miss) == 0L,
                 issues = issues),
            class = "meso_design_check")
}

#' @export
print.meso_design_check <- function(x, ...) {
  cat("<meso_design_check>\n")
  cat("Replicates (mesocosms sampled) per treatment and time:\n")
  print(x$replicates, row.names = FALSE)
  vs <- x$value_status
  cat(sprintf("\nValues: %d recorded, %d zero-filled, %d missing, %d not recorded\n",
              vs$recorded, vs$zero_filled, vs$missing, vs$not_recorded))
  if (length(x$issues) == 0L) {
    cat("No problems found.\n")
  } else {
    cat("Problems:\n")
    cat(paste0("* ", x$issues, "\n"), sep = "")
  }
  if (nrow(x$missing_samples) > 0L) {
    cat("Samples not taken:\n")
    print(utils::head(x$missing_samples, 10L), row.names = FALSE)
    if (nrow(x$missing_samples) > 10L) {
      cat(sprintf("... and %d more\n", nrow(x$missing_samples) - 10L))
    }
  }
  invisible(x)
}

#' Transform the values to the scale of the analyses
#'
#' `transform_response()` sets the scale on which the tests, the MDD and the
#' principal response curves are computed: \eqn{y = \ln(A x + 1)}. The
#' original values are kept in `value`; the transformed ones go to
#' `response`, and the factor `A` is recorded, because the MDD needs it to
#' return to the original scale.
#'
#' @param x A [meso_data()] object.
#' @param A Positive factor of the transformation, in the inverse of the unit
#'   of the values. No default: van den Brink et al. (2000) choose it from the
#'   lowest non-zero abundance of the data set (see Brock et al. 2015, p.
#'   1162); Amador et al. (2024) used 33.3 for phytoplankton, 10 for
#'   zooplankton and 2 for macroinvertebrates (also in [s11_thresholds]).
#'
#' @return `x` with `response = log(A * value + 1)` and
#'   `settings$transform = list(method = "log", A = A)`. Applying it again
#'   replaces the previous transformation (it is always computed from
#'   `value`).
#'
#' @section Missing values and zeros:
#' Zeros become zero; `NA` stays `NA`. Negative values are an error.
#'
#' @references
#' van den Brink PJ, Hattink J, Bransen F, Van Donk E, Brock TCM (2000).
#' Impact of the fungicide carbendazim in freshwater microcosms. II.
#' Zooplankton, primary producers and final conclusions. *Aquatic Toxicology*
#' 48, 251–264. \doi{10.1016/S0166-445X(99)00037-5}
#'
#' Brock TCM, Hammers-Wirtz M, Hommen U, Preuss TG, Ratte H-T, Roessink I,
#' Strauss T, Van den Brink PJ (2015). The minimum detectable difference (MDD)
#' and the interpretation of treatment-related effects of pesticides in
#' experimental ecosystems. *Environmental Science and Pollution Research* 22,
#' 1160–1174. \doi{10.1007/s11356-014-3398-2}
#'
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @export
#' @examples
#' counts <- data.frame(
#'   mesocosm = c("C1", "C2", "T1", "T2"), time = 7,
#'   treatment = c("control", "control", "high", "high"),
#'   taxon = "Cloeon", value = c(10, 12, 2, 0)
#' )
#' d <- transform_response(meso_data(counts, control = "control", zero_fill = TRUE), A = 2)
#' d$data[c("value", "response")]
transform_response <- function(x, A) {
  check_meso_data(x)
  if (missing(A)) {
    abort_input(paste(
      "`A` is missing and has no default: choose the factor of ln(A x + 1)",
      "for your data (for example 2 for counts of macroinvertebrates)."
    ))
  }
  check_number(A, lower = 0, include_lower = FALSE)
  v <- x$data$value
  neg <- which(!is.na(v) & v < 0)
  if (length(neg) > 0L) {
    abort_input(sprintf(
      "ln(A x + 1) needs non-negative values; %d value(s) are negative (taxa %s).",
      length(neg), format_items(unique(x$data$taxon[neg]))
    ))
  }
  x$data$response <- log1p(A * v)
  x$settings$transform <- list(method = "log", A = A)
  x
}

#' Taxa that can be evaluated
#'
#' `taxon_eligibility()` summarises how abundant and frequent each taxon is
#' in the controls, and marks as eligible those that meet the criteria given.
#' The criteria use the controls only, so a taxon is not excluded because it
#' declined in the treatments.
#'
#' @param x A [meso_data()] object.
#' @param min_control_mean `NULL` (default, not used) or the minimum mean value
#'   in the control samples over all sampling times, in the unit of the data.
#' @param min_control_frequency `NULL` (default, not used) or the minimum
#'   proportion (0 to 1) of control samples in which the taxon was found.
#'
#' @return A data frame with one row per taxon (C-locale order) and columns
#'   `taxon`, `n_samples` (samples with a known value), `frequency`
#'   (proportion of those with a value above zero, all treatments),
#'   `control_mean` (mean of the known values in the controls, all times),
#'   `control_max_time_mean` (largest of the control means of each time),
#'   `control_frequency`, `n_missing` (unknown values) and `eligible`
#'   (logical: all the criteria given are met; `TRUE` for all taxa when no
#'   criterion is given; `NA` if a criterion cannot be evaluated because the
#'   taxon has no known value in the controls).
#'
#' @section Missing values and zeros:
#' Unknown values (`NA`) are left out of the means and frequencies and counted
#' in `n_missing`. Zeros count as observed absences.
#' @export
#' @examples
#' counts <- data.frame(
#'   mesocosm = rep(c("C1", "C2", "T1", "T2"), each = 2), time = 7,
#'   treatment = rep(c("control", "control", "high", "high"), each = 2),
#'   taxon = rep(c("Cloeon", "Chironomini"), 4),
#'   value = c(10, 0, 12, 1, 2, 0, 0, 0)
#' )
#' d <- meso_data(counts, control = "control", zero_fill = TRUE)
#' taxon_eligibility(d, min_control_mean = 3)
taxon_eligibility <- function(x, min_control_mean = NULL, min_control_frequency = NULL) {
  check_meso_data(x)
  check_number(min_control_mean, lower = 0, allow_null = TRUE)
  check_number(min_control_frequency, lower = 0, upper = 1, allow_null = TRUE)
  d <- x$data
  is_control <- d$treatment == x$settings$control
  taxa <- sort(unique(d$taxon), method = "radix")
  rows <- split(seq_len(nrow(d)), factor(d$taxon, levels = taxa))
  stat <- lapply(rows, function(i) {
    v <- d$value[i]
    ctl <- is_control[i]
    known <- !is.na(v)
    vc <- v[ctl & known]
    tc <- d$time[i][ctl & known]
    time_means <- if (length(vc) > 0L) tapply(vc, tc, mean) else numeric()
    c(n_samples = sum(known),
      frequency = if (any(known)) mean(v[known] > 0) else NA_real_,
      control_mean = if (length(vc) > 0L) mean(vc) else NA_real_,
      control_max_time_mean = if (length(time_means) > 0L) max(time_means) else NA_real_,
      control_frequency = if (length(vc) > 0L) mean(vc > 0) else NA_real_,
      n_missing = sum(!known))
  })
  m <- do.call(rbind, stat)
  out <- data.frame(taxon = taxa, n_samples = as.integer(m[, "n_samples"]),
                    frequency = m[, "frequency"], control_mean = m[, "control_mean"],
                    control_max_time_mean = m[, "control_max_time_mean"],
                    control_frequency = m[, "control_frequency"],
                    n_missing = as.integer(m[, "n_missing"]), stringsAsFactors = FALSE)
  eligible <- rep(TRUE, nrow(out))
  if (!is.null(min_control_mean)) {
    eligible <- eligible & out$control_mean >= min_control_mean
  }
  if (!is.null(min_control_frequency)) {
    eligible <- eligible & out$control_frequency >= min_control_frequency
  }
  out$eligible <- eligible
  rownames(out) <- NULL
  out
}
