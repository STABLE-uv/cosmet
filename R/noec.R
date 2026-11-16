#' No observed effect concentrations
#'
#' `derive_noec()` turns the tests of [test_control()] into a lowest observed
#' effect concentration (LOEC) and a no observed effect concentration (NOEC)
#' for every taxon and sampling time, and, optionally, marks the effects that
#' meet criteria for a consistent, treatment-related response.
#'
#' @param x A `meso_test` object from [test_control()].
#' @param pre_times `NULL` (default) or the sampling times before the first
#'   application. Their results are reported but never selected, and the
#'   criterion of consistency with the pre-treatment period uses them.
#' @param criteria `NULL` (default: an effect is selected when it is
#'   significant, in a single direction and after the application) or a set
#'   of criteria from [noec_criteria()] or [s11_criteria()].
#'
#' @return An object of class `meso_noec`, a list with `noec` (a data frame
#'   with one row per taxon and time, in the order of `x$cells`) and
#'   `settings`. The columns of `noec` are `taxon`, `time`, `status` (that of
#'   [test_control()], or `"ambiguous_direction"` when effects in both
#'   directions are significant, or `"undefined_step"` when a step of
#'   Williams' procedure could not be computed), `reason`, `direction` (of
#'   the effect detected), `p_value` (Williams: adjusted p-value of the test
#'   of the highest concentration, the smallest of the two directions;
#'   Dunnett: smallest adjusted p-value), `p_loec`, `loec`, `noec` (numeric,
#'   concentration unit; `NA` when there is an effect at the lowest
#'   concentration), `noec_label` (`">200"` when there is no effect at the
#'   highest concentration, `"<2"` when there is one at the lowest, otherwise
#'   the NOEC), `consistent` (every concentration above the LOEC also shows
#'   the effect; always `TRUE` with Williams), `mean_control_value`,
#'   `mean_for_threshold`, `concentration_for_threshold`, `threshold` and
#'   `threshold_unit` (criterion b, in the unit of the data), `crit_a` to
#'   `crit_d` (logical; `NA` when not evaluated), `is_post` and `selected`.
#'
#' @details
#' The LOEC is the lowest concentration at which the effect is significant
#' in the direction detected, and the NOEC the tested concentration just
#' below it. With Williams' test the direction is that of the procedure that
#' rejects at the highest concentration. With Dunnett's test, a significant
#' concentration with a non-significant one above it gives
#' `consistent = FALSE`.
#'
#' The criteria (see [noec_criteria()]) are:
#' a. a significant effect in a single direction;
#' b. enough abundance: for decreases, the mean of the control; for
#'    increases, the largest mean of the concentrations with an effect that
#'    exceed the control, at least `min_abundance`;
#' c. means on the original scale that follow the direction of the effect
#'    from the control to the highest concentration (adjacent reversals
#'    smaller than `monotonic_tol` allowed) and differ at the extremes;
#' d. no effect in the same direction at the pre-treatment samplings
#'    (`NA`, so not selected, if one of them could not be tested).
#'
#' @section Missing values and zeros:
#' Taxa and times that could not be tested keep their `status` and have `NA`
#' NOEC and criteria.
#'
#' @references
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828} (criteria a to d: Supplementary
#' Text S11).
#'
#' @seealso [noec_table()], [plot_noec()].
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" &
#'                   habitat == "eutrophic_lake" & taxon %in% c("Chironomini", "Cloeon"))
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance",
#'                unit = "individuals/sample")
#' d <- transform_response(d, A = 2)
#' tw <- test_control(d, method = "williams", n_sim = 9999)
#' noec <- derive_noec(tw, pre_times = -7, criteria = s11_criteria("macroinvertebrates"))
#' noec
#' noec_table(noec)
derive_noec <- function(x, pre_times = NULL, criteria = NULL) {
  check_meso_test(x)
  s <- x$settings
  if (all(is.na(s$concentrations))) {
    abort_input(paste(
      "NOECs need concentrations: give `concentration_col` in meso_data().",
      "Without them, read the comparisons with the control in the tests of",
      "test_control()."
    ))
  }
  if (anyDuplicated(s$concentrations)) {
    abort_input("NOECs need one treatment per concentration; some treatments share a concentration.")
  }
  times <- sort(unique(x$cells$time))
  if (!is.null(pre_times)) {
    if (!is.numeric(pre_times) || length(pre_times) == 0L || anyNA(pre_times)) {
      abort_input("`pre_times` must be NULL or a vector of sampling times.")
    }
    pre_times <- sort(unique(as.numeric(pre_times)))
    unknown <- pre_times[!pre_times %in% times]
    if (length(unknown) > 0L) {
      abort_input(sprintf("`pre_times` has times without data: %s. Sampling times are %s.",
                          format_items(unknown), format_items(times, max = 20L)))
    }
  }
  if (!is.null(criteria) && !inherits(criteria, "cosmet_criteria")) {
    abort_input("`criteria` must be NULL or the result of noec_criteria() or s11_criteria().")
  }
  if (!is.null(criteria) && isTRUE(criteria$pretreatment) && is.null(pre_times)) {
    abort_input(paste(
      "The criterion of consistency with the pre-treatment period needs",
      "`pre_times`, the sampling times before the first application."
    ))
  }
  threshold <- criteria_threshold(criteria, s$unit)

  cells <- x$cells
  tests <- x$tests
  key_t <- paste(tests$taxon, format_time_key(tests$time), sep = "\u001f")
  key_c <- paste(cells$taxon, format_time_key(cells$time), sep = "\u001f")
  rows_of <- split(seq_len(nrow(tests)), factor(key_t, levels = key_c))
  out <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    out[[i]] <- noec_cell(tests[rows_of[[i]], , drop = FALSE], cells$status[i],
                          cells$reason[i], s, criteria, threshold)
  }
  noec <- cbind(cells[c("taxon", "time")], do.call(rbind, out))
  rownames(noec) <- NULL
  pre_end <- if (is.null(pre_times)) -Inf else max(pre_times)
  noec$is_post <- noec$time > pre_end
  if (!is.null(criteria) && isTRUE(criteria$pretreatment)) {
    noec$crit_d <- pretreatment_criterion(noec, tests, pre_times, s$method)
  }
  use <- noec$is_post & noec$status == "ok" & noec$crit_a %in% TRUE
  if (!is.null(criteria)) {
    if (!is.null(criteria$min_abundance)) use <- use & noec$crit_b %in% TRUE
    if (isTRUE(criteria$monotonic)) use <- use & noec$crit_c %in% TRUE
    if (isTRUE(criteria$pretreatment)) use <- use & noec$crit_d %in% TRUE
  }
  noec$selected <- use
  settings <- c(s[c("method", "direction", "alpha", "transform", "unit", "concentration_unit",
                    "control", "treatments", "concentrations", "applications")],
                list(pre_times = pre_times, criteria = criteria))
  structure(list(noec = noec, settings = settings), class = "meso_noec")
}

check_meso_test <- function(x, arg = arg_name(substitute(x)), call = sys.call(-1)) {
  force(arg)
  if (!inherits(x, "meso_test")) {
    abort_input(sprintf("`%s` must be the result of test_control(), not %s.",
                        arg, describe_value(x)), call = call)
  }
  invisible(x)
}

# One taxon and time. `tab`: its rows of the test table.
noec_cell <- function(tab, status, reason, s, criteria, threshold) {
  conc <- sort(unique(tab$concentration))
  out <- list(status = status, reason = reason, direction = NA_character_,
              p_value = NA_real_, p_loec = NA_real_, loec = NA_real_, noec = NA_real_,
              noec_label = NA_character_, consistent = NA,
              mean_control_value = tab$mean_control_value[1L],
              mean_for_threshold = NA_real_, concentration_for_threshold = NA_real_,
              threshold = if (is.null(threshold)) NA_real_ else threshold$value,
              threshold_unit = if (is.null(threshold)) NA_character_ else threshold$unit,
              crit_a = NA, crit_b = NA, crit_c = NA, crit_d = NA)
  if (identical(status, "ok")) {
    eff <- detected_effects(tab, s$method, conc)
    out$p_value <- eff$p_value
    out$crit_a <- eff$crit_a
    if (identical(eff$crit_a, FALSE)) {
      out$noec <- max(conc)
      out$noec_label <- paste0(">", format_conc(max(conc)))
    } else if (length(eff$directions) > 1L) {
      out$status <- "ambiguous_direction"
      out$reason <- "significant effects in both directions"
    } else if (length(eff$rejected) > 0L) {
      dir <- eff$directions
      out$direction <- dir
      path <- eff$path
      loec <- min(eff$rejected)
      out$p_loec <- path$p_adjusted[match(loec, path$concentration)]
      if (any(path$tested & is.na(path$significant))) {
        out$status <- "undefined_step"
        out$reason <- "a step of the procedure has no residual variance"
        out$p_loec <- NA_real_
      } else {
        out$loec <- loec
        out$consistent <- all(conc[conc >= loec] %in% eff$rejected)
        if (loec == min(conc)) {
          out$noec_label <- paste0("<", format_conc(min(conc)))
        } else {
          out$noec <- conc[match(loec, conc) - 1L]
          out$noec_label <- format_conc(out$noec)
        }
        if (!is.null(criteria)) {
          means <- c(path$mean_control_value[1L],
                     path$mean_treatment_value[match(conc, path$concentration)])
          all_conc <- c(0, conc)
          if (!is.null(threshold)) {
            if (dir == "decrease") {
              out$mean_for_threshold <- means[1L]
              out$concentration_for_threshold <- 0
            } else {
              ix <- which(all_conc %in% eff$rejected & means > means[1L])
              if (length(ix) > 0L) {
                j <- ix[which.max(means[ix])]
                out$mean_for_threshold <- means[j]
                out$concentration_for_threshold <- all_conc[j]
              }
            }
            out$crit_b <- is.finite(out$mean_for_threshold) &&
              out$mean_for_threshold >= threshold$value
          }
          if (isTRUE(criteria$monotonic)) {
            sgn <- if (dir == "increase") 1 else -1
            out$crit_c <- all(sgn * diff(means) >= -criteria$monotonic_tol) &&
              sgn * (means[length(means)] - means[1L]) > 0
          }
        }
      }
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Effects detected in one taxon and time.
detected_effects <- function(tab, method, conc) {
  if (method == "williams") {
    top <- tab[tab$concentration == max(conc), , drop = FALSE]
    valid <- all(!is.na(top$p_adjusted))
    dirs <- top$direction[top$significant %in% TRUE]
    p_value <- if (valid) min(top$p_adjusted) else NA_real_
    rejected <- numeric()
    path <- NULL
    if (length(dirs) == 1L) {
      path <- tab[tab$direction == dirs, , drop = FALSE]
      rejected <- path$concentration[path$significant %in% TRUE]
    }
  } else {
    valid <- all(!is.na(tab$p_adjusted))
    effect_dir <- ifelse(tab$direction == "both",
                         ifelse(tab$difference > 0, "increase", "decrease"),
                         tab$direction)
    sig <- tab$significant %in% TRUE
    dirs <- unique(effect_dir[sig])
    p_value <- if (valid) min(tab$p_adjusted) else NA_real_
    rejected <- numeric()
    path <- tab
    if (length(dirs) == 1L) {
      rejected <- tab$concentration[sig & effect_dir == dirs]
    }
  }
  list(crit_a = if (valid) length(dirs) > 0L else NA, directions = dirs,
       rejected = rejected, path = path, p_value = p_value)
}

# Criterion d: no effect in the same direction before the application.
pretreatment_criterion <- function(noec, tests, pre_times, method) {
  crit_d <- rep(NA, nrow(noec))
  eligible <- which(noec$is_post & noec$crit_a %in% TRUE & !is.na(noec$direction))
  for (i in eligible) {
    same <- vapply(pre_times, function(t) {
      tab <- tests[tests$taxon == noec$taxon[i] & tests$time == t, , drop = FALSE]
      if (nrow(tab) == 0L) return(NA)
      dir <- noec$direction[i]
      if (method == "williams") {
        top <- tab[tab$concentration == max(tab$concentration) & tab$direction == dir, ,
                   drop = FALSE]
        if (nrow(top) != 1L) return(NA)
        top$significant
      } else {
        effect_dir <- ifelse(tab$direction == "both",
                             ifelse(tab$difference > 0, "increase", "decrease"),
                             tab$direction)
        sig <- tab$significant
        if (anyNA(sig)) return(NA)
        any(sig & effect_dir == dir)
      }
    }, logical(1))
    crit_d[i] <- if (any(same, na.rm = TRUE)) {
      FALSE
    } else if (anyNA(same)) {
      NA
    } else {
      TRUE
    }
  }
  crit_d
}

format_conc <- function(x) {
  format(x, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)
}

#' Criteria for consistent, treatment-related effects
#'
#' `noec_criteria()` collects the criteria that [derive_noec()] applies to
#' select the effects to report. `s11_criteria()` gives the four criteria of
#' Supplementary Text S11 of Amador et al. (2024), with the abundance
#' threshold of a community from [s11_thresholds].
#'
#' @param min_abundance `NULL` (default, criterion b not applied) or the
#'   minimum mean abundance (control mean for decreases, mean of the affected
#'   concentration for increases).
#' @param abundance_unit Unit of `min_abundance`, or `NULL` for the unit of
#'   the data. When both are known and differ, `min_abundance` is converted
#'   (`"cells/mL"` and `"cells/L"`; `"individuals/mL"`, `"individuals/L"` and
#'   `"individuals/m3"`); other differences are an error.
#' @param monotonic If `TRUE`, criterion c: means on the original scale must
#'   follow the direction of the effect. Default `FALSE`.
#' @param monotonic_tol Largest reversal allowed between adjacent means in
#'   criterion c, in the unit of the data (default 0, none).
#' @param pretreatment If `TRUE`, criterion d: no effect in the same direction
#'   at the pre-treatment samplings (`pre_times` of [derive_noec()]). Default
#'   `FALSE`.
#' @param community One of `"phytoplankton_large"`, `"phytoplankton_small"`,
#'   `"microzooplankton"`, `"macrozooplankton"` and `"macroinvertebrates"`
#'   (see [s11_thresholds]).
#'
#' @return An object of class `cosmet_criteria`: a list with the arguments
#'   and `source`.
#' @references
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @export
#' @examples
#' noec_criteria(min_abundance = 5, monotonic = TRUE)
#' s11_criteria("macroinvertebrates")
noec_criteria <- function(min_abundance = NULL, abundance_unit = NULL, monotonic = FALSE,
                          monotonic_tol = 0, pretreatment = FALSE) {
  check_number(min_abundance, lower = 0, allow_null = TRUE)
  check_string(abundance_unit, allow_null = TRUE)
  if (!is.null(abundance_unit) && is.null(min_abundance)) {
    abort_input("`abundance_unit` needs `min_abundance`.")
  }
  check_flag(monotonic)
  check_number(monotonic_tol, lower = 0)
  check_flag(pretreatment)
  structure(list(min_abundance = min_abundance, abundance_unit = abundance_unit,
                 monotonic = monotonic, monotonic_tol = monotonic_tol,
                 pretreatment = pretreatment, source = "user"),
            class = "cosmet_criteria")
}

#' @rdname noec_criteria
#' @export
s11_criteria <- function(community) {
  if (missing(community)) {
    abort_input(sprintf("`community` is missing: use one of %s.",
                        format_items(quote_items(cosmet::s11_thresholds$community), max = Inf)))
  }
  check_choice(community, cosmet::s11_thresholds$community)
  row <- cosmet::s11_thresholds[cosmet::s11_thresholds$community == community, , drop = FALSE]
  out <- noec_criteria(min_abundance = row$min_abundance, abundance_unit = row$unit,
                       monotonic = TRUE, monotonic_tol = 0, pretreatment = TRUE)
  out$source <- sprintf("Amador et al. (2024), Text S11: %s", community)
  out
}

#' @export
print.cosmet_criteria <- function(x, ...) {
  cat("<cosmet_criteria>", if (identical(x$source, "user")) "given by the user" else x$source)
  cat("\n")
  cat("a. significant effect in a single direction\n")
  cat(sprintf("b. minimum abundance: %s\n", if (is.null(x$min_abundance)) "not applied" else
    paste(format(x$min_abundance), if (is.null(x$abundance_unit)) "(unit of the data)" else x$abundance_unit)))
  cat(sprintf("c. monotonic means: %s\n", if (x$monotonic) sprintf("yes (tolerance %s)", format(x$monotonic_tol)) else "not applied"))
  cat(sprintf("d. no effect before the application: %s\n", if (x$pretreatment) "yes" else "not applied"))
  invisible(x)
}

# Threshold of criterion b in the unit of the data.
criteria_threshold <- function(criteria, data_unit, call = sys.call(-1)) {
  if (is.null(criteria) || is.null(criteria$min_abundance)) return(NULL)
  from <- criteria$abundance_unit
  if (is.null(from)) {
    return(list(value = criteria$min_abundance,
                unit = if (is.null(data_unit)) NA_character_ else data_unit))
  }
  if (is.null(data_unit)) {
    abort_input(sprintf(
      paste(
        "The abundance threshold is in %s but the unit of the data is unknown.",
        "Give `unit` in meso_data() so that it can be checked."
      ),
      from
    ), call = call)
  }
  list(value = convert_abundance(criteria$min_abundance, from, data_unit, call = call),
       unit = data_unit)
}

# Density units that can be converted into each other: counted thing and
# volume (litres).
abundance_units <- data.frame(
  unit = c("cells/mL", "cells/L", "individuals/mL", "individuals/L", "individuals/m3",
           "individuals/sample"),
  what = c("cells", "cells", "individuals", "individuals", "individuals", "individuals"),
  per = c("volume", "volume", "volume", "volume", "volume", "sample"),
  litres = c(1e-3, 1, 1e-3, 1, 1e3, NA),
  stringsAsFactors = FALSE
)

convert_abundance <- function(x, from, to, call = sys.call(-1)) {
  if (identical(from, to)) return(x)
  i <- match(from, abundance_units$unit)
  j <- match(to, abundance_units$unit)
  ok <- !is.na(i) && !is.na(j) && abundance_units$what[i] == abundance_units$what[j] &&
    abundance_units$per[i] == "volume" && abundance_units$per[j] == "volume"
  if (!ok) {
    abort_input(sprintf(
      paste(
        "Cannot convert an abundance threshold in %s to %s, the unit of the",
        "data. Give `min_abundance` in %s with noec_criteria()."
      ),
      from, to, to
    ), call = call)
  }
  x * abundance_units$litres[j] / abundance_units$litres[i]
}

#' Tables of NOECs for reports
#'
#' `noec_table()` lays out the NOECs of [derive_noec()] or [community_noec()]
#' with one row per taxon and one column per sampling time.
#'
#' @param x A `meso_noec` object.
#' @param selected If `FALSE` (default), every NOEC; if `TRUE`, only the
#'   effects selected (the other cells are `NA`) and only the taxa with at
#'   least one.
#' @param arrows If `TRUE` (default), an arrow after the NOEC of the effects
#'   with a direction: up for increases, down for decreases.
#'
#' @return A data frame with a column `taxon` and one character column per
#'   sampling time (named with the time). Cells are `NA` where no NOEC could
#'   be derived.
#' @export
#' @examples
#' macro <- subset(azoxy, community == "macroinvertebrates" &
#'                   habitat == "eutrophic_lake" & taxon %in% c("Chironomini", "Cloeon"))
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' d <- transform_response(d, A = 2)
#' noec <- derive_noec(test_control(d, method = "dunnett"), pre_times = -7)
#' noec_table(noec)
noec_table <- function(x, selected = FALSE, arrows = TRUE) {
  if (!inherits(x, "meso_noec")) {
    abort_input("`x` must be the result of derive_noec() or community_noec().")
  }
  check_flag(selected)
  check_flag(arrows)
  n <- x$noec
  label <- n$noec_label
  if (arrows) {
    arrow <- rep("", nrow(n))
    arrow[n$direction %in% "increase"] <- "\u2191"
    arrow[n$direction %in% "decrease"] <- "\u2193"
    has <- !is.na(label)
    label[has] <- paste0(label[has], arrow[has])
  }
  if (selected) {
    label[!n$selected] <- NA_character_
  }
  taxa <- unique(n$taxon)
  times <- sort(unique(n$time))
  out <- data.frame(taxon = taxa, stringsAsFactors = FALSE)
  for (t in times) {
    rows <- n$time == t
    out[[format(t, trim = TRUE)]] <- label[rows][match(taxa, n$taxon[rows])]
  }
  if (selected) {
    keep <- rowSums(!is.na(out[-1L])) > 0L
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' @export
print.meso_noec <- function(x, ...) {
  s <- x$settings
  n <- x$noec
  cat(sprintf("<meso_noec> from %s (direction: %s, alpha = %s)\n",
              if (s$method == "williams") "Williams' test" else "Dunnett's test", s$direction,
              format(s$alpha)))
  cat(sprintf("Criteria: %s\n", if (is.null(s$criteria)) "none (significant effect only)" else
    if (identical(s$criteria$source, "user")) "given by the user" else s$criteria$source))
  st <- table(factor(n$status, levels = unique(c("ok", sort(unique(n$status))))))
  cat("Status: ", paste(sprintf("%s %d", names(st), as.integer(st)), collapse = ", "), "\n", sep = "")
  cat(sprintf("Effects detected: %d; selected: %d\n", sum(n$crit_a %in% TRUE), sum(n$selected)))
  if (any(n$selected)) {
    cat("Selected NOECs (concentration unit of the data):\n")
    print(noec_table(x, selected = TRUE), row.names = FALSE)
  }
  invisible(x)
}

#' @export
as.data.frame.meso_noec <- function(x, ...) {
  x$noec
}
