# Colours: the control in grey and the treatments on the viridis "mako"
# scale (perceptually uniform, readable with colour-vision deficiencies),
# lighter for lower concentrations. The lightest step (0.75 on the scale)
# still has a contrast of about 2:1 with a white background.
treatment_colours <- function(levels_trt, colours = NULL, call = sys.call(-1)) {
  n <- length(levels_trt)
  if (!is.null(colours)) {
    if (!is.character(colours) || anyNA(colours)) {
      abort_input("`colours` must be a character vector of colours.", call = call)
    }
    if (!is.null(names(colours))) {
      missing_trt <- setdiff(levels_trt, names(colours))
      if (length(missing_trt) > 0L) {
        abort_input(sprintf("`colours` has no colour for the treatments %s.",
                            format_items(quote_items(missing_trt))), call = call)
      }
      return(colours[levels_trt])
    }
    if (length(colours) != n) {
      abort_input(sprintf("`colours` must have one colour per treatment (%d), control first.", n),
                  call = call)
    }
    return(stats::setNames(colours, levels_trt))
  }
  k <- n - 1L
  steps <- if (k == 1L) {
    viridisLite::mako(3, begin = 0.15, end = 0.75)[2L]
  } else {
    viridisLite::mako(k, begin = 0.15, end = 0.75, direction = -1)
  }
  stats::setNames(c("grey45", steps), levels_trt)
}

# Sequential fills for ordered classes, darkest first, and the text colour
# that reads on each.
class_fills <- function(k) {
  if (k == 1L) viridisLite::mako(3, begin = 0.15, end = 0.8)[2L] else
    viridisLite::mako(k, begin = 0.15, end = 0.8)
}

text_on <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  lin <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  lum <- 0.2126 * lin[1L, ] + 0.7152 * lin[2L, ] + 0.0722 * lin[3L, ]
  ifelse(lum < 0.3, "#ffffff", "#0b0b0b")
}

treatment_labels <- function(levels_trt, conc, control, conc_unit, labels = NULL,
                             call = sys.call(-1)) {
  if (!is.null(labels)) {
    if (!is.character(labels) || anyNA(labels)) {
      abort_input("`labels` must be a character vector.", call = call)
    }
    if (!is.null(names(labels))) {
      missing_trt <- setdiff(levels_trt, names(labels))
      if (length(missing_trt) > 0L) {
        abort_input(sprintf("`labels` has no label for the treatments %s.",
                            format_items(quote_items(missing_trt))), call = call)
      }
      return(labels[levels_trt])
    }
    if (length(labels) != length(levels_trt)) {
      abort_input(sprintf("`labels` must have one label per treatment (%d), control first.",
                          length(levels_trt)), call = call)
    }
    return(stats::setNames(labels, levels_trt))
  }
  out <- levels_trt
  if (!all(is.na(conc))) {
    out <- paste0(format_conc(conc), if (is.null(conc_unit)) "" else paste0(" ", conc_unit))
  }
  out[levels_trt == control] <- "Control"
  stats::setNames(out, levels_trt)
}

#' Plot principal response curves
#'
#' `plot_prc()` draws the PRC of [fit_prc()]: the coefficient \eqn{c_{dt}} of
#' every treatment over time, and, in a panel on its right, the weights
#' \eqn{b_k} of the taxa with the largest absolute weights. Needs the optional
#' packages 'ggplot2' and 'patchwork'.
#'
#' @param x A `meso_prc` object.
#' @param min_weight Taxa with an absolute weight of at least this value are
#'   labelled (default 0.25).
#' @param weight_scale `"pseudo_log"` (default) spaces the weights with
#'   \eqn{\mathrm{asinh}(b_k / \sigma)}{asinh(bk / sigma)}, which keeps a few
#'   large weights from squeezing the others; `"linear"` uses the weights as
#'   they are. The axis is labelled with the weights in both cases.
#' @param sigma Scale of the pseudo-logarithm (default 0.35); smaller values
#'   compress more.
#' @param colours,labels `NULL` (defaults: the control in grey and the
#'   treatments on the viridis "mako" scale, from light to dark as the
#'   concentration increases; labels with the concentrations) or one colour or
#'   label per treatment, control first, optionally named by treatment.
#' @param xlab Label of the time axis (default `"Time"`).
#' @param label_size Size of the taxon labels (default 3).
#'
#' @return A 'patchwork' object (print it to draw it).
#' @seealso [fit_prc()].
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance",
#'                concentration_unit = "µg/L")
#' prc <- fit_prc(transform_response(d, A = 2))
#' plot_prc(prc, xlab = "Day")
plot_prc <- function(x, min_weight = 0.25, weight_scale = "pseudo_log", sigma = 0.35,
                     colours = NULL, labels = NULL, xlab = "Time", label_size = 3) {
  if (!inherits(x, "meso_prc")) abort_input("`x` must be the result of fit_prc().")
  check_installed(c("ggplot2", "patchwork", "viridisLite"), "plot principal response curves")
  check_number(min_weight, lower = 0)
  check_choice(weight_scale, c("pseudo_log", "linear"))
  check_number(sigma, lower = 0, include_lower = FALSE)
  check_string(xlab)
  check_number(label_size, lower = 0, include_lower = FALSE)
  s <- x$settings
  lev <- s$treatments
  col <- treatment_colours(lev, colours)
  lab <- treatment_labels(lev, s$concentrations, s$control, s$concentration_unit, labels)
  curves <- x$curves
  curves$treatment <- factor(curves$treatment, levels = lev)
  times <- sort(unique(curves$time))

  p_prc <- ggplot2::ggplot(curves, ggplot2::aes(x = .data$time, y = .data$cdt,
                                                colour = .data$treatment,
                                                group = .data$treatment)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70",
                        linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = col, breaks = lev, labels = lab, name = NULL,
                                 drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = times) +
    ggplot2::labs(x = xlab, y = expression("Community response (" * c[dt] * ")")) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "top",
                   axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
                   plot.margin = ggplot2::margin(t = 8, r = 5, b = 8, l = 8))

  w <- x$weights
  w <- w[!is.na(w$bk) & abs(w$bk) >= min_weight, , drop = FALSE]
  if (nrow(w) == 0L) {
    abort_input(sprintf("No taxon has an absolute weight of at least %s; lower `min_weight`.",
                        format(min_weight)))
  }
  w <- w[order(w$bk), , drop = FALSE]
  tf <- if (weight_scale == "pseudo_log") function(b) asinh(b / sigma) else function(b) b
  w$y <- tf(w$bk)
  w$label <- sprintf("%s  (%+.2f)", w$taxon, w$bk)
  breaks_raw <- pretty(range(c(w$bk, 0)), n = 6)
  label_x <- 0.10
  p_bk <- ggplot2::ggplot(w, ggplot2::aes(x = 0, y = .data$y)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = tf(0), linetype = "dashed", colour = "grey70",
                        linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = label_x - 0.006, y = .data$y,
                                       yend = .data$y),
                          colour = "grey65", linewidth = 0.35) +
    ggplot2::geom_point(size = 2.4, colour = "black") +
    ggplot2::geom_text(ggplot2::aes(x = label_x, y = .data$y, label = .data$label),
                       hjust = 0, vjust = 0.5, size = label_size) +
    ggplot2::scale_x_continuous(limits = c(-0.02, 0.18), breaks = NULL,
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(breaks = tf(breaks_raw), labels = round(breaks_raw, 2),
                                expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
    ggplot2::labs(x = NULL, y = expression("Taxon weight (" * b[k] * ")")) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 90,
                                                        margin = ggplot2::margin(r = 8)),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.margin = ggplot2::margin(t = 8, r = 100, b = 8, l = 4))

  patchwork::wrap_plots(A = patchwork::guide_area(), B = p_prc, C = p_bk,
                        design = "AA\nBC", widths = c(5, 1), heights = c(0.12, 1),
                        guides = "collect") &
    ggplot2::theme(legend.position = "top", legend.justification = "center",
                   legend.box.just = "center")
}

#' Plot the NOECs of every taxon and time
#'
#' `plot_noec()` draws the NOECs of [derive_noec()] or [community_noec()] as
#' a grid of taxa by sampling times. Cells with an effect are coloured by
#' their NOEC on the viridis "mako" scale (darker for lower concentrations)
#' and labelled with it, with a
#' triangle for the direction (up for increases, down for decreases); cells
#' without an effect at the highest concentration are light grey; cells that
#' could not be tested are blank. Needs the optional package 'ggplot2'.
#'
#' @param x A `meso_noec` object.
#' @param selected If `TRUE`, only the selected effects are coloured and
#'   labelled, and only the taxa with at least one are shown (default
#'   `FALSE`).
#' @param xlab Label of the time axis (default `"Time"`).
#'
#' @return A ggplot object.
#' @seealso [noec_table()].
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' noec <- derive_noec(test_control(transform_response(d, A = 2), method = "dunnett"),
#'                     pre_times = -7)
#' plot_noec(noec, xlab = "Day")
plot_noec <- function(x, selected = FALSE, xlab = "Time") {
  if (!inherits(x, "meso_noec")) {
    abort_input("`x` must be the result of derive_noec() or community_noec().")
  }
  check_installed(c("ggplot2", "viridisLite"), "plot NOECs")
  check_flag(selected)
  check_string(xlab)
  n <- x$noec
  conc <- sort(x$settings$concentrations[-1L])
  has_effect <- !is.na(n$loec)
  no_effect <- n$crit_a %in% FALSE
  if (selected) {
    has_effect <- has_effect & n$selected
    no_effect <- rep(FALSE, nrow(n))
    n <- n[n$taxon %in% n$taxon[n$selected], , drop = FALSE]
    has_effect <- has_effect[x$noec$taxon %in% n$taxon]
    no_effect <- no_effect[x$noec$taxon %in% n$taxon]
    if (nrow(n) == 0L) {
      abort_input("No effect was selected; use selected = FALSE to see all NOECs.")
    }
  }
  n$label <- ifelse(has_effect | no_effect, n$noec_label, NA_character_)
  # Direction as a triangle (drawn, so it does not depend on the fonts of the
  # graphics device): up for increases, down for decreases.
  n$shape <- NA_real_
  n$shape[has_effect & n$direction %in% "increase"] <- 24
  n$shape[has_effect & n$direction %in% "decrease"] <- 25
  # NOEC classes: below the lowest concentration, each tested concentration
  # below the highest, and no effect.
  lev_effect <- c(paste0("<", format_conc(conc[1L])), format_conc(conc[-length(conc)]))
  lev <- c(lev_effect, "No effect")
  n$fill <- NA_character_
  n$fill[has_effect] <- n$noec_label[has_effect]
  n$fill[no_effect] <- "No effect"
  n$fill <- factor(n$fill, levels = lev)
  # Lower NOEC (more sensitive) = darker; no effect = light grey.
  fills <- stats::setNames(c(class_fills(length(lev_effect)), "#f0efec"), lev)
  text_col <- stats::setNames(text_on(fills), lev)
  n$text_colour <- ifelse(is.na(n$fill), "#0b0b0b", text_col[as.character(n$fill)])
  n$taxon <- factor(n$taxon, levels = rev(unique(n$taxon)))
  n$time_f <- factor(format(n$time, trim = TRUE),
                     levels = format(sort(unique(n$time)), trim = TRUE))
  n$x <- as.integer(n$time_f)
  n$label_x <- n$x - ifelse(is.na(n$shape), 0, 0.1)
  arrows <- n[!is.na(n$shape), , drop = FALSE]
  arrows$x <- arrows$x + 0.28
  ggplot2::ggplot(n, ggplot2::aes(x = .data$x, y = .data$taxon)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "#ffffff", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(x = .data$label_x, label = .data$label,
                                    colour = .data$text_colour),
                       size = 3, na.rm = TRUE) +
    ggplot2::geom_point(data = arrows, ggplot2::aes(x = .data$x, shape = .data$shape,
                                                    colour = .data$text_colour),
                        fill = arrows$text_colour, size = 1.6) +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(n$time_f)), labels = levels(n$time_f),
                                expand = ggplot2::expansion(add = 0.05)) +
    ggplot2::scale_fill_manual(values = fills, breaks = lev, drop = FALSE,
                               na.value = "#ffffff", name = "NOEC") +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_shape_identity() +
    ggplot2::labs(x = xlab, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "top")
}

#' Plot the minimum detectable differences
#'
#' `plot_mdd()` draws the MDD classes of [calc_mdd()] for one concentration
#' as a grid of taxa by sampling times, labelled with the MDD as % of the
#' control. Needs the optional package 'ggplot2'.
#'
#' @param x A data frame from [calc_mdd()].
#' @param concentration The concentration to show; `NULL` (default) for the
#'   highest.
#' @param xlab Label of the time axis (default `"Time"`).
#'
#' @return A ggplot object.
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance")
#' mdd <- calc_mdd(test_control(transform_response(d, A = 2), method = "dunnett"))
#' plot_mdd(mdd, xlab = "Day")
plot_mdd <- function(x, concentration = NULL, xlab = "Time") {
  check_data_frame(x)
  check_columns(x, c("taxon", "time", "concentration", "mdd_pct", "mdd_class"))
  check_installed(c("ggplot2", "viridisLite"), "plot minimum detectable differences")
  check_number(concentration, allow_null = TRUE)
  check_string(xlab)
  conc_all <- sort(unique(x$concentration))
  if (is.null(concentration)) concentration <- max(conc_all)
  if (!concentration %in% conc_all) {
    abort_input(sprintf("`concentration` must be one of %s.", format_items(conc_all, max = 20L)))
  }
  d <- x[x$concentration == concentration, , drop = FALSE]
  lev <- c("IV", "III", "II", "I", "0")
  fills <- stats::setNames(c(class_fills(4L), "#f0efec"), lev)
  text_col <- stats::setNames(text_on(fills), lev)
  d$fill <- factor(d$mdd_class, levels = lev)
  d$label <- ifelse(is.na(d$mdd_pct), NA_character_, sprintf("%.0f", d$mdd_pct))
  d$text_colour <- ifelse(is.na(d$mdd_class), "#0b0b0b", text_col[d$mdd_class])
  d$taxon <- factor(d$taxon, levels = rev(sort(unique(x$taxon), method = "radix")))
  d$time_f <- factor(format(d$time, trim = TRUE),
                     levels = format(sort(unique(x$time)), trim = TRUE))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$time_f, y = .data$taxon)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "#ffffff", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$text_colour),
                       size = 3, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = fills, breaks = lev, drop = FALSE,
                               na.value = "#ffffff",
                               labels = c("IV (< 50%)", "III (50-70%)", "II (70-90%)",
                                          "I (90-100%)", "0 (> 100%)"),
                               name = "MDD class") +
    ggplot2::scale_colour_identity() +
    ggplot2::labs(x = xlab, y = NULL,
                  subtitle = sprintf("MDD (%% of the control) at concentration %s",
                                     format_conc(concentration))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "top")
}

#' Plot taxa over time
#'
#' `plot_timecourse()` draws, for each taxon, the mean of every treatment at
#' every sampling time, with the application times as dashed vertical lines.
#' Needs the optional package 'ggplot2'.
#'
#' @param x A [meso_data()] object.
#' @param taxa The taxa (or variables) to plot.
#' @param scale `"response"` (default, the scale of the analyses) or
#'   `"value"` (the original values).
#' @param show_replicates If `TRUE`, also the value of each mesocosm (default
#'   `FALSE`).
#' @param colours,labels As in [plot_prc()].
#' @param xlab Label of the time axis (default `"Time"`).
#'
#' @return A ggplot object.
#'
#' @section Missing values and zeros:
#' Means are those of the known values; unknown values are left out of them.
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
#' d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
#'                time_col = "day", treatment_col = "concentration",
#'                concentration_col = "concentration", value_col = "abundance",
#'                applications = c(0, 14))
#' plot_timecourse(transform_response(d, A = 2), taxa = c("Chironomini", "Cloeon"))
plot_timecourse <- function(x, taxa, scale = "response", show_replicates = FALSE,
                            colours = NULL, labels = NULL, xlab = "Time") {
  check_meso_data(x)
  check_installed(c("ggplot2", "viridisLite"), "plot taxa over time")
  if (missing(taxa)) abort_input("`taxa` is missing: give the taxa to plot.")
  check_choice(scale, c("response", "value"))
  check_flag(show_replicates)
  check_string(xlab)
  d <- select_taxa(x, taxa)
  s <- x$settings
  lev <- s$treatments
  conc <- x$design$concentration[match(lev, x$design$treatment)]
  col <- treatment_colours(lev, colours)
  lab <- treatment_labels(lev, conc, s$control, s$concentration_unit, labels)
  d$y <- d[[scale]]
  keep <- !is.na(d$y)
  grp <- interaction(d$taxon[keep], d$treatment[keep], format_time_key(d$time[keep]),
                     drop = TRUE, sep = "\u001f")
  means <- data.frame(
    taxon = tapply(d$taxon[keep], grp, function(v) v[1L]),
    treatment = tapply(d$treatment[keep], grp, function(v) v[1L]),
    time = tapply(d$time[keep], grp, function(v) v[1L]),
    y = tapply(d$y[keep], grp, mean),
    stringsAsFactors = FALSE
  )
  means$treatment <- factor(means$treatment, levels = lev)
  d$treatment <- factor(d$treatment, levels = lev)
  ylab <- if (scale == "value") {
    if (is.null(s$unit)) "Value" else paste0("Value (", s$unit, ")")
  } else {
    if (identical(s$transform$method, "log")) {
      sprintf("ln(%s x + 1)", format_conc(s$transform$A))
    } else {
      "Value"
    }
  }
  p <- ggplot2::ggplot(means, ggplot2::aes(x = .data$time, y = .data$y,
                                           colour = .data$treatment, group = .data$treatment))
  if (!is.null(s$applications)) {
    p <- p + ggplot2::geom_vline(xintercept = s$applications, linetype = "dashed",
                                 colour = "#898781", linewidth = 0.4)
  }
  if (show_replicates) {
    p <- p + ggplot2::geom_point(data = d[keep, , drop = FALSE],
                                 ggplot2::aes(x = .data$time, y = .data$y,
                                              colour = .data$treatment),
                                 size = 1.2, alpha = 0.45, inherit.aes = FALSE)
  }
  p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$taxon), scales = "free_y") +
    ggplot2::scale_colour_manual(values = col, breaks = lev, labels = lab, name = NULL,
                                 drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$time))) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top", strip.background = ggplot2::element_blank())
}
