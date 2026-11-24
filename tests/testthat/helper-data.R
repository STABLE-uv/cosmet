# Small experiments used across the tests. Values chosen so that the
# conclusions do not depend on Monte Carlo error.

# Control and three concentrations, three replicates, one time.
# Means: control 11, 1 -> 10.667, 3 -> 7, 10 -> 2 (clear monotonic decrease).
make_monotone <- function(times = 7) {
  one <- data.frame(
    mesocosm = c("C1", "C2", "C3", "A1", "A2", "A3", "B1", "B2", "B3", "D1", "D2", "D3"),
    concentration = rep(c(0, 1, 3, 10), each = 3),
    value = c(10, 12, 11, 9, 13, 10, 6, 7, 8, 2, 1, 3),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, lapply(times, function(t) cbind(one, time = t)))
  out$taxon <- "Cloeon"
  out
}

meso_monotone <- function(times = 7, ...) {
  meso_data(make_monotone(times), control = "0", zero_fill = FALSE,
            treatment_col = "concentration", concentration_col = "concentration", ...)
}

# Two taxa, two times, with a taxon missing from some samples (presence-only).
make_presence <- function() {
  data.frame(
    mesocosm = c("C1", "C1", "C2", "T1", "T2", "C1", "C2", "T1", "T2"),
    time = c(7, 7, 7, 7, 7, 14, 14, 14, 14),
    treatment = c("ctrl", "ctrl", "ctrl", "high", "high", "ctrl", "ctrl", "high", "high"),
    taxon = c("Cloeon", "Caenis", "Cloeon", "Cloeon", "Caenis", "Cloeon", "Cloeon", "Cloeon",
              "Cloeon"),
    value = c(10, 3, 12, 2, 1, 9, 11, 1, 2),
    stringsAsFactors = FALSE
  )
}

azoxy_macro <- function(habitat = "eutrophic_lake", taxa = NULL, A = 2) {
  macro <- azoxy[azoxy$community == "macroinvertebrates" & azoxy$habitat == habitat, ]
  if (!is.null(taxa)) macro <- macro[macro$taxon %in% taxa, ]
  d <- meso_data(macro, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
                 time_col = "day", treatment_col = "concentration",
                 concentration_col = "concentration", value_col = "abundance",
                 unit = "individuals/sample", applications = c(0, 14))
  if (is.null(A)) d else transform_response(d, A = A)
}

# Measured exposure of a small experiment: a control and two treatments, two
# mesocosms each, applications on days 0 and 14, samples on days 0, 7, 14
# (before and after the second application), 21 and 35. The values follow
# first-order dissipation exactly, with rate constants `k` (one per
# application) and initial concentrations `c0` (one per application, times
# the nominal concentration).
make_exposure <- function(k = c(0.1, 0.05), c0 = c(1, 1.2)) {
  occ <- data.frame(time = c(0, 7, 14, 14, 21, 35),
                    timing = c("after", NA, "before", "after", NA, NA),
                    period = c(1, 1, 1, 2, 2, 2),
                    elapsed = c(0, 7, 14, 0, 7, 21))
  mes <- data.frame(mesocosm = c("C1", "C2", "L1", "L2", "H1", "H2"),
                    concentration = c(0, 0, 10, 10, 100, 100))
  out <- merge(mes, occ, by = NULL)
  out$value <- out$concentration * c0[out$period] * exp(-k[out$period] * out$elapsed)
  out$value[out$concentration == 0] <- NA
  out <- out[order(out$mesocosm, out$time, out$period), ]
  out[c("mesocosm", "concentration", "time", "timing", "value")]
}

meso_exposure_small <- function(data = make_exposure(), ...) {
  meso_exposure(data, control = "0", applications = c(0, 14),
                treatment_col = "concentration", concentration_col = "concentration",
                timing_col = "timing", unit = "ug/L", ...)
}

# Means of the two treated mesocosms: 10, 6, 2 (period 1) and 12, 5, 1
# (period 2), on days 0, 7, 14 and 14, 21, 56. T1 has a zero on day 56.
make_twa <- function() {
  data.frame(
    mesocosm = rep(c("C1", "T1", "T2"), each = 6),
    treatment = rep(c("0", "5", "5"), each = 6),
    time = rep(c(0, 7, 14, 14, 21, 56), 3),
    timing = rep(c("after", NA, "before", "after", NA, NA), 3),
    value = c(NA, NA, NA, NA, NA, NA,
              9, 5, 1, 11, 4, 0,
              11, 7, 3, 13, 6, 2),
    stringsAsFactors = FALSE
  )
}

meso_twa <- function(data = make_twa()) {
  meso_exposure(data, control = "0", applications = c(0, 14), timing_col = "timing")
}

azoxy_exposure_lake <- function() {
  expo <- azoxy_exposure[azoxy_exposure$habitat == "eutrophic_lake", ]
  meso_exposure(expo, control = "0", applications = c(0, 14), mesocosm_col = "mesocosm",
                time_col = "day", treatment_col = "concentration",
                concentration_col = "concentration", value_col = "azoxystrobin",
                timing_col = "timing", unit = "ug/L")
}

# A schedule for each treatment: "single" receives 20 on day 0 and "weekly"
# 5 on days 0, 7, 14 and 21 (or from day `first` on). The values follow
# first-order dissipation with rate constant `k` and add up the doses
# (perfect dosing), so after each application the curve is first-order with
# C0 equal to the dose plus what remained. The control has no applications
# but its samples on the days of the weekly applications are labelled before
# and after them, like the others.
make_schedule <- function(k = 0.1, first = 0) {
  weekly <- first + c(0, 7, 14, 21)
  schedule <- data.frame(treatment = c("single", rep("weekly", 4)),
                         time = c(0, weekly), nominal = c(20, 5, 5, 5, 5))
  times <- sort(unique(c(0, 3, 10, 17, 24, 28, weekly)))
  occ <- data.frame(time = c(times, weekly[-1]),
                    timing = c(ifelse(times %in% c(0, weekly), "after", NA),
                               rep("before", 3)))
  mes <- data.frame(mesocosm = c("C1", "S1", "S2", "W1", "W2"),
                    treatment = c("control", "single", "single", "weekly", "weekly"),
                    nominal_conc = c(0, 20, 20, 5, 5), stringsAsFactors = FALSE)
  out <- merge(mes, occ, by = NULL)
  out$value <- vapply(seq_len(nrow(out)), function(i) {
    if (out$treatment[i] == "control") return(NA_real_)
    s <- schedule[schedule$treatment == out$treatment[i], ]
    dosed <- s$time < out$time[i] |
      (s$time == out$time[i] & !(out$timing[i] %in% "before"))
    sum(s$nominal[dosed] * exp(-k * (out$time[i] - s$time[dosed])))
  }, numeric(1))
  out <- out[order(out$mesocosm, out$time, match(out$timing, c("before", "after"))), ]
  rownames(out) <- NULL
  list(data = out, schedule = schedule)
}

meso_schedule <- function(x = make_schedule(), ...) {
  meso_exposure(x$data, control = "control", applications = x$schedule,
                concentration_col = "nominal_conc", timing_col = "timing", unit = "ug/L", ...)
}
