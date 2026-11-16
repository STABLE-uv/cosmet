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
