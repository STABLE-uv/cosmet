# Measured azoxystrobin concentrations of the mesocosm experiment of
# Amador et al. (2024), the data set `azoxy_exposure`.
#
# Source: data-raw/azoxystrobin/exposure.csv, written from the original
# laboratory table of the experiment (2022): one row per mesocosm and sampling
# occasion, 126 samples. The water was sampled on days 0, 7, 14 (before and
# after the second application), 21 and 56 after the first application; the
# controls only on days 0, 14 and 56. Samples not detected ("not_detected"),
# reported below 0.1 ug/L ("below_0.1") or lost ("lost") have no value. The
# flag "repeated_analysis" marks samples analysed twice, with the value of the
# repeated analysis. See inst/COPYRIGHTS for the licence.
#
# Packages used only here: usethis. Run from the package root, in a UTF-8
# locale: Rscript data-raw/azoxy_exposure.R

src <- file.path("data-raw", "azoxystrobin", "exposure.csv")
stopifnot(unname(tools::md5sum(src)) == "ba32cd5c5476044a499b2d26cf8d2c97")

raw <- utils::read.csv(src, stringsAsFactors = FALSE, na.strings = "")
azoxy_exposure <- data.frame(
  habitat = raw$habitat,
  mesocosm = raw$mesocosm,
  concentration = as.numeric(raw$concentration_ug_L),
  day = as.numeric(raw$day),
  timing = raw$timing,
  azoxystrobin = raw$azoxystrobin_ug_L,
  flag = raw$flag,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(azoxy_exposure) == 126L,
  !anyNA(azoxy_exposure[c("habitat", "mesocosm", "concentration", "day")]),
  all(azoxy_exposure$timing %in% c(NA, "before", "after")),
  all(is.na(azoxy_exposure$azoxystrobin) ==
        azoxy_exposure$flag %in% c("not_detected", "below_0.1", "lost"))
)
usethis::use_data(azoxy_exposure, overwrite = TRUE, compress = "xz")
