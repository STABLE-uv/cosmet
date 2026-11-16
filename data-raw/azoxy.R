# Data of the package: `azoxy` (counts of the azoxystrobin mesocosm experiment)
# and `s11_thresholds` (criteria of Supplementary Text S11 of the same study).
#
# Source of the counts: inst/extdata/azoxystrobin/ of the R package bmassa
# (github.com/STABLE-uv/bmassa, commit 05a7ee8, 2026-11-10), written by its
# data-raw/azoxystrobin.R from the data workbook of
#   Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#   Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#   fungicide azoxystrobin in two habitats representative of Mediterranean
#   coastal wetlands: a mesocosm experiment. Aquatic Toxicology 267, 106828.
#   https://doi.org/10.1016/j.aquatox.2023.106828
# Licence: CC BY 4.0 (the authors' choice); cite the article.
# Copied to data-raw/azoxystrobin/ on 2026-11-16 (tools::md5sum):
#   design.csv              a17ab15d2ce215ef28862267a21189d5
#   phytoplankton.csv       8e059a7f3bbfb599b3a0cabf70f88b4a
#   zooplankton.csv         3ad7ad298158dd06036bf8a94b2a50aa
#   macroinvertebrates.csv  3a55b805a38c3407e0d56b681503a304
#
# Units: phytoplankton cells per millilitre (not per litre, as the
# documentation of bmassa 0.6.1 says), zooplankton individuals per litre,
# macroinvertebrates individuals in the two traps of each mesocosm and date
# (one sample). The tables have no tank identifier: within a day, habitat and
# concentration, the three rows are the three replicates, always in the same
# order (data-raw/azoxystrobin.R of bmassa), numbered 1 to 3 by design.csv.
# Taxon names are as in the workbook (with its spelling), as in bmassa.
#
# Thresholds of Text S11 (Supplementary Information of Amador et al. 2024), as
# used in the analysis script of the article; the factors A of ln(Ax + 1) are
# those of the Methods of the article.
#
# Packages used only here: usethis. Run from the package root, in a UTF-8
# locale: Rscript data-raw/azoxy.R

raw_dir <- file.path("data-raw", "azoxystrobin")
md5 <- c(design.csv = "a17ab15d2ce215ef28862267a21189d5",
         phytoplankton.csv = "8e059a7f3bbfb599b3a0cabf70f88b4a",
         zooplankton.csv = "3ad7ad298158dd06036bf8a94b2a50aa",
         macroinvertebrates.csv = "3a55b805a38c3407e0d56b681503a304")
stopifnot(identical(unname(tools::md5sum(file.path(raw_dir, names(md5)))), unname(md5)))

design <- read.csv(file.path(raw_dir, "design.csv"), stringsAsFactors = FALSE)

to_long <- function(file, community, unit) {
  d <- read.csv(file.path(raw_dir, file), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(names(d)[1:3], c("Day", "Treatment", "Concentration")), !anyNA(d))
  replicate <- stats::ave(seq_len(nrow(d)), d$Day, d$Treatment, d$Concentration,
                          FUN = seq_along)
  tank <- match(paste(d$Treatment, d$Concentration, replicate),
                paste(design$treatment, design$concentration_ug_L, design$replicate))
  stopifnot(!anyNA(tank))
  taxa <- setdiff(names(d), c("Day", "Treatment", "Concentration"))
  k <- length(taxa)
  data.frame(
    habitat = rep(design$habitat[tank], k),
    tank = rep(design$tank[tank], k),
    concentration = rep(d$Concentration, k),
    day = rep(as.numeric(d$Day), k),
    community = community,
    taxon = rep(taxa, each = nrow(d)),
    abundance = as.numeric(unlist(d[taxa], use.names = FALSE)),
    unit = unit,
    stringsAsFactors = FALSE
  )
}

azoxy <- rbind(
  to_long("phytoplankton.csv", "phytoplankton", "cells/mL"),
  to_long("zooplankton.csv", "zooplankton", "individuals/L"),
  to_long("macroinvertebrates.csv", "macroinvertebrates", "individuals/sample")
)
azoxy <- azoxy[order(azoxy$community, azoxy$habitat, azoxy$tank, azoxy$day, azoxy$taxon,
                     method = "radix"), ]
rownames(azoxy) <- NULL
stopifnot(!anyDuplicated(azoxy[c("community", "tank", "day", "taxon")]),
          nrow(azoxy) == 120 * (45 + 34 + 29))

s11_thresholds <- data.frame(
  community = c("phytoplankton_large", "phytoplankton_small", "microzooplankton",
                "macrozooplankton", "macroinvertebrates"),
  A = c(33.3, 33.3, 10, 10, 2),
  min_abundance = c(10, 100, 10, 3, 3),
  unit = c("cells/mL", "cells/mL", "individuals/L", "individuals/L", "individuals/sample"),
  stringsAsFactors = FALSE
)

usethis::use_data(azoxy, s11_thresholds, overwrite = TRUE, compress = "xz")
