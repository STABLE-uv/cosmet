# data-raw

Script and tables that build the data sets of *cosmet*. This folder is
excluded from the package tarball (`.Rbuildignore`).

| File | Content |
|---|---|
| `azoxy.R` | Builds `data/azoxy.rda` (the counts in long format, one row per tank, sampling day and taxon) and `data/s11_thresholds.rda` (the factors of the transformation and the abundance thresholds of Supplementary Text S11 of Amador et al. 2024). Run from the package root, in a UTF-8 locale: `Rscript data-raw/azoxy.R` (R package *usethis*, Wickham et al. 2023) |
| `azoxy_exposure.R` | Builds `data/azoxy_exposure.rda` (the measured azoxystrobin concentrations, one row per mesocosm and sampling occasion) from `azoxystrobin/exposure.csv`. Run from the package root: `Rscript data-raw/azoxy_exposure.R` (R package *usethis*) |
| `azoxystrobin/` | `design.csv` (habitat, treatment, concentration and replicate of each tank) and the counts of `phytoplankton.csv` (cells/mL), `zooplankton.csv` (individuals/L) and `macroinvertebrates.csv` (individuals per sample: the two traps of a tank on one day), as distributed in `inst/extdata/azoxystrobin/` of the R package *bmassa*; and `exposure.csv`, the measured azoxystrobin (µg/L) from the laboratory table of the experiment. The scripts check their MD5 sums |

The counts and concentrations come from the mesocosm experiment of Amador P, Vega C, Navarro
Pacheco NI, Moratalla-López J, Palacios J, Crettaz Minaglia MC, López I, Díaz
M, Rico A (2024). Effects of the fungicide azoxystrobin in two habitats
representative of Mediterranean coastal wetlands: a mesocosm experiment.
*Aquatic Toxicology* 267, 106828.
<https://doi.org/10.1016/j.aquatox.2023.106828>. Licence: CC BY 4.0; cite the
article when you use them.
