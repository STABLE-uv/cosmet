#' Counts of the azoxystrobin mesocosm experiment
#'
#' Abundances of phytoplankton, zooplankton and macroinvertebrates in the 24
#' outdoor mesocosms (1230 L) of Amador et al. (2024), with the communities of
#' two habitats of the Albufera Natural Park (Valencia, Spain), a eutrophic
#' lake and a water spring, exposed to the fungicide azoxystrobin: controls
#' and 2, 20 and 200 µg/L, applied twice, on days 0 and 14, three mesocosms
#' each. Samples were taken 7 days before the first application and 7, 21, 35
#' and 56 days after it. Each habitat is a separate experiment.
#'
#' @format A data frame with 12960 rows (one per community, mesocosm, day and
#'   taxon) and 8 columns:
#' \describe{
#'   \item{habitat}{`"eutrophic_lake"` or `"water_spring"` (character).}
#'   \item{tank}{Mesocosm identifier: habitat code (`W` lake, `S` spring),
#'     concentration and replicate (character).}
#'   \item{concentration}{Nominal concentration of azoxystrobin (µg/L).}
#'   \item{day}{Day relative to the first application.}
#'   \item{community}{`"phytoplankton"`, `"zooplankton"` or
#'     `"macroinvertebrates"`.}
#'   \item{taxon}{Taxon as recorded in the data workbook, with its spelling
#'     (character).}
#'   \item{abundance}{Abundance in `unit`. A zero is a taxon looked for and
#'     not found; there are no unknown values.}
#'   \item{unit}{`"cells/mL"` (phytoplankton), `"individuals/L"`
#'     (zooplankton) or `"individuals/sample"` (macroinvertebrates: the
#'     individuals caught in the two traps, 30 x 10 cm each, of a mesocosm
#'     and date).}
#' }
#' @source The data workbook of the authors, as distributed in
#'   `inst/extdata/azoxystrobin` of the R package 'bmassa'; see
#'   `data-raw/azoxy.R` in the source of 'cosmet'. Licence CC BY 4.0: cite
#'   the article when you use them.
#' @references
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @examples
#' table(azoxy$community, azoxy$habitat)
#' head(azoxy[azoxy$abundance > 0, ])
"azoxy"

#' Thresholds of the criteria of Text S11
#'
#' Abundance thresholds of criterion b, by community, of the criteria for
#' consistent, treatment-related population effects of Supplementary Text
#' S11 of Amador et al. (2024), with the factor A of ln(A x + 1) that the
#' study used for each community. [s11_criteria()] reads this table.
#'
#' @format A data frame with 5 rows and 4 columns:
#' \describe{
#'   \item{community}{`"phytoplankton_large"`, `"phytoplankton_small"`,
#'     `"microzooplankton"`, `"macrozooplankton"` or
#'     `"macroinvertebrates"`.}
#'   \item{A}{Factor of ln(A x + 1): 33.3 for phytoplankton, 10 for
#'     zooplankton and 2 for macroinvertebrates (Methods of the article).}
#'   \item{min_abundance}{Minimum mean abundance of criterion b.}
#'   \item{unit}{Unit of `min_abundance`.}
#' }
#' @source Supplementary Text S11 and Methods of Amador et al. (2024); see
#'   `data-raw/azoxy.R`.
#' @references
#' Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#' Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#' fungicide azoxystrobin in two habitats representative of Mediterranean
#' coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267, 106828.
#' \doi{10.1016/j.aquatox.2023.106828}
#' @examples
#' s11_thresholds
"s11_thresholds"
