
<!-- README.md is generated from README.Rmd. Edit README.Rmd and render it with devtools::build_readme() or rmarkdown::render("README.Rmd") after installing the package. -->

# *cosmet*

<!-- badges: start -->

[![R-CMD-check](https://github.com/STABLE-uv/cosmet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/STABLE-uv/cosmet/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

*cosmet* analyses aquatic mesocosm and microcosm experiments with
chemicals, the classic way: it validates the design, compares every
treatment with the control for each taxon and sampling date (Williams’
or Dunnett’s test), derives no observed effect concentrations (NOEC)
with optional criteria for consistent effects, computes minimum
detectable differences (MDD), and fits and tests principal response
curves (PRC). Results are data frames that keep every intermediate
decision (status of each taxon and date, criteria, critical values), and
plots follow the usual figures of mesocosm studies.

## Installation

From GitHub, with its dependencies:

``` r
# install.packages("remotes")
remotes::install_github("STABLE-uv/cosmet", dependencies = TRUE)
```

## Example

The macroinvertebrates of the eutrophic-lake mesocosms of the
azoxystrobin experiment of Amador et al. (2024), included in the
package:

``` r
library(cosmet)

macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
d <- meso_data(
  macro,
  mesocosm_col = "tank", time_col = "day", treatment_col = "concentration",
  concentration_col = "concentration", taxon_col = "taxon", value_col = "abundance",
  control = "0", zero_fill = FALSE, unit = "individuals/sample",
  concentration_unit = "µg/L", applications = c(0, 14)
)
d <- transform_response(d, A = 2)

# Population NOECs: Williams' test and the criteria of Text S11
tw <- test_control(d, method = "williams")
noec <- derive_noec(tw, pre_times = -7, criteria = s11_criteria("macroinvertebrates"))
noec_table(noec, selected = TRUE)
#>                       taxon   -7    7   21   35   56
#> 1               Chironomini <NA> <NA> <NA> <NA>  20↓
#> 2 Colymbetes fuscus (adult) <NA> <NA> <NA>   2↓ <NA>
#> 3  Colymbetes fuscus larvae <NA> <NA> <NA> <NA>   2↓
#> 4                Micronecta <NA> <NA>  20↓ <NA> <NA>

# Minimum detectable differences at the highest concentration, day 56
mdd <- calc_mdd(tw)
head(mdd[mdd$concentration == 200 & mdd$time == 56 & !is.na(mdd$mdd_pct),
         c("taxon", "effect_pct", "mdd_pct", "mdd_class")])
#>                         taxon effect_pct   mdd_pct mdd_class
#> 15                 Anisoptera  -80.52679 126.79815         0
#> 30                   Caenidae  -85.43269  82.65584        II
#> 60                Chironomini  -87.04420  85.51850        II
#> 75                     Cloeon -100.00000 112.59202         0
#> 90  Colymbetes fuscus (adult) -100.00000  88.43059        II
#> 105  Colymbetes fuscus larvae -100.00000  84.68090        II

# Principal response curves and their test
prc <- fit_prc(d)
prc
#> <meso_prc> principal response curves
#> 60 samples, 29 taxa; scale: ln(2 x + 1); scaling: symmetric
#> Variance explained by time: 26.9%, by treatment: 27.7%
#> First axis: 34.7% of the treatment variance (9.6% of the total)
#> Taxa with the largest weights:
#>                      taxon     bk
#>                 Micronecta  1.717
#>                   Caenidae  1.466
#>             Orthocladineae  1.340
#>                Chironomini  1.093
#>   Colymbetes fuscus larvae  0.716
#>                     Cloeon  0.586
#>  Colymbetes fuscus (adult)  0.511
#>       Hydroglyphus geminus -0.388
test_prc(prc, n_perm = 499)
#> <meso_prc_test>
#> First PRC axis: F = 8.473, p = 0.002 (499 permutations of whole time series)
#> Tests by time (RDA on treatment):
#>  time  n r_squared adj_r_squared    F p_value n_perm
#>    -7 12     0.346        0.1006 1.41   0.104    499
#>     7 12     0.306        0.0455 1.17   0.268    499
#>    21 12     0.363        0.1246 1.52   0.082    499
#>    35 12     0.455        0.2508 2.23   0.004    499
#>    56 12     0.425        0.2090 1.97   0.004    499
```

``` r
plot_prc(prc, xlab = "Day")
```

<img src="man/figures/README-prc-1.png" width="100%" />

## Methods and decisions

| Step      | Function                                      | Method                                                                                                                                               | Decisions (argument)                                                                                                              |
|-----------|-----------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| Data      | `meso_data()`                                 | One row per mesocosm, time and taxon; observed zeros, unknown values and samples not taken are kept apart                                            | Whether a taxon without a row is a zero (`zero_fill`, no default); analysed samples (`samples`)                                   |
| Scale     | `transform_response()`                        | ln(A x + 1) (van den Brink et al. 2000)                                                                                                              | `A` (no default)                                                                                                                  |
| Tests     | `test_control()`                              | Step-down Williams test (Williams 1971, 1972) with a simulated null distribution, or Dunnett’s test (Dunnett 1955) by numerical integration          | `method` (no default), `direction`, `alpha`, `step_variance` (Williams: variance of all groups or of the reduced design), `na_rm` |
| NOEC      | `derive_noec()`                               | LOEC and NOEC per taxon and date; criteria of consistent effects                                                                                     | `criteria` (`noec_criteria()`, `s11_criteria()`), `pre_times`                                                                     |
| MDD       | `calc_mdd()`                                  | Brock et al. (2015), eq. 2, 5 and 6, and their MDD classes                                                                                           | —                                                                                                                                 |
| Community | `fit_prc()`, `test_prc()`, `community_noec()` | PRC (van den Brink and ter Braak 1999) with ‘vegan’; permutations of whole time series; Williams’ test on the first principal component of each date | `scaling`, `orientation`, `n_perm`, `seed`                                                                                        |

The vignette (`vignette("cosmet")`) goes through the whole analysis and
explains each decision; the help pages give the equations and
references.

## Status

- The NOECs, LOECs, criteria and selections of
  `test_control(step_variance =   "reduced")` and `derive_noec()` are
  identical to those of the original analysis script of the authors on
  its simulated data and on the azoxystrobin macroinvertebrates of both
  habitats; the PRC coefficients, weights and permutation tests are
  identical to those of the original PRC script.
- The simulated null distribution of Williams’ statistic agrees with the
  tables of Williams (1972) within Monte Carlo error; Dunnett’s
  probabilities agree with ‘mvtnorm’ and with nested numerical
  integration.
- `R CMD check` runs on GitHub Actions on Linux, macOS and Windows
  (badge above). The minimum R version (3.6.0) has not been tested.

## Data

`azoxy` comes from Amador et al. (2024) under a CC BY 4.0 licence: cite
the article when you use it.

Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
fungicide azoxystrobin in two habitats representative of Mediterranean
coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267,
106828. <https://doi.org/10.1016/j.aquatox.2023.106828>
