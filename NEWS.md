# cosmet 0.1.0

First version.

## New features

* `meso_data()` validates a long table of a mesocosm experiment (mesocosm,
  time, treatment, concentration, taxon, value) and keeps observed zeros,
  unknown values and samples not taken apart. `zero_fill` (no default) says
  whether a taxon without a row in an analysed sample is a zero.
  `check_design()` reports replication and missing samples.
* `transform_response()` sets the scale of the analyses, ln(A x + 1), and
  keeps the original values.
* `taxon_eligibility()` summarises abundance and frequency in the controls.
* `test_control()` compares every treatment with the control by taxon and
  time with the step-down Williams test (simulated null distribution;
  `step_variance = "all"`, the procedure of Williams 1971 and 1972, or
  `"reduced"`, that of the original analysis script of the authors) or with
  Dunnett's test (adjusted p-values and simultaneous intervals by numerical
  integration, without randomness).
* `derive_noec()` gives NOECs and LOECs with optional criteria for
  consistent effects (`noec_criteria()`, and `s11_criteria()` with the
  thresholds of Text S11 of Amador et al. 2024, converted between
  compatible units); `noec_table()` lays them out for reports.
* `calc_mdd()` computes minimum detectable differences following Brock et
  al. (2015), with their classes, next to the estimated effect and its
  interval.
* `fit_prc()`, `test_prc()` and `community_noec()` fit and test principal
  response curves and derive community NOECs from the first principal
  component of each time.
* `plot_prc()`, `plot_noec()`, `plot_mdd()` and `plot_timecourse()` (need
  'ggplot2'; `plot_prc()` also 'patchwork').
* `azoxy`: counts of the azoxystrobin mesocosm experiment of Amador et al.
  (2024), and `s11_thresholds`.

## Changes with respect to the original scripts

* Williams' test uses by default the variance and isotonic means of all
  groups (`step_variance = "all"`); the original script re-analysed the
  reduced design at each step, which `step_variance = "reduced"` reproduces
  (identical NOECs, LOECs and criteria on the simulated data of the script
  and on the azoxystrobin macroinvertebrates). In the `azoxy` data, 49 of
  656 NOECs differ between the two (`dev/validation/williams_variants.md`
  in the source).
* The null distribution of Williams' statistic is simulated with the seed
  for every replication pattern, so the result of a taxon no longer depends
  on the other taxa analysed; p-values differ from the original script within
  Monte Carlo error when more than one pattern is simulated.
* Unknown values (`NA`) are no longer left out silently (NOECs) or turned
  into zeros (PRC): see `na_rm` and `zero_fill`.
* The PRC uses an explicit scaling of the taxon weights (`"symmetric"`), which
  does not change with the default of `vegan::summary.prc()`.
