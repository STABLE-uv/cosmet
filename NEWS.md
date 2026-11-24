# cosmet 0.1.0.9000

## New features

* Measured exposure. `meso_exposure()` validates the concentrations measured
  in each mesocosm and assigns every sample to the application it follows;
  samples taken on the day of an application are marked as taken before or
  after it (`timing_col`). It takes a single application, several
  applications at the same times in all the treatments, or a schedule for
  each treatment with the nominal concentration of each application
  (`applications` as a data frame).
* `fit_dissipation()` fits first-order dissipation after each application,
  by default as a regression of ln C on the time since the application with
  one rate constant per application, and gives the rate constants with
  their confidence intervals, the DT50 and the DT90. `method = "nonlinear"`
  fits the untransformed concentrations (the recommendation of FOCUS 2014),
  `k_by = "series"` one rate constant for all the applications and
  `by = "mesocosm"` one fit per mesocosm. `predict()` gives the fitted
  curves.
* `calc_twa()` computes time-weighted average concentrations over windows
  that start at the first application of each treatment, from the measured
  concentrations (trapezoidal rule) or from the fitted curves, also as a
  percentage of the nominal concentration (`twa_pct_nominal`).
* `compare_nominal()` compares the concentration measured right after each
  application with its nominal concentration (`pct_nominal`) and, with the
  sample taken right before it, the increase it produced (`pct_added`).
  `fit_dissipation()` gives the fitted initial concentrations as a
  percentage of the nominal one (`c0_pct_nominal`).
* `plot_exposure()` draws the measured concentrations over time, one panel
  per treatment, with its applications, its nominal concentration and the
  fitted curves (needs 'ggplot2').
* `azoxy_exposure`: measured azoxystrobin of the experiment of Amador et al.
  (2024).

## Changes with respect to the original scripts

* The article describes the dissipation rate constants as a regression of
  the ln-transformed concentrations, but the original script averaged the
  negative rates ln(C2 / C1) / (t2 - t1) between consecutive sampling days,
  over both applications. `fit_dissipation()` fits the regression, one rate
  constant per application. In the azoxystrobin data the dissipation slowed
  down after the second application, most of all in the water spring (DT50
  of 17 to 29 days after the first application and 40 to 43 days after the
  second, against 21 to 31 days from the original script).
* The original script placed the sample taken before the second application
  on day 13, with a one-day ramp up to the sample taken after it on day 14.
  `calc_twa()` keeps both on day 14, where the profile jumps. The 7-day
  time-weighted averages are identical; those of 14 to 56 days are 0.4 % to
  3.0 % lower.

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
* The PRC needs 'vegan' 2.6-6 or later (and so R 4.1.0 or later). Since that
  version, 'vegan' scales the PRC coefficients so that the product of weight
  and coefficient is the difference from the control fitted by the first
  axis, as van den Brink and ter Braak (1999) define it. The original script
  ran with 'vegan' 2.6-4, which gave coefficients smaller by a constant factor
  (0.30 for the macroinvertebrates of the eutrophic lake); the taxon weights,
  the variance explained and the permutation tests do not change.
