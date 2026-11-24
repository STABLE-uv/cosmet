# meso_exposure() -------------------------------------------------------------

test_that("samples are assigned to the application they follow", {
  e <- meso_exposure_small()
  h <- e$data[e$data$mesocosm == "H1", ]
  expect_identical(h$time, c(0, 7, 14, 14, 21, 35))
  expect_identical(h$timing, c("after", NA, "before", "after", NA, NA))
  expect_identical(h$period, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_identical(h$application, c(0, 0, 0, 14, 14, 14))
  expect_identical(h$elapsed, c(0, 7, 14, 0, 7, 21))
  expect_identical(e$settings$treatments, c("0", "10", "100"))
  expect_false(e$settings$timing_assumed)
  expect_identical(e$diagnostics$n_missing, 12L)
  expect_identical(unique(e$data$value_status[e$data$treatment == "0"]), "missing")
})

test_that("the result does not depend on the order of the rows", {
  dat <- make_exposure()
  set.seed(1)
  shuffled <- dat[sample(nrow(dat)), ]
  expect_identical(meso_exposure_small(shuffled)$data, meso_exposure_small(dat)$data)
})

test_that("tibbles are accepted and zeros are kept", {
  skip_if_not_installed("tibble")
  dat <- make_exposure()
  dat$value[dat$mesocosm == "L1" & dat$time == 35] <- 0
  e <- meso_exposure_small(tibble::as_tibble(dat))
  expect_identical(e$diagnostics$n_zero, 1L)
  expect_identical(e$data$value_status[e$data$mesocosm == "L1" & e$data$time == 35], "zero")
})

test_that("without a timing column, samples at an application count as after it", {
  dat <- make_exposure()
  dat <- dat[!(dat$time == 14 & dat$timing %in% "before"), ]
  e <- meso_exposure(dat, control = "0", applications = c(0, 14),
                     treatment_col = "concentration", concentration_col = "concentration")
  expect_true(e$settings$timing_assumed)
  h <- e$data[e$data$mesocosm == "H1", ]
  expect_identical(h$period, c(1L, 1L, 2L, 2L, 2L))
  expect_identical(h$timing, c("after", NA, "after", NA, NA))
  expect_error(
    meso_exposure(make_exposure(), control = "0", applications = c(0, 14),
                  treatment_col = "concentration"),
    "timing_col"
  )
})

test_that("invalid exposure data give informative errors", {
  dat <- make_exposure()
  expect_error(meso_exposure(dat, applications = 0, treatment_col = "concentration"),
               "`control` is missing", class = "cosmet_error_input")
  expect_error(meso_exposure(dat, control = "0", treatment_col = "concentration"),
               "`applications` is missing", class = "cosmet_error_input")
  expect_error(meso_exposure_small(dat[0, ]), "no rows", class = "cosmet_error_input")
  expect_error(meso_exposure(dat, control = "5", applications = c(0, 14),
                             treatment_col = "concentration", timing_col = "timing"),
               "not a level", class = "cosmet_error_input")
  bad <- dat
  bad$value[3] <- -1
  expect_error(meso_exposure_small(bad), "value", class = "cosmet_error_input")
  bad <- dat
  bad$timing[bad$time == 14 & bad$mesocosm == "H1"] <- NA
  expect_error(meso_exposure_small(bad), "need \"before\" or \"after\"",
               class = "cosmet_error_input")
  bad <- dat
  bad$timing[bad$time == 7] <- "after"
  expect_error(meso_exposure_small(bad), "not taken at an application time",
               class = "cosmet_error_input")
  bad <- dat
  bad$timing[bad$timing %in% "before"] <- "pre"
  expect_error(meso_exposure_small(bad), "'pre'", class = "cosmet_error_input")
  bad <- dat
  bad$time <- as.Date("2022-05-01") + bad$time
  expect_error(meso_exposure_small(bad), "dates", class = "cosmet_error_input")
  expect_error(meso_exposure(dat, control = "0", applications = NA,
                             treatment_col = "concentration"),
               "finite numbers", class = "cosmet_error_input")
  expect_error(meso_exposure_small(rbind(dat, dat[1, ])), "one sample per time",
               class = "cosmet_error_input")
})

test_that("the exposure object prints a summary", {
  expect_snapshot(print(meso_exposure_small()))
})

# fit_dissipation() -----------------------------------------------------------

test_that("first-order data give back their rate constants and initial concentrations", {
  e <- meso_exposure_small()
  for (method in c("log_linear", "nonlinear")) {
    est <- fit_dissipation(e, method = method)$estimates
    expect_identical(est$treatment, c("10", "10", "100", "100"))
    expect_identical(est$period, c(1L, 2L, 1L, 2L))
    expect_equal(est$k, c(0.1, 0.05, 0.1, 0.05), tolerance = 1e-8)
    expect_equal(est$c0, c(10, 12, 100, 120), tolerance = 1e-8)
    expect_equal(est$dt50, log(2) / est$k)
    expect_equal(est$dt90, log(10) / est$k)
    expect_identical(est$n, rep(6L, 4))
    expect_identical(est$df, rep(4L, 4))
    expect_identical(est$status, rep("ok", 4))
  }
  same_k <- meso_exposure_small(make_exposure(k = c(0.08, 0.08)))
  for (method in c("log_linear", "nonlinear")) {
    est <- fit_dissipation(same_k, method = method, k_by = "series")$estimates
    expect_equal(est$k, rep(0.08, 4), tolerance = 1e-8)
    expect_equal(est$c0, c(10, 12, 100, 120), tolerance = 1e-8)
    expect_identical(est$df, rep(9L, 4))
  }
  est <- fit_dissipation(e, by = "mesocosm")$estimates
  expect_identical(est$mesocosm, rep(c("L1", "L2", "H1", "H2"), each = 2))
  expect_equal(est$k, rep(c(0.1, 0.05), 4), tolerance = 1e-8)
})

test_that("the log-linear fit is the regression of ln C on the time since the application", {
  e <- azoxy_exposure_lake()
  dis <- fit_dissipation(e, level = 0.9)
  d <- e$data[e$data$treatment == "200" & e$data$period == 1L, ]
  ref <- stats::lm(log(value) ~ elapsed, data = d)
  est <- dis$estimates[dis$estimates$treatment == "200" & dis$estimates$period == 1L, ]
  expect_equal(est$k, -unname(stats::coef(ref)[2]))
  expect_equal(est$k_se, unname(summary(ref)$coefficients[2, 2]))
  expect_equal(c(est$k_upper, est$k_lower), -unname(stats::confint(ref, level = 0.9)[2, ]))
  expect_equal(est$c0, exp(unname(stats::coef(ref)[1])))
  expect_equal(c(est$dt50_lower, est$dt50_upper), log(2) / c(est$k_upper, est$k_lower))
  expect_identical(est$df, 7L)

  series <- fit_dissipation(e, k_by = "series")$estimates
  d <- e$data[e$data$treatment == "20" & e$data$period > 0L, ]
  ref <- stats::lm(log(value) ~ 0 + factor(period) + elapsed, data = d)
  est <- series[series$treatment == "20", ]
  expect_equal(est$k, -rep(unname(stats::coef(ref)["elapsed"]), 2))
  expect_equal(est$c0, exp(unname(stats::coef(ref)[1:2])))
  expect_identical(est$df, c(15L, 15L))
})

test_that("the nonlinear fit agrees with nls() on the concentrations", {
  e <- azoxy_exposure_lake()
  dis <- fit_dissipation(e, method = "nonlinear")
  d <- e$data[e$data$treatment == "200" & e$data$period == 2L, ]
  ref <- stats::nls(value ~ c0 * exp(-k * elapsed), data = d,
                    start = list(c0 = 200, k = 0.05))
  est <- dis$estimates[dis$estimates$treatment == "200" & dis$estimates$period == 2L, ]
  expect_equal(est$k, unname(stats::coef(ref)["k"]), tolerance = 1e-5)
  expect_equal(est$c0, unname(stats::coef(ref)["c0"]), tolerance = 1e-5)
  expect_equal(est$k_se, unname(summary(ref)$coefficients["k", 2]), tolerance = 1e-4)
})

test_that("estimates that cannot be computed are NA with a reason", {
  dat <- make_exposure()
  # L: a single sampling time in period 2; H2: only two values in period 1
  dat <- dat[!(dat$mesocosm %in% c("L1", "L2") & dat$time %in% c(21, 35)), ]
  dat <- dat[!(dat$mesocosm == "H2" & dat$time == 7), ]
  e <- meso_exposure_small(dat)
  est <- fit_dissipation(e)$estimates
  expect_identical(est$status, c("ok", "not_estimable", "ok", "ok"))
  expect_match(est$reason[2], "fewer than two sampling times")
  expect_true(all(is.na(unlist(est[2, c("k", "c0", "dt50", "k_lower")]))))
  m <- fit_dissipation(e, by = "mesocosm")$estimates
  h2 <- m[m$mesocosm == "H2" & m$period == 1L, ]
  expect_identical(h2$status, "no_uncertainty")
  expect_identical(h2$df, 0L)
  expect_equal(h2$k, 0.1)
  expect_true(is.na(h2$k_lower) && is.na(h2$dt50_upper))
  # Period 2 of L has no curve of its own, but the series fit uses its value
  s <- fit_dissipation(e, k_by = "series")$estimates
  expect_identical(s$status[1:2], c("ok", "ok"))
  expect_equal(s$k[1:2], c(0.1, 0.1))

  rising <- make_exposure(k = c(-0.02, 0.05))
  rising$value[rising$mesocosm %in% c("H1", "L1") & rising$time == 7] <-
    rising$value[rising$mesocosm %in% c("H1", "L1") & rising$time == 7] * 0.9
  est <- fit_dissipation(meso_exposure_small(rising))$estimates
  expect_identical(est$status[est$period == 1L], c("no_decline", "no_decline"))
  expect_true(all(is.na(est$dt50[est$period == 1L])))
  upper <- est$dt50_upper[est$period == 1L]
  expect_true(all(upper == Inf | is.na(upper)))

  wide <- make_exposure()
  wide$value[wide$mesocosm == "H1" & wide$time == 7] <- 80
  est <- fit_dissipation(meso_exposure_small(wide), by = "mesocosm")$estimates
  row <- est[est$mesocosm == "H1" & est$period == 1L, ]
  expect_gt(row$k, 0)
  expect_identical(row$status, "ok")
  expect_lt(row$k_lower, 0)
  expect_identical(row$dt50_upper, Inf)
})

test_that("zeros, missing values and controls are left out of the fits", {
  dat <- make_exposure()
  dat$value[dat$mesocosm == "H1" & dat$time == 21] <- 0
  dat$value[dat$mesocosm == "H2" & dat$time == 7] <- NA
  dis <- fit_dissipation(meso_exposure_small(dat))
  expect_identical(dis$diagnostics$n_zero, 1L)
  expect_identical(dis$diagnostics$n_missing, 1L)
  expect_identical(dis$estimates$n, c(6L, 6L, 5L, 5L))
  expect_equal(dis$estimates$k, c(0.1, 0.05, 0.1, 0.05), tolerance = 1e-8)
  expect_false(any(dis$fitted$treatment == "0"))
  expect_equal(dis$fitted$fitted, dis$fitted$value, tolerance = 1e-8)
})

test_that("fit options are checked", {
  e <- meso_exposure_small()
  expect_error(fit_dissipation(e, method = "log"), "must be one of", class = "cosmet_error_input")
  expect_error(fit_dissipation(e, periods = 3), "from 1 to 2", class = "cosmet_error_input")
  expect_error(fit_dissipation(e, periods = 0), "from 1 to 2", class = "cosmet_error_input")
  expect_error(fit_dissipation(e, level = 1), "level", class = "cosmet_error_input")
  expect_error(fit_dissipation(e$data), "meso_exposure", class = "cosmet_error_input")
  one <- fit_dissipation(e, periods = 2)
  expect_identical(one$estimates$period, c(2L, 2L))
})

test_that("the fit prints a summary", {
  e <- meso_exposure_small()
  expect_snapshot({
    print(fit_dissipation(e))
    print(fit_dissipation(e, method = "nonlinear", k_by = "series", by = "mesocosm"))
  })
})

test_that("predictions follow the curve of each application", {
  dis <- fit_dissipation(meso_exposure_small())
  p <- predict(dis, times = c(-1, 0, 7, 14, 28))
  h <- p[p$treatment == "100", ]
  expect_identical(h$period, c(0L, 1L, 1L, 2L, 2L))
  expect_equal(h$fit, c(NA, 100, 100 * exp(-0.7), 120, 120 * exp(-0.7)), tolerance = 1e-8)
  grid <- predict(dis, n = 5)
  g <- grid[grid$treatment == "10", ]
  expect_identical(g$time, c(seq(0, 14, length.out = 5), seq(14, 35, length.out = 5)))
  expect_error(predict(dis, times = NA), "finite numbers", class = "cosmet_error_input")
})

# calc_twa() ------------------------------------------------------------------

test_that("trapezoidal time-weighted averages match the hand calculation", {
  e <- meso_twa()
  twa <- calc_twa(e, windows = c(56, 7, 14, 21, 10))
  expect_identical(twa$window, c(7, 10, 14, 21, 56))
  expect_identical(twa$start, rep(0, 5))
  # Means 10, 6, 2 | 12, 5, 1 on days 0, 7, 14 | 14, 21, 56
  area <- c(56, 56 + (6 + (6 - 4 * 3 / 7)) / 2 * 3, 84, 84 + 59.5, 143.5 + 105)
  expect_equal(twa$twa, area / c(7, 10, 14, 21, 56))
  expect_identical(twa$status, rep("ok", 5))
  expect_identical(twa$n_samples, c(2L, 2L, 3L, 5L, 6L))
  expect_identical(twa$method, rep("trapezoid", 5))

  m <- calc_twa(e, windows = c(7, 56), by = "mesocosm")
  expect_identical(m$mesocosm, c("T1", "T1", "T2", "T2"))
  expect_equal(m$twa[c(1, 3)], c(7, 9))
  expect_equal(mean(m$twa[c(2, 4)]), area[5] / 56)

  later <- calc_twa(e, windows = 7, start = 14)
  expect_equal(later$twa, (12 + 5) / 2)
})

test_that("missing values are not replaced and gaps make a window not estimable", {
  dat <- make_twa()
  dat$value[dat$mesocosm == "T2" & dat$time == 21] <- NA
  twa <- calc_twa(meso_twa(dat), windows = c(14, 21, 60))
  expect_equal(twa$twa[1:2], c(84 / 14, (84 + (12 + 4) / 2 * 7) / 21))
  expect_identical(twa$n_missing, c(0L, 1L, 1L))
  expect_identical(twa$status[3], "not_estimable")
  expect_match(twa$reason[3], "after the last measured value \\(time 56\\)")

  no_before <- make_twa()
  no_before <- no_before[!(no_before$timing %in% "before"), ]
  twa <- calc_twa(meso_twa(no_before), windows = c(7, 14))
  expect_identical(twa$status, c("ok", "not_estimable"))
  expect_match(twa$reason[2], "from time 7 to 14")
  twa <- calc_twa(meso_twa(), windows = 7, start = -2)
  expect_match(twa$reason, "starts before the first measured value")
})

test_that("model time-weighted averages integrate the fitted curves", {
  e <- meso_exposure_small()
  dis <- fit_dissipation(e)
  twa <- calc_twa(e, windows = c(7, 21, 35, 56), method = "model", fit = dis)
  expect_equal(twa$twa[1], 100 / 0.1 * (1 - exp(-0.7)) / 7 / 10, tolerance = 1e-8)
  curve <- function(t) predict(dis, times = t)$fit[predict(dis, times = t)$treatment == "100"]
  num <- function(b) {
    (stats::integrate(Vectorize(curve), 0, min(b, 14), rel.tol = 1e-10)$value +
       if (b > 14) stats::integrate(Vectorize(curve), 14, b, rel.tol = 1e-10)$value else 0) / b
  }
  h <- twa[twa$treatment == "100", ]
  expect_equal(h$twa, vapply(c(7, 21, 35, 56), num, numeric(1)), tolerance = 1e-8)
  expect_identical(h$status, c("ok", "ok", "ok", "extrapolated"))
  expect_match(h$reason[4], "time 35")

  expect_match(calc_twa(e, windows = 7, start = -1, method = "model", fit = dis)$reason[1],
               "start at the first application")
  first <- fit_dissipation(e, periods = 1)
  expect_match(calc_twa(e, windows = 21, method = "model", fit = first)$reason[1],
               "period 2 was not fitted")
  m <- calc_twa(e, windows = 7, method = "model", fit = fit_dissipation(e, by = "mesocosm"))
  expect_identical(m$mesocosm, c("L1", "L2", "H1", "H2"))
})

test_that("time-weighted average options are checked", {
  e <- meso_exposure_small()
  dis <- fit_dissipation(e)
  expect_error(calc_twa(e), "`windows` is missing", class = "cosmet_error_input")
  expect_error(calc_twa(e, windows = c(7, -1)), "positive", class = "cosmet_error_input")
  expect_error(calc_twa(e, windows = 7, method = "model"), "`fit` is missing",
               class = "cosmet_error_input")
  expect_error(calc_twa(e, windows = 7, fit = dis), "only used with method",
               class = "cosmet_error_input")
  expect_error(calc_twa(e, windows = 7, method = "model", fit = dis, by = "mesocosm"),
               "one fit per treatment", class = "cosmet_error_input")
  expect_error(calc_twa(meso_twa(), windows = 7, method = "model", fit = dis),
               "does not come from `x`", class = "cosmet_error_input")
})

# plot_exposure() and data --------------------------------------------------

test_that("exposure plots are built", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("viridisLite")
  e <- azoxy_exposure_lake()
  dis <- fit_dissipation(e)
  p <- plot_exposure(e, fit = dis)
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_no_warning(ggplot2::ggplot_build(plot_exposure(e, show = "both", scale = "log")))
  expect_no_warning(ggplot2::ggplot_build(
    plot_exposure(e, fit = fit_dissipation(e, by = "mesocosm"), show = "replicates",
                  show_control = TRUE)
  ))
  expect_error(plot_exposure(e, show = "points"), "must be one of",
               class = "cosmet_error_input")
  expect_error(plot_exposure(meso_twa(), fit = dis), "does not come from `x`",
               class = "cosmet_error_input")
})

test_that("the exposure plot explains how to install a missing package", {
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_error(plot_exposure(meso_exposure_small()), class = "cosmet_error_missing_package")
})

test_that("the azoxystrobin exposure data are complete", {
  expect_identical(dim(azoxy_exposure), c(126L, 7L))
  expect_identical(sort(unique(azoxy_exposure$day)), c(0, 7, 14, 21, 56))
  expect_identical(as.vector(table(azoxy_exposure$habitat)), c(63L, 63L))
  expect_true(all(is.na(azoxy_exposure$azoxystrobin) ==
                    azoxy_exposure$flag %in% c("not_detected", "below_0.1", "lost")))
  expect_true(all(azoxy_exposure$azoxystrobin[azoxy_exposure$concentration > 0] > 0,
                  na.rm = TRUE))
})

# Designs: schedules, single applications, nominal concentrations -----------

test_that("each treatment follows its own schedule of applications", {
  e <- meso_schedule()
  expect_false(e$settings$shared_schedule)
  # Treatments in the order of their nominal concentrations
  expect_identical(e$settings$treatments, c("control", "weekly", "single"))
  expect_identical(e$applications$treatment, c(rep("weekly", 4), "single"))
  expect_identical(e$applications$application, c(1:4, 1L))
  expect_identical(e$applications$nominal, c(5, 5, 5, 5, 20))
  s1 <- e$data[e$data$mesocosm == "S1", ]
  expect_identical(unique(s1$period), 1L)
  expect_identical(s1$elapsed, s1$time)
  w1 <- e$data[e$data$mesocosm == "W1", ]
  expect_identical(w1$period, c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L))
  expect_identical(w1$application, rep(c(0, 7, 14, 21), each = 3))
  # The control has no applications; its labels are kept
  c1 <- e$data[e$data$mesocosm == "C1", ]
  expect_identical(unique(c1$period), 0L)
  expect_identical(sum(c1$timing %in% "before"), 3L)
  expect_snapshot(print(e))

  shared <- meso_exposure_small()
  expect_true(shared$settings$shared_schedule)
  expect_identical(shared$applications$treatment, rep(c("0", "10", "100"), each = 2))
  expect_identical(shared$applications$nominal, c(0, 0, 10, 10, 100, 100))
})

test_that("fits, predictions and averages use the schedule of each treatment", {
  e <- meso_schedule()
  dis <- fit_dissipation(e)
  est <- dis$estimates
  expect_identical(est$treatment, c(rep("weekly", 4), "single"))
  expect_identical(est$period, c(1:4, 1L))
  expect_identical(est$application, c(0, 7, 14, 21, 0))
  expect_equal(est$k, rep(0.1, 5), tolerance = 1e-8)
  build <- 5 * cumsum(exp(-0.7 * (0:3)))
  expect_equal(est$c0, c(build, 20), tolerance = 1e-8)
  expect_equal(est$c0_pct_nominal, c(100 * build / 5, 100), tolerance = 1e-8)
  expect_identical(fit_dissipation(e, periods = 2)$estimates$treatment, "weekly")
  expect_error(fit_dissipation(e, periods = 5), "from 1 to 4", class = "cosmet_error_input")
  expect_equal(fit_dissipation(e, method = "nonlinear")$estimates$c0, c(build, 20),
               tolerance = 1e-6)

  p <- predict(dis, times = c(7, 21))
  expect_equal(p$fit[p$treatment == "single"], 20 * exp(-0.1 * c(7, 21)), tolerance = 1e-8)
  expect_identical(p$period[p$treatment == "weekly"], c(2L, 4L))

  # Model averages: the integral of the sum of the doses
  twa <- calc_twa(e, windows = c(14, 28), method = "model", fit = dis)
  dose_integral <- function(times, doses, b) sum(doses * (1 - exp(-0.1 * (b - times))) / 0.1)
  expect_equal(twa$twa[twa$treatment == "single"],
               c(dose_integral(0, 20, 14) / 14, dose_integral(0, 20, 28) / 28), tolerance = 1e-8)
  expect_equal(twa$twa[twa$treatment == "weekly"],
               c(dose_integral(c(0, 7), c(5, 5), 14) / 14,
                 dose_integral(c(0, 7, 14, 21), rep(5, 4), 28) / 28), tolerance = 1e-8)
  expect_equal(twa$twa_pct_nominal, 100 * twa$twa / twa$concentration)
  trap <- calc_twa(e, windows = 28)
  expect_equal(trap$twa_pct_nominal, 100 * trap$twa / c(5, 20))
  expect_true(all(trap$twa > twa$twa[twa$window == 28]))
})

test_that("windows start at the first application of each treatment", {
  e <- meso_schedule(make_schedule(first = 2))
  twa <- calc_twa(e, windows = 7)
  expect_identical(twa$treatment, c("weekly", "single"))
  expect_identical(twa$start, c(2, 0))
  expect_identical(twa$end, c(9, 7))
  expect_identical(calc_twa(e, windows = 7, start = 1)$start, c(1, 1))
})

test_that("a single application needs no timing column", {
  dat <- make_exposure()
  dat <- dat[dat$time %in% c(0, 7) & !(dat$timing %in% "before"), ]
  e <- meso_exposure(dat[c("mesocosm", "concentration", "time", "value")], control = "0",
                     applications = 0, treatment_col = "concentration",
                     concentration_col = "concentration")
  expect_true(e$settings$timing_assumed)
  est <- fit_dissipation(e)$estimates
  expect_equal(est$k, c(0.1, 0.1), tolerance = 1e-8)
  expect_identical(est$df, c(2L, 2L))
  cn <- compare_nominal(e)
  expect_equal(cn$pct_nominal, c(100, 100))
  expect_identical(cn$status, c("no_before", "no_before"))
  expect_equal(calc_twa(e, windows = 7)$twa_pct_nominal, rep(100 * (1 + exp(-0.7)) / 2, 2))
})

test_that("measured concentrations are compared with the nominal ones", {
  e <- meso_schedule()
  cn <- compare_nominal(e)
  expect_identical(cn$treatment, c(rep("weekly", 4), "single"))
  expect_identical(cn$application, c(1:4, 1L))
  build <- 5 * cumsum(exp(-0.7 * (0:3)))
  expect_equal(cn$after, c(build, 20), tolerance = 1e-12)
  expect_equal(cn$before, c(NA, build[2:4] - 5, NA), tolerance = 1e-12)
  expect_equal(cn$pct_nominal, c(100 * build / 5, 100), tolerance = 1e-12)
  expect_equal(cn$pct_added, c(NA, 100, 100, 100, NA), tolerance = 1e-12)
  expect_identical(cn$status, c("no_before", "ok", "ok", "ok", "no_before"))
  expect_identical(cn$n_after, rep(2L, 5))
  expect_identical(cn$time_before, c(NA, 7, 14, 21, NA))

  # Means of the mesocosms with a value; missing values are not replaced
  x <- make_schedule()
  x$data$value[x$data$mesocosm == "W1" & x$data$time == 7 & x$data$timing == "after"] <- NA
  x$data$value[x$data$mesocosm == "W2" & x$data$time == 7 & x$data$timing == "after"] <- 9
  cn <- compare_nominal(meso_schedule(x))
  expect_identical(cn$n_after[2], 1L)
  expect_equal(cn$pct_nominal[2], 180)
  m <- compare_nominal(meso_schedule(x), by = "mesocosm")
  w1 <- m[m$mesocosm == "W1" & m$application == 2L, ]
  expect_identical(w1$status, "not_estimable")
  expect_match(w1$reason, "no value right after")

  # Samples taken a while after the application count with `lag`
  late <- make_schedule()
  first <- late$data$time == 0
  late$data$time[first] <- 1 / 24
  late$data$timing[first] <- NA
  late$data$value[first] <- late$data$value[first] * exp(-0.1 / 24)
  e <- meso_exposure(late$data, control = "control", applications = late$schedule,
                     timing_col = "timing")
  expect_identical(compare_nominal(e)$status[1], "not_estimable")
  cn <- compare_nominal(e, lag = 1 / 24)
  first <- cn$application == 1L
  expect_identical(cn$time_after[first], c(1 / 24, 1 / 24))
  expect_equal(cn$pct_nominal[first], 100 * exp(-0.1 / 24) * c(1, 1))

  unknown <- meso_exposure(make_exposure(), control = "0", applications = c(0, 14),
                           treatment_col = "concentration", timing_col = "timing")
  cn <- compare_nominal(unknown)
  expect_identical(unique(cn$status), "not_estimable")
  expect_match(cn$reason[1], "unknown")
  expect_true(all(is.na(calc_twa(unknown, windows = 7)$twa_pct_nominal)))
})

test_that("invalid schedules give informative errors", {
  x <- make_schedule()
  ex <- function(schedule) {
    meso_exposure(x$data, control = "control", applications = schedule, timing_col = "timing")
  }
  expect_error(ex(cbind(x$schedule, dose = 1)), "not used: 'dose'", class = "cosmet_error_input")
  expect_error(ex(x$schedule[-1, ]), "without applications.*'single'",
               class = "cosmet_error_input")
  expect_error(ex(rbind(x$schedule, x$schedule[1, ])), "repeats an application",
               class = "cosmet_error_input")
  bad <- x$schedule
  bad$treatment[1] <- "double"
  expect_error(ex(bad), "not in column 'treatment'", class = "cosmet_error_input")
  bad <- rbind(x$schedule, data.frame(treatment = "control", time = 0, nominal = 1))
  expect_error(ex(bad), "cannot receive", class = "cosmet_error_input")
  bad <- x$schedule
  bad$nominal[2] <- -5
  expect_error(ex(bad), "nominal", class = "cosmet_error_input")
  expect_error(ex(x$schedule["treatment"]), "'time'", class = "cosmet_error_input")
  ok <- rbind(x$schedule, data.frame(treatment = "control", time = 0, nominal = 0))
  expect_identical(nrow(ex(ok)$applications), 6L)
})

test_that("treatments with the same nominal concentration keep their own labels", {
  lab <- exposure_labels(c("0", "pulse", "split"), c(0, 10, 10), "0", "ug/L")
  expect_identical(unname(lab), c("Control", "pulse (10 ug/L)", "split (10 ug/L)"))
  lab <- exposure_labels(c("0", "2", "20"), c(0, 2, 20), "0", "ug/L")
  expect_identical(unname(lab), c("Control", "2 ug/L", "20 ug/L"))
  lab <- exposure_labels(c("C", "low", "high"), c(NA, NA, NA), "C", NULL)
  expect_identical(unname(lab), c("Control", "low", "high"))
  expect_error(exposure_labels(c("0", "a", "b"), c(0, 1, 2), "0", NULL, labels = c("C", "x", "x")),
               "different for each treatment", class = "cosmet_error_input")
})

test_that("the plot follows the schedule and shows the nominal concentration", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("viridisLite")
  e <- meso_schedule()
  p <- plot_exposure(e, fit = fit_dissipation(e))
  b <- ggplot2::ggplot_build(p)
  vlines <- b$data[[1]]
  expect_identical(as.vector(table(vlines$PANEL)), c(4L, 1L))
  hlines <- b$data[[2]]
  expect_equal(sort(hlines$yintercept), c(5, 20))
  expect_no_warning(ggplot2::ggplot_build(plot_exposure(e, show_nominal = FALSE,
                                                        scale = "log")))
})
