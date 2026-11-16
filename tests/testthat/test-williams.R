test_that("the isotonic mean of the highest dose is the largest tail mean", {
  # Tail means of (1, 3, 2) with equal weights: 2, 2.5 and 2 -> 2.5
  expect_equal(williams_last_iso(c(1, 3, 2), c(2, 2, 2)), 2.5)
  # Weights matter: (3 * 1 + 2 * 3) / 4 = 2.25
  expect_equal(williams_last_iso(c(3, 2), c(1, 3)), 2.25)
})

test_that("the simulated null distribution matches the tables of Williams (1972)", {
  # One-sided 5% critical values for equal replication, from the table of
  # Williams (1972) as transcribed in PMCMRplus (TabCrit$williams.tk005):
  # columns = number of doses, rows = degrees of freedom.
  cases <- list(
    list(n = rep(3, 3), table = 2.058),  # 2 doses, df 6
    list(n = rep(3, 4), table = 1.997),  # 3 doses, df 8
    list(n = rep(3, 5), table = 1.956),  # 4 doses, df 10
    list(n = rep(4, 4), table = 1.903)   # 3 doses, df 12
  )
  for (cs in cases) {
    null <- withr::with_seed(1, williams_null(cs$n, 5e5))
    # Monte Carlo standard error of the quantile: about 0.003
    expect_lt(abs(stats::quantile(null, 0.95, names = FALSE) - cs$table), 0.012)
  }
})

test_that("with a single dose Williams' test is Student's t test", {
  x <- make_monotone()
  x <- x[x$concentration %in% c(0, 10), ]
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  tw <- test_control(d, method = "williams", direction = "decrease")
  tt <- stats::t.test(x$value[x$concentration == 10], x$value[x$concentration == 0],
                      var.equal = TRUE, alternative = "less")
  expect_equal(tw$tests$p_value, tt$p.value)
  expect_equal(tw$tests$statistic, -unname(tt$statistic))
  expect_equal(tw$tests$critical_value, stats::qt(0.95, 4))
  both <- test_control(d, method = "williams", direction = "both")
  dec <- both$tests[both$tests$direction == "decrease", ]
  expect_equal(dec$p_adjusted, 2 * tt$p.value)
})

test_that("the statistic matches a hand calculation and the procedure steps down", {
  d <- meso_monotone()
  tw <- test_control(d, method = "williams", direction = "decrease", step_variance = "reduced")
  top <- tw$tests[tw$tests$concentration == 10, ]
  # Decrease: means 11 (control), 10.667, 7, 2; isotonic mean of the highest
  # dose = 2; pooled SS = 2 + 8.667 + 2 + 2 over 8 df.
  ss <- 2 + sum((c(9, 13, 10) - 32 / 3)^2) + 2 + 2
  expect_equal(top$statistic, (11 - 2) / sqrt(ss / 8 * (1 / 3 + 1 / 3)))
  expect_identical(tw$tests$tested, c(TRUE, TRUE, TRUE))
  expect_identical(tw$tests$significant, c(FALSE, TRUE, TRUE))
  # The step of dose 3 uses only the control and doses 1 and 3
  s3 <- tw$tests[tw$tests$concentration == 3, ]
  ss3 <- 2 + sum((c(9, 13, 10) - 32 / 3)^2) + 2
  expect_equal(s3$se, sqrt(ss3 / 6) * sqrt(2 / 3))
  expect_identical(s3$df, 6)
  # With the variance of all groups every step has the 8 df of the full design
  ta <- test_control(d, method = "williams", direction = "decrease", step_variance = "all")
  expect_identical(ta$tests$df, c(8, 8, 8))
  expect_equal(ta$tests$se, rep(sqrt(ss / 8) * sqrt(2 / 3), 3))
  expect_equal(ta$tests$statistic[3], (11 - 2) / sqrt(ss / 8 * (2 / 3)))
  expect_identical(ta$tests$significant, c(FALSE, TRUE, TRUE))
})

test_that("isotonic means of all doses match a pool-adjacent-violators fit", {
  # Equal weights: the max-min estimates equal stats::isoreg()
  mu <- c(3, 1, 4, 1.5, 5, 9, 2.5)
  expect_equal(williams_iso_all(mu, rep(2, 7)), stats::isoreg(mu)$yf)
  # Weighted: pooling (3, 1) with weights (1, 3) gives 1.5 for both
  expect_equal(williams_iso_all(c(3, 1, 4), c(1, 3, 2)), c(1.5, 1.5, 4))
})

test_that("critical values of each step follow the tables of Williams (1972)", {
  d <- meso_monotone()
  # 4 groups of 3: 8 df with the variance of all groups; 6 and 4 df reduced.
  ta <- test_control(d, method = "williams", direction = "decrease", step_variance = "all",
                     n_sim = 199999)
  expect_lt(abs(ta$tests$critical_value[3] - 1.997), 0.015)  # 3 doses, 8 df
  expect_lt(abs(ta$tests$critical_value[2] - 1.962), 0.015)  # 2 doses, 8 df
  expect_equal(ta$tests$critical_value[1], stats::qt(0.95, 8))
  tr <- test_control(d, method = "williams", direction = "decrease", step_variance = "reduced",
                     n_sim = 199999)
  expect_lt(abs(tr$tests$critical_value[2] - 2.058), 0.015)  # 2 doses, 6 df
  expect_equal(tr$tests$critical_value[1], stats::qt(0.95, 4))
})

test_that("the step-down procedure stops at the first non-significant dose", {
  x <- make_monotone()
  x$value[x$concentration == 3] <- c(10, 12, 11)  # dose 3 like the control
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  tw <- test_control(d, method = "williams", direction = "decrease", step_variance = "reduced")
  expect_identical(tw$tests$tested, c(FALSE, TRUE, TRUE))
  expect_identical(tw$tests$significant, c(NA, FALSE, TRUE))
  expect_true(is.na(tw$tests$p_value[1]))
  # MDD quantities exist for every dose, tested or not
  expect_false(anyNA(tw$tests$se))
})

test_that("a seed makes results reproducible and leaves the global stream alone", {
  d <- meso_monotone()
  a <- test_control(d, method = "williams", n_sim = 999, seed = 5)
  b <- test_control(d, method = "williams", n_sim = 999, seed = 5)
  expect_identical(a$tests, b$tests)
  set.seed(10)
  r1 <- stats::runif(1)
  set.seed(10)
  invisible(test_control(d, method = "williams", n_sim = 999, seed = 5))
  r2 <- stats::runif(1)
  expect_identical(r1, r2)
  # Without a seed, the current stream is used
  set.seed(3)
  c1 <- test_control(d, method = "williams", n_sim = 999, seed = NULL)
  set.seed(3)
  c2 <- test_control(d, method = "williams", n_sim = 999, seed = NULL)
  expect_identical(c1$tests, c2$tests)
})

test_that("the result of a taxon does not depend on the other taxa analysed", {
  d <- azoxy_macro()
  all <- test_control(d, method = "williams", n_sim = 9999)
  one <- test_control(d, method = "williams", n_sim = 9999, taxa = "Micronecta")
  expect_identical(one$tests, all$tests[all$tests$taxon == "Micronecta", ], ignore_attr = TRUE)
})
