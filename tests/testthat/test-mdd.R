test_that("MDD follows eq. 2 and 3 of Brock et al. (2015) without transformation", {
  x <- make_monotone()
  x <- x[x$concentration %in% c(0, 10), ]
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  m <- calc_mdd(test_control(d, method = "dunnett"))
  # One treatment: Dunnett = t test. s^2 = (2 + 2) / 4 = 1; se = sqrt(2/3)
  mdd <- stats::qt(0.975, 4) * sqrt(2 / 3)
  expect_equal(m$mdd, mdd, tolerance = 1e-4)
  expect_equal(m$mdd_value, m$mdd)
  expect_equal(m$mdd_pct, 100 * mdd / 11, tolerance = 1e-4)
  expect_equal(m$control_value, 11)
  expect_equal(m$effect_pct, 100 * (2 - 11) / 11)
})

test_that("MDD is back-transformed with eq. 5 and 6 of Brock et al. (2015)", {
  d <- transform_response(meso_monotone(), A = 2)
  tw <- test_control(d, method = "williams")
  m <- calc_mdd(tw)
  y0 <- mean(log(2 * c(10, 12, 11) + 1))
  mdd <- m$critical_value * m$se
  expect_equal(m$mdd, mdd)
  expect_equal(m$mdd_value, (exp(y0) - exp(y0 - mdd)) / 2)
  expect_equal(m$control_value, rep((exp(y0) - 1) / 2, 3))
  expect_equal(m$mdd_pct, 100 * m$mdd_value / ((exp(y0) - 1) / 2))
  # One row per dose even though Williams tested two directions
  expect_identical(nrow(m), 3L)
  expect_identical(m$concentration, c(1, 3, 10))
})

test_that("MDD classes follow Table 1 of Brock et al. (2015)", {
  pct <- c(10, 49.9, 50, 69.9, 70, 89.9, 90, 100, 100.1, NA)
  expect_identical(mdd_class(pct),
                   c("IV", "IV", "III", "III", "II", "II", "I", "I", "0", NA))
})

test_that("a difference larger than the MDD is significant (Dunnett, two-sided)", {
  d <- transform_response(azoxy_macro(taxa = c("Chironomini", "Micronecta", "Caenidae")), A = 2)
  td <- test_control(d, method = "dunnett")
  m <- calc_mdd(td)
  ok <- m$status == "ok"
  sig <- td$tests$significant[ok]
  expect_identical(sig, abs(m$difference[ok]) > m$mdd[ok])
})

test_that("taxa not tested or with a zero control mean get NA", {
  d <- transform_response(azoxy_macro(taxa = "Cloeon"), A = 2)
  m <- calc_mdd(test_control(d, method = "dunnett"))
  ctrl_zero <- m$control_value %in% 0 | is.na(m$control_value)
  expect_true(all(is.na(m$mdd_pct[ctrl_zero])))
  expect_true(all(is.na(m$mdd_class[m$status != "ok"])))
})
