test_that("NOEC and LOEC follow the step-down Williams test", {
  tw <- test_control(meso_monotone(), method = "williams")
  n <- derive_noec(tw)$noec
  expect_identical(n$direction, "decrease")
  expect_identical(n$loec, 3)
  expect_identical(n$noec, 1)
  expect_identical(n$noec_label, "1")
  expect_true(n$crit_a)
  expect_true(n$consistent)
  expect_true(n$selected)
})

test_that("no effect gives '>max' and an effect at the lowest gives '<min'", {
  x <- make_monotone()
  x$value <- c(10, 12, 11, 10, 12, 11, 11, 10, 12, 12, 11, 10)
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  n <- derive_noec(test_control(d, method = "williams"))$noec
  expect_identical(n$noec_label, ">10")
  expect_identical(n$noec, 10)
  expect_false(n$crit_a)
  expect_false(n$selected)
  x$value <- c(10, 12, 11, 3, 4, 2, 2, 3, 2, 1, 2, 1)
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  n <- derive_noec(test_control(d, method = "williams"))$noec
  expect_identical(n$noec_label, "<1")
  expect_true(is.na(n$noec))
  expect_identical(n$loec, 1)
})

test_that("effects in both directions are ambiguous, not forced into one NOEC", {
  x <- make_monotone()
  x$value <- c(10, 10.5, 9.5, 30, 30.5, 29.5, 10, 10.5, 9.5, 0, 0.5, 0.2)
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  n <- derive_noec(test_control(d, method = "williams"))$noec
  expect_identical(n$status, "ambiguous_direction")
  expect_true(is.na(n$noec_label))
  expect_false(n$selected)
})

test_that("criterion b compares the right mean with the threshold", {
  tw <- test_control(meso_monotone(), method = "williams")
  # Decrease: the control mean (11) is compared with the threshold.
  n <- derive_noec(tw, criteria = noec_criteria(min_abundance = 12))$noec
  expect_identical(n$mean_for_threshold, 11)
  expect_false(n$crit_b)
  expect_false(n$selected)
  n <- derive_noec(tw, criteria = noec_criteria(min_abundance = 11))$noec
  expect_true(n$crit_b)
  expect_true(n$selected)
})

test_that("criterion c requires means that follow the direction of the effect", {
  x <- make_monotone()
  x$value[x$concentration == 1] <- c(12, 13, 14)  # 1 above the control: not monotonic
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  tw <- test_control(d, method = "williams")
  n <- derive_noec(tw, criteria = noec_criteria(monotonic = TRUE))$noec
  expect_identical(n$direction, "decrease")
  expect_false(n$crit_c)
  expect_false(n$selected)
  # A tolerance as large as the reversal (2 units) lets it pass
  n2 <- derive_noec(tw, criteria = noec_criteria(monotonic = TRUE, monotonic_tol = 2))$noec
  expect_true(n2$crit_c)
})

test_that("criterion d rejects effects already present before the application", {
  x <- make_monotone(times = c(-7, 7))
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  tw <- test_control(d, method = "williams")
  n <- derive_noec(tw, pre_times = -7, criteria = noec_criteria(pretreatment = TRUE))$noec
  expect_identical(n$is_post, c(FALSE, TRUE))
  expect_identical(n$crit_d, c(NA, FALSE))
  expect_identical(n$selected, c(FALSE, FALSE))
  # Without the criterion, the post-application effect is selected
  n0 <- derive_noec(tw, pre_times = -7)$noec
  expect_identical(n0$selected, c(FALSE, TRUE))
  expect_error(derive_noec(tw, criteria = noec_criteria(pretreatment = TRUE)), "pre_times")
  expect_error(derive_noec(tw, pre_times = -1), "without data")
})

test_that("criterion d is unknown when a pre-treatment test is not available", {
  x <- make_monotone(times = c(-7, 7))
  x$value[x$time == -7] <- 0  # absent before the application
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  n <- derive_noec(test_control(d, method = "williams"), pre_times = -7,
                   criteria = noec_criteria(pretreatment = TRUE))$noec
  expect_identical(n$status, c("absent", "ok"))
  expect_identical(n$crit_d, c(NA, NA))
  expect_false(n$selected[2])
})

test_that("abundance thresholds are converted between compatible units", {
  expect_equal(convert_abundance(10, "cells/mL", "cells/L"), 10000)
  expect_equal(convert_abundance(3, "individuals/L", "individuals/m3"), 3000)
  expect_error(convert_abundance(3, "individuals/L", "cells/L"), "Cannot convert")
  expect_error(convert_abundance(3, "individuals/sample", "individuals/L"), "Cannot convert")
  tw <- test_control(meso_monotone(unit = "cells/L"), method = "williams")
  n <- derive_noec(tw, criteria = noec_criteria(min_abundance = 0.01,
                                                abundance_unit = "cells/mL"))$noec
  expect_equal(n$threshold, 10)
  expect_identical(n$threshold_unit, "cells/L")
  tw0 <- test_control(meso_monotone(), method = "williams")
  expect_error(derive_noec(tw0, criteria = noec_criteria(min_abundance = 3,
                                                        abundance_unit = "individuals/sample")),
               "unit of the data is unknown")
})

test_that("s11_criteria() uses the thresholds of Text S11", {
  s <- s11_criteria("macroinvertebrates")
  expect_s3_class(s, "cosmet_criteria")
  expect_identical(s$min_abundance, 3)
  expect_identical(s$abundance_unit, "individuals/sample")
  expect_true(s$monotonic)
  expect_true(s$pretreatment)
  expect_identical(s11_criteria("phytoplankton_small")$min_abundance, 100)
  expect_error(s11_criteria("fish"), "must be one of")
  expect_error(s11_criteria(), "missing")
  expect_snapshot(print(s))
})

test_that("Dunnett NOECs flag concentrations without effect above the LOEC", {
  x <- make_monotone()
  x$value <- c(10, 10.5, 9.5, 2, 2.5, 1.5, 10, 10.5, 9.5, 9, 10.5, 10)
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  n <- derive_noec(test_control(d, method = "dunnett"))$noec
  expect_identical(n$direction, "decrease")
  expect_identical(n$noec_label, "<1")
  expect_false(n$consistent)
})

test_that("NOECs need concentrations", {
  x <- make_presence()
  d <- meso_data(x, control = "ctrl", zero_fill = TRUE)
  td <- test_control(d, method = "dunnett")
  expect_error(derive_noec(td), "concentrations")
  expect_error(test_control(d, method = "williams"), "concentration")
  expect_error(test_control(d), "no default")
})

test_that("noec_table() lays out one row per taxon and one column per time", {
  x <- make_monotone(times = c(-7, 7))
  x$value[x$time == -7] <- c(10, 12, 11, 10, 12, 11, 11, 10, 12, 12, 11, 10)
  d <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                 concentration_col = "concentration")
  noec <- derive_noec(test_control(d, method = "williams"), pre_times = -7)
  tab <- noec_table(noec)
  expect_identical(names(tab), c("taxon", "-7", "7"))
  expect_identical(tab[["-7"]], ">10")
  expect_identical(tab[["7"]], "1\u2193")
  expect_identical(noec_table(noec, arrows = FALSE)[["7"]], "1")
  sel <- noec_table(noec, selected = TRUE)
  expect_identical(sel[["-7"]], NA_character_)
  expect_snapshot(print(noec))
})
