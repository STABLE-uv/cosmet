test_that("meso_data() returns a complete, ordered long table", {
  d <- meso_monotone(times = c(14, 7))
  expect_s3_class(d, "meso_data")
  expect_identical(names(d$data), c("mesocosm", "time", "treatment", "concentration", "taxon",
                                    "value", "response", "value_status"))
  expect_identical(nrow(d$data), 24L)
  expect_identical(d$data$mesocosm, rep(sort(unique(d$data$mesocosm), method = "radix"),
                                        each = 2))
  expect_identical(d$data$time, rep(c(7, 14), 12))
  expect_identical(d$settings$treatments, c("0", "1", "3", "10"))
  expect_identical(d$data$response, d$data$value)
  expect_true(all(d$data$value_status == "recorded"))
  expect_identical(nrow(d$design), 12L)
})

test_that("row order of the input does not change the result", {
  x <- make_monotone(c(7, 14))
  ref <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                   concentration_col = "concentration")
  set.seed(4)
  shuffled <- x[sample.int(nrow(x)), ]
  out <- meso_data(shuffled, control = "0", zero_fill = FALSE,
                   treatment_col = "concentration", concentration_col = "concentration")
  expect_identical(out$data, ref$data)
  expect_identical(out$design, ref$design)
})

test_that("zero_fill decides whether a taxon without a row is zero or unknown", {
  x <- make_presence()
  z <- meso_data(x, control = "ctrl", zero_fill = TRUE)
  u <- meso_data(x, control = "ctrl", zero_fill = FALSE)
  pick <- function(d, m, t, tx) d$data[d$data$mesocosm == m & d$data$time == t & d$data$taxon == tx, ]
  expect_identical(pick(z, "C2", 7, "Caenis")$value, 0)
  expect_identical(pick(z, "C2", 7, "Caenis")$value_status, "zero_filled")
  expect_true(is.na(pick(u, "C2", 7, "Caenis")$value))
  expect_identical(pick(u, "C2", 7, "Caenis")$value_status, "not_recorded")
  expect_identical(z$diagnostics$n_zero_filled, 7L)
  expect_identical(u$diagnostics$n_not_recorded, 7L)
})

test_that("samples not taken are not filled, and `samples` adds empty samples", {
  x <- make_presence()
  x <- x[!(x$mesocosm == "T2" & x$time == 14), ]
  d <- meso_data(x, control = "ctrl", zero_fill = TRUE)
  expect_false(any(d$data$mesocosm == "T2" & d$data$time == 14))
  chk <- check_design(d)
  expect_identical(chk$missing_samples$mesocosm, "T2")
  expect_false(chk$balanced)
  # The sample was analysed and nothing was found: listing it gives zeros.
  smp <- data.frame(mesocosm = rep(c("C1", "C2", "T1", "T2"), 2), time = rep(c(7, 14), each = 4))
  d2 <- meso_data(x, control = "ctrl", zero_fill = TRUE, samples = smp)
  t2 <- d2$data[d2$data$mesocosm == "T2" & d2$data$time == 14, ]
  expect_identical(t2$value, c(0, 0))
  expect_true(check_design(d2)$balanced)
})

test_that("`samples` must cover every sample with data and only known mesocosms", {
  x <- make_presence()
  smp <- data.frame(mesocosm = c("C1", "C2", "T1"), time = 7)
  expect_error(meso_data(x, control = "ctrl", zero_fill = TRUE, samples = smp),
               class = "cosmet_error_input")
  smp2 <- data.frame(mesocosm = c("C1", "C2", "T1", "T2", "X9"), time = 7)
  expect_error(meso_data(x[x$time == 7, ], control = "ctrl", zero_fill = TRUE, samples = smp2),
               "treatment is unknown")
})

test_that("a recorded NA stays unknown", {
  x <- make_monotone()
  x$value[2] <- NA
  d <- meso_data(x, control = "0", zero_fill = TRUE, treatment_col = "concentration",
                 concentration_col = "concentration")
  expect_identical(sum(d$data$value_status == "missing"), 1L)
  expect_identical(sum(is.na(d$data$value)), 1L)
})

test_that("input errors are classed and point to the fix", {
  x <- make_monotone()
  expect_error(meso_data(x, zero_fill = TRUE, treatment_col = "concentration"),
               class = "cosmet_error_input")
  expect_error(meso_data(x, control = "0", treatment_col = "concentration"),
               "zero_fill")
  expect_error(meso_data(x, control = "5", zero_fill = TRUE, treatment_col = "concentration"),
               "not a level")
  dup <- rbind(x, x[1, ])
  expect_error(meso_data(dup, control = "0", zero_fill = TRUE, treatment_col = "concentration"),
               class = "cosmet_error_input")
  moving <- make_monotone(c(7, 14))
  moving$concentration[moving$mesocosm == "A1" & moving$time == 14] <- 3
  expect_error(meso_data(moving, control = "0", zero_fill = TRUE,
                         treatment_col = "concentration", concentration_col = "concentration"),
               "single treatment")
  bad_ctrl <- x
  bad_ctrl$treatment <- paste0("T", bad_ctrl$concentration)
  bad_ctrl$treatment[bad_ctrl$concentration == 1] <- "C"
  expect_error(meso_data(bad_ctrl, control = "C", zero_fill = TRUE,
                         concentration_col = "concentration"),
               "control concentration must be 0")
  dates <- x
  dates$time <- as.Date("2024-05-01")
  expect_error(meso_data(dates, control = "0", zero_fill = TRUE, treatment_col = "concentration"),
               "dates")
  fac <- x
  fac$value <- factor(fac$value)
  expect_error(meso_data(fac, control = "0", zero_fill = TRUE, treatment_col = "concentration"),
               "factor")
  expect_error(meso_data(x[0, ], control = "0", zero_fill = TRUE, treatment_col = "concentration"),
               "no rows")
  expect_error(meso_data(x, control = "0", zero_fill = TRUE, taxon_col = "species",
                         treatment_col = "concentration"), "taxon_col")
})

test_that("a numeric control is matched as text", {
  d <- meso_data(make_monotone(), control = 0, zero_fill = FALSE,
                 treatment_col = "concentration", concentration_col = "concentration")
  expect_identical(d$settings$control, "0")
})

test_that("tibbles are accepted and grouped data frames are not", {
  skip_if_not_installed("tibble")
  x <- make_monotone()
  ref <- meso_data(x, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                   concentration_col = "concentration")
  tb <- meso_data(tibble::as_tibble(x), control = "0", zero_fill = FALSE,
                  treatment_col = "concentration", concentration_col = "concentration")
  expect_identical(tb$data, ref$data)
  grouped <- structure(x, class = c("grouped_df", "tbl_df", "tbl", "data.frame"))
  expect_error(meso_data(grouped, control = "0", zero_fill = FALSE,
                         treatment_col = "concentration"), "grouped")
})

test_that("check_design() counts replicates and reports problems", {
  d <- meso_monotone(times = c(7, 14))
  chk <- check_design(d)
  expect_s3_class(chk, "meso_design_check")
  expect_identical(chk$replicates$treatment, c("0", "1", "3", "10"))
  expect_identical(chk$replicates[["7"]], rep(3L, 4))
  expect_identical(chk$issues, character())
  expect_snapshot(print(chk))
})

test_that("transform_response() computes ln(A x + 1) and records A", {
  d <- transform_response(meso_monotone(), A = 2)
  expect_equal(d$data$response, log(2 * d$data$value + 1))
  expect_identical(d$settings$transform, list(method = "log", A = 2))
  # Re-applying starts again from the values
  d2 <- transform_response(d, A = 10)
  expect_equal(d2$data$response, log(10 * d$data$value + 1))
  expect_error(transform_response(meso_monotone()), "no default")
  expect_error(transform_response(meso_monotone(), A = 0), class = "cosmet_error_input")
  neg <- make_monotone()
  neg$value[1] <- -1
  dn <- meso_data(neg, control = "0", zero_fill = FALSE, treatment_col = "concentration")
  expect_error(transform_response(dn, A = 2), "non-negative")
  expect_snapshot(print(d))
})

test_that("taxon_eligibility() uses the controls only", {
  x <- make_presence()
  d <- meso_data(x, control = "ctrl", zero_fill = TRUE)
  e <- taxon_eligibility(d, min_control_mean = 2)
  # Caenis in controls: 3 (C1, t7), 0 (C2, t7), 0 (C1, t14), 0 (C2, t14) -> mean 0.75
  expect_identical(e$taxon, c("Caenis", "Cloeon"))
  expect_equal(e$control_mean, c(0.75, 10.5))
  expect_equal(e$control_max_time_mean, c(1.5, 11))
  expect_equal(e$control_frequency, c(0.25, 1))
  expect_identical(e$eligible, c(FALSE, TRUE))
  expect_true(all(taxon_eligibility(d)$eligible))
  expect_identical(taxon_eligibility(d, min_control_frequency = 0.5)$eligible, c(FALSE, TRUE))
})
