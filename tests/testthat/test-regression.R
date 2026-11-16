# Regression against the original analysis scripts of the authors
# (population_noec.R and PRC_macroinvertebrados.R): the fixtures hold their
# output on the same data.
# The script implements Williams' test with step_variance = "reduced".
# NOECs, LOECs, criteria and selections must be identical. Williams'
# p-values are identical for steps with one dose (Student's t) and for the
# first replication pattern simulated (all groups); for other patterns the
# random draws differ, so they agree within Monte Carlo error.

wide_to_meso <- function(dat, taxa_cols, unit = NULL) {
  long <- data.frame(
    mesocosm = rep(dat$tank_id, length(taxa_cols)), time = rep(dat$time, length(taxa_cols)),
    concentration = rep(dat$conc, length(taxa_cols)), taxon = rep(taxa_cols, each = nrow(dat)),
    value = unlist(dat[taxa_cols], use.names = FALSE), stringsAsFactors = FALSE
  )
  meso_data(long, control = "0", zero_fill = FALSE, treatment_col = "concentration",
            concentration_col = "concentration", unit = unit)
}

expect_same_noec <- function(new, ref) {
  new <- new[match(paste(ref$taxon, ref$time), paste(new$taxon, new$time)), ]
  expect_identical(new$status, ref$status)
  expect_identical(new$direction, ref$direction)
  expect_identical(new$noec_label, ref$NOEC)
  expect_identical(new$loec, ref$LOEC)
  expect_identical(new$crit_a, ref$crit_a)
  expect_identical(new$crit_b, ref$crit_b)
  expect_identical(new$crit_c, ref$crit_c)
  expect_identical(new$crit_d, ref$crit_d)
  expect_identical(new$selected, ref$selected)
  expect_equal(new$p_value, ref$p_williams)
}

test_that("NOECs reproduce the original script on its simulated data", {
  ref <- readRDS(test_path("fixtures", "original_noec_dummy.rds"))
  taxa_cols <- grep("^Taxon_", names(ref$dat), value = TRUE)
  d <- transform_response(wide_to_meso(ref$dat, taxa_cols), A = 33.3)
  tw <- test_control(d, method = "williams", step_variance = "reduced")
  crit <- noec_criteria(min_abundance = 10, monotonic = TRUE, pretreatment = TRUE)
  expect_same_noec(derive_noec(tw, pre_times = -7, criteria = crit)$noec, ref$summary)
  # Every step: identical statistics; p-values identical for 2 and 4 groups
  dt <- ref$dose_tests
  tt <- tw$tests
  m <- match(paste(dt$taxon, dt$time, dt$direction, dt$dose),
             paste(tt$taxon, tt$time, tt$direction, tt$concentration))
  expect_identical(tt$tested[m], dt$tested)
  both <- !is.na(dt$statistic)
  expect_equal(tt$statistic[m][both], dt$statistic[both])
  exact <- both & dt$dose %in% c(1, 10)
  expect_equal(tt$p_value[m][exact], dt$p_one_sided[exact])
  mc <- both & dt$dose == 3
  p <- dt$p_one_sided[mc]
  z <- (tt$p_value[m][mc] - p) / sqrt(2 * p * (1 - p) / 99999 + 1e-12)
  expect_lt(max(abs(z)), 4)
})

test_that("NOECs reproduce the original script on the azoxystrobin data", {
  ref <- readRDS(test_path("fixtures", "original_noec_azoxy.rds"))
  for (h in c("W", "S")) {
    hab <- if (h == "W") "eutrophic_lake" else "water_spring"
    d <- azoxy_macro(habitat = hab)
    tw <- test_control(d, method = "williams", step_variance = "reduced")
    s11 <- derive_noec(tw, pre_times = -7, criteria = s11_criteria("macroinvertebrates"))
    expect_same_noec(s11$noec, ref[[paste0(h, "_s11")]])
    if (h == "W") {
      expect_same_noec(derive_noec(tw, pre_times = -7)$noec, ref$W_williams)
    }
  }
})

test_that("the PRC reproduces the original script", {
  ref <- readRDS(test_path("fixtures", "original_prc_azoxy_W.rds"))
  prc <- fit_prc(azoxy_macro())
  r <- ref$prc_df[ref$prc_df$Treatment != "1", ]
  lev <- c("2" = "2", "3" = "20", "4" = "200")
  new <- prc$curves[match(paste(lev[as.character(r$Treatment)], r$Day),
                          paste(prc$curves$treatment, prc$curves$time)), ]
  expect_equal(new$cdt, r$Cdt)
  expect_equal(prc$weights$bk[match(names(ref$bk), prc$weights$taxon)], unname(ref$bk))
  expect_equal(prc$variance$axis1_pct_total, ref$pct_total)
  expect_equal(prc$variance$axis1_pct_treatment, ref$pct_treat)
  expect_setequal(prc$weights$taxon[abs(prc$weights$bk) >= 0.25], ref$taxa_plot)
  tp <- test_prc(prc, n_perm = 999, seed = 123)
  expect_equal(tp$global$F, ref$F_global)
  expect_identical(tp$global$p_value, ref$p_global)
  bt <- tp$by_time[match(as.numeric(ref$per_date$Day), tp$by_time$time), ]
  expect_identical(round(bt$F, 2), ref$per_date$F)
  expect_identical(bt$p_value, ref$per_date$p)
  expect_identical(round(100 * bt$r_squared, 1), ref$per_date$Var_expl_pct)
})
