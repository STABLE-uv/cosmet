# The flow of dev/usage_examples.R, with checks on each step.

test_that("the complete workflow runs on the example data", {
  macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")
  d <- meso_data(macro, mesocosm_col = "tank", time_col = "day",
                 treatment_col = "concentration", concentration_col = "concentration",
                 taxon_col = "taxon", value_col = "abundance", control = "0",
                 zero_fill = FALSE, unit = "individuals/sample", applications = c(0, 14))
  expect_true(check_design(d)$balanced)
  elig <- taxon_eligibility(d, min_control_mean = 3)
  expect_identical(nrow(elig), 29L)
  d <- transform_response(d, A = 2)

  tw <- test_control(d, method = "williams")
  expect_identical(nrow(tw$cells), 145L)
  expect_identical(nrow(tw$tests), 145L * 3L * 2L)
  noec <- derive_noec(tw, pre_times = -7, criteria = s11_criteria("macroinvertebrates"))
  sel <- noec_table(noec, selected = TRUE)
  # The same four taxa as the original analysis script (which used
  # step_variance = "reduced"; see test-regression.R and
  # dev/validation/williams_variants.md)
  expect_identical(sel$taxon, c("Chironomini", "Colymbetes fuscus (adult)",
                                "Colymbetes fuscus larvae", "Micronecta"))
  mdd <- calc_mdd(tw)
  expect_identical(nrow(mdd), 145L * 3L)
  expect_true(all(mdd$mdd_class[!is.na(mdd$mdd_pct)] %in% c("0", "I", "II", "III", "IV")))

  td <- test_control(d, method = "dunnett")
  expect_identical(nrow(td$tests), 145L * 3L)
  expect_s3_class(derive_noec(td, pre_times = -7), "meso_noec")

  prc <- fit_prc(d)
  expect_equal(prc$variance$axis1_pct_treatment, 34.6985, tolerance = 1e-5)
  tp <- test_prc(prc, n_perm = 199)
  expect_lt(tp$global$p_value, 0.05)
  cn <- community_noec(d, pre_times = -7)
  expect_identical(names(noec_table(cn)), c("taxon", "-7", "7", "21", "35", "56"))
})

test_that("print methods are stable", {
  d <- azoxy_macro(taxa = c("Chironomini", "Micronecta", "Cloeon"))
  tw <- test_control(d, method = "williams")
  expect_snapshot({
    print(tw)
    print(derive_noec(tw, pre_times = -7, criteria = s11_criteria("macroinvertebrates")))
  })
})
