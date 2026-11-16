# Complete workflow of a user with the example data of cosmet.
# Run with Rscript --vanilla in an empty folder; it is the basis of the README,
# the vignette and tests/testthat/test-workflow.R.

library(cosmet)

# 1. Data: macroinvertebrates of the eutrophic-lake mesocosms of the
#    azoxystrobin experiment (Amador et al. 2024); control and 2, 20 and
#    200 ug/L, three replicates, applications on days 0 and 14.
macro <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake")

d <- meso_data(
  macro,
  mesocosm_col = "tank", time_col = "day", treatment_col = "concentration",
  concentration_col = "concentration", taxon_col = "taxon", value_col = "abundance",
  control = "0", zero_fill = FALSE, unit = "individuals/sample",
  applications = c(0, 14)
)
d
check_design(d)

# 2. Taxa that can be evaluated (mean in the controls) and ln(2x + 1)
elig <- taxon_eligibility(d, min_control_mean = 3)
head(elig)
d <- transform_response(d, A = 2)

# 3. Williams' test by taxon and date, NOECs with the criteria of Text S11, MDD
tw <- test_control(d, method = "williams")
tw
noec <- derive_noec(tw, pre_times = -7, criteria = s11_criteria("macroinvertebrates"))
noec
noec_table(noec, selected = TRUE)
mdd <- calc_mdd(tw)
head(mdd[mdd$concentration == 200, c("taxon", "time", "mdd_pct", "mdd_class")])

# Dunnett's test, which does not assume a monotonic response
td <- test_control(d, method = "dunnett")
head(as.data.frame(derive_noec(td, pre_times = -7)))

# 4. Community: PRC, its test and community NOECs
prc <- fit_prc(d)
prc
prc_test <- test_prc(prc, n_perm = 199)
prc_test
cnoec <- community_noec(d, pre_times = -7)
noec_table(cnoec)

# 5. Plots (need ggplot2, patchwork and viridisLite) and export
if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE) &&
    requireNamespace("viridisLite", quietly = TRUE)) {
  p1 <- plot_prc(prc)
  p2 <- plot_noec(noec)
  p3 <- plot_timecourse(d, taxa = c("Chironomini", "Cloeon"))
  out <- tempfile(fileext = ".pdf")
  grDevices::pdf(out, width = 10, height = 6)
  print(p1)
  print(p2)
  print(p3)
  grDevices::dev.off()
  unlink(out)
}
utils::write.csv(noec_table(noec), tmp <- tempfile(fileext = ".csv"), row.names = FALSE)
unlink(tmp)
