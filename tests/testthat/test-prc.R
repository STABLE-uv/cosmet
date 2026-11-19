test_that("fit_prc() returns the coefficients and weights of vegan::prc()", {
  d <- azoxy_macro()
  prc <- fit_prc(d, orientation = "none")
  m <- community_matrix(d, NULL)
  trt <- factor(m$samples$treatment, levels = c("0", "2", "20", "200"))
  time <- factor(m$samples$time, levels = sort(unique(m$samples$time)))
  ref <- summary(vegan::prc(m$y, trt, time), scaling = "symmetric")
  cdt <- prc$curves[prc$curves$treatment != "0", ]
  expect_equal(cdt$cdt, as.vector(ref$coefficients))
  expect_equal(prc$weights$bk, unname(drop(ref$sp)[prc$weights$taxon]))
  expect_identical(prc$curves$cdt[prc$curves$treatment == "0"], rep(0, 5))
  v <- prc$variance
  expect_equal(v$time + v$treatment + v$residual, v$total)
})

test_that("b_k c_dt is the difference from the control fitted by the first axis", {
  d <- azoxy_macro()
  # The fitted first axis computed directly: changes over time removed,
  # treatment x time fitted by least squares, rank-1 part of the fit.
  m <- community_matrix(d, NULL)
  tm <- factor(m$samples$time)
  trt <- factor(m$samples$treatment, levels = d$settings$treatments)
  z <- stats::model.matrix(~ tm - 1)
  x <- stats::model.matrix(~ trt:tm - 1)
  residual <- function(a) a - z %*% qr.solve(z, a)
  fit <- qr.fitted(qr(residual(x)), residual(m$y))
  s <- svd(fit)
  rank1 <- s$d[1] * s$u[, 1, drop = FALSE] %*% t(s$v[, 1, drop = FALSE])
  key <- paste(m$samples$treatment, m$samples$time)
  means <- rowsum(rank1, key) / as.vector(table(key)[sort(unique(key))])
  for (sc in c("symmetric", "sites", "species")) {
    prc <- fit_prc(d, scaling = sc)
    cur <- prc$curves[prc$curves$treatment != d$settings$control, ]
    bk <- prc$weights$bk[match(colnames(m$y), prc$weights$taxon)]
    effect <- means[paste(cur$treatment, cur$time), , drop = FALSE] -
      means[paste(d$settings$control, cur$time), , drop = FALSE]
    expect_equal(unname(outer(cur$cdt, bk)), unname(effect), tolerance = 1e-8)
  }
})

test_that("the orientation puts the highest treatment downwards", {
  d <- azoxy_macro()
  up <- fit_prc(d, orientation = "none")
  down <- fit_prc(d)
  hi <- down$curves[down$curves$treatment == "200", ]
  expect_lt(hi$cdt[which.max(abs(hi$cdt))], 0)
  sgn <- if (down$settings$reversed) -1 else 1
  expect_equal(down$weights$bk, sgn * up$weights$bk)
})

test_that("multivariate analyses need complete tables", {
  x <- make_presence()
  d <- meso_data(x, control = "ctrl", zero_fill = FALSE)
  expect_error(fit_prc(d), "zero_fill")
  expect_error(community_noec(d), "zero_fill")
  one_time <- azoxy_macro()
  one_time <- meso_data(subset(azoxy, community == "macroinvertebrates" & day == 7 &
                                 habitat == "eutrophic_lake"),
                        control = "0", zero_fill = FALSE, mesocosm_col = "tank",
                        time_col = "day", treatment_col = "concentration",
                        concentration_col = "concentration", value_col = "abundance")
  expect_error(fit_prc(one_time), "two sampling times")
})

test_that("test_prc() permutes whole time series and is reproducible", {
  prc <- fit_prc(azoxy_macro())
  a <- test_prc(prc, n_perm = 99, seed = 1)
  b <- test_prc(prc, n_perm = 99, seed = 1)
  expect_identical(a$global, b$global)
  expect_identical(a$by_time, b$by_time)
  expect_identical(nrow(a$by_time), 5L)
  expect_true(a$global$p_value >= 1 / 100)
  set.seed(2)
  r1 <- stats::runif(1)
  set.seed(2)
  invisible(test_prc(prc, n_perm = 99, by_time = FALSE))
  expect_identical(stats::runif(1), r1)
})

test_that("test_prc() refuses unbalanced designs", {
  x <- subset(azoxy, community == "macroinvertebrates" & habitat == "eutrophic_lake" &
                !(tank == "W-0-1" & day == 56))
  d <- meso_data(x, control = "0", zero_fill = FALSE, mesocosm_col = "tank", time_col = "day",
                 treatment_col = "concentration", concentration_col = "concentration",
                 value_col = "abundance")
  prc <- fit_prc(transform_response(d, A = 2))
  expect_error(test_prc(prc, n_perm = 9), "every mesocosm")
})

test_that("community_noec() tests the first principal component of each time", {
  cn <- community_noec(azoxy_macro(), pre_times = -7, n_sim = 9999)
  expect_s3_class(cn, "meso_noec")
  expect_identical(unique(cn$noec$taxon), "PC1")
  expect_identical(cn$pca$time, c(-7, 7, 21, 35, 56))
  expect_true(all(cn$pca$pc1_pct > 0 & cn$pca$pc1_pct < 100))
  expect_true(all(cn$noec$direction %in% c(NA, "higher_scores", "lower_scores")))
  # Scores of a principal component have an arbitrary sign: no arrows
  expect_false(any(grepl("\u2191|\u2193", unlist(noec_table(cn)[-1]))))
  # The test on PC1 equals Williams' test on the scores computed by hand
  m <- community_matrix(azoxy_macro(), NULL)
  idx <- m$samples$time == 56
  pc1 <- vegan::scores(vegan::rda(m$y[idx, apply(m$y[idx, ], 2, stats::var) > 0]),
                       display = "sites", choices = 1, scaling = "sites")
  s <- data.frame(mesocosm = m$samples$mesocosm[idx], time = 56,
                  concentration = m$samples$concentration[idx], taxon = "PC1",
                  value = as.vector(pc1))
  dd <- meso_data(s, control = "0", zero_fill = FALSE, treatment_col = "concentration",
                  concentration_col = "concentration")
  ref <- derive_noec(test_control(dd, method = "williams", n_sim = 9999))
  expect_identical(cn$noec$noec_label[cn$noec$time == 56], ref$noec$noec_label)
})

test_that("print methods are stable", {
  prc <- fit_prc(azoxy_macro())
  expect_snapshot(print(prc))
})
