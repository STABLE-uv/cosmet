test_that("with one treatment Dunnett's test is Student's t test", {
  lambda <- sqrt(3 / 6)
  for (df in c(4, 10, 30)) {
    expect_lt(abs(dunnett_critical(lambda, df, 0.05, two_sided = TRUE) - stats::qt(0.975, df)),
              1e-4)
    expect_lt(abs(dunnett_critical(lambda, df, 0.05, two_sided = FALSE) - stats::qt(0.95, df)),
              1e-4)
    expect_lt(abs(dunnett_p(2.3, lambda, df, two_sided = TRUE) - 2 * stats::pt(-2.3, df)), 1e-5)
  }
})

test_that("Dunnett probabilities agree with mvtnorm", {
  skip_if_not_installed("mvtnorm")
  n <- c(4, 3, 3, 5, 3)  # control first, unequal replication
  lambda <- sqrt(n[-1] / (n[-1] + n[1]))
  df <- sum(n - 1)
  r <- outer(lambda, lambda)
  diag(r) <- 1
  alg <- mvtnorm::GenzBretz(maxpts = 2e5, abseps = 1e-6, releps = 0)
  for (cc in c(1.5, 2.2, 2.8)) {
    two <- withr::with_seed(1, mvtnorm::pmvt(lower = rep(-cc, 4), upper = rep(cc, 4), df = df,
                                              corr = r, algorithm = alg))
    one <- withr::with_seed(1, mvtnorm::pmvt(lower = rep(-Inf, 4), upper = rep(cc, 4), df = df,
                                              corr = r, algorithm = alg))
    # Within three times the error that mvtnorm reports for its own estimate
    expect_lt(abs(dunnett_prob(cc, lambda, df, TRUE) - two), 3 * attr(two, "error") + 1e-8)
    expect_lt(abs(dunnett_prob(cc, lambda, df, FALSE) - one), 3 * attr(one, "error") + 1e-8)
  }
})

test_that("Dunnett critical values have the stated coverage under mvtnorm", {
  skip_if_not_installed("mvtnorm")
  r <- matrix(0.5, 4, 4)
  diag(r) <- 1
  crit <- dunnett_critical(rep(sqrt(0.5), 4), 20, 0.05, TRUE)
  cover <- withr::with_seed(1, mvtnorm::pmvt(lower = rep(-crit, 4), upper = rep(crit, 4),
                                             df = 20, corr = r,
                                             algorithm = mvtnorm::GenzBretz(maxpts = 2e5,
                                                                            abseps = 1e-6)))
  expect_lt(abs(as.numeric(cover) - 0.95), 1e-5)
})

test_that("test_control() with Dunnett: p-values, intervals and decisions agree", {
  d <- meso_monotone()
  td <- test_control(d, method = "dunnett")
  tt <- td$tests
  expect_identical(tt$direction, rep("both", 3))
  expect_identical(tt$ci_type, rep("simultaneous", 3))
  # Two-sided: significant exactly when the interval excludes zero
  expect_identical(tt$significant, tt$lower > 0 | tt$upper < 0)
  # Same statistic as a pooled t statistic
  ss <- 2 + sum((c(9, 13, 10) - 32 / 3)^2) + 2 + 2
  s <- sqrt(ss / 8)
  expect_equal(tt$statistic, (c(32 / 3, 7, 2) - 11) / (s * sqrt(2 / 3)))
  expect_identical(tt$df, rep(8, 3))
  # One-sided decreases: lower bound -Inf
  dec <- test_control(d, method = "dunnett", direction = "decrease")$tests
  expect_true(all(dec$lower == -Inf))
  expect_true(all(dec$p_value <= tt$p_value + 1e-12))
})

test_that("Dunnett adjusted p-values match mvtnorm in a real case", {
  skip_if_not_installed("mvtnorm")
  d <- azoxy_macro(taxa = "Chironomini")
  td <- test_control(d, method = "dunnett")
  row <- td$tests[td$tests$time == 56, ]
  n <- c(row$n_control[1], row$n_treatment)
  lambda <- sqrt(n[-1] / (n[-1] + n[1]))
  r <- outer(lambda, lambda)
  diag(r) <- 1
  ref <- vapply(abs(row$statistic), function(t) {
    1 - withr::with_seed(1, mvtnorm::pmvt(lower = rep(-t, 3), upper = rep(t, 3), df = row$df[1],
                                          corr = r,
                                          algorithm = mvtnorm::GenzBretz(maxpts = 2e5,
                                                                         abseps = 1e-6)))
  }, numeric(1))
  expect_lt(max(abs(row$p_value - ref)), 2e-5)
})
