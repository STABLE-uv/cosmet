# Broad validation of the numerical core (outside R CMD check).
# Use, from the root of the package:
#   PMCMRPLUS_SYSDATA=<path> Rscript dev/validation/validate_core.R
# Writes dev/validation/validate_core.md.
#
# 1. Simulated null distribution of Williams' statistic against the table of
#    Williams (1972) (one-sided alpha = 0.05, equal replication), as
#    transcribed in PMCMRplus (TabCrit$williams.tk005 in R/sysdata.rda of the
#    CRAN package, read from github.com/cran/PMCMRplus). PMCMRPLUS_SYSDATA
#    holds the path to that file.
# 2. Dunnett probabilities by quadrature against nested numerical integration
#    with stats::integrate() (relative tolerance 1e-10) and against mvtnorm.

pkgload::load_all(".", quiet = TRUE)
out <- c("# Validation of the numerical core", "",
         sprintf("Written by `dev/validation/validate_core.R` on %s with R %s.",
                 format(Sys.Date()), getRversion()), "")

# 1. Williams ------------------------------------------------------------------
sysdata <- Sys.getenv("PMCMRPLUS_SYSDATA")
if (!nzchar(sysdata) || !file.exists(sysdata)) {
  stop("Set PMCMRPLUS_SYSDATA to the path of R/sysdata.rda of the PMCMRplus sources.")
}
e <- new.env()
load(sysdata, envir = e)
tab <- e$TabCrit$williams.tk005
n_sim <- 1e6
res <- NULL
for (nrep in 2:6) for (k in c(2, 3, 4, 5, 6, 8)) {
  df <- (k + 1) * (nrep - 1)
  if (!as.character(df) %in% rownames(tab)) next
  q <- withr::with_seed(2026, stats::quantile(williams_null(rep(nrep, k + 1), n_sim), 0.95,
                                               names = FALSE))
  res <- rbind(res, data.frame(n = nrep, doses = k, df = df, simulated = q,
                               table = tab[as.character(df), as.character(k)]))
}
res$difference <- res$simulated - res$table
out <- c(out, "## 1. Null distribution of Williams' statistic against the table of Williams (1972)", "",
         sprintf("%d combinations of replication (n per group) and number of doses with df in the table; %s draws each. Monte Carlo error of the 0.95 quantile: about 0.003.",
                 nrow(res), format(n_sim, big.mark = " ", scientific = FALSE)), "",
         "| n | doses | df | simulated | table | difference |", "|---|---|---|---|---|---|",
         sprintf("| %d | %d | %d | %.4f | %.3f | %+.4f |", res$n, res$doses, res$df, res$simulated,
                 res$table, res$difference), "",
         sprintf("Largest absolute difference: %.4f; mean: %+.4f. The table has 3 decimals.",
                 max(abs(res$difference)), mean(res$difference)), "")

# 2. Dunnett -------------------------------------------------------------------
nested <- function(c, lambda, df, two_sided) {
  s <- sqrt(1 - lambda^2)
  inner <- function(v) vapply(v, function(vv) stats::integrate(function(z) {
    o <- stats::dnorm(z)
    for (i in seq_along(lambda)) {
      up <- stats::pnorm((c * vv + lambda[i] * z) / s[i])
      if (two_sided) up <- up - stats::pnorm((-c * vv + lambda[i] * z) / s[i])
      o <- o * up
    }
    o
  }, -Inf, Inf, rel.tol = 1e-10, abs.tol = 0)$value, numeric(1))
  logf <- function(v) (df / 2) * log(df) - (df / 2 - 1) * log(2) - lgamma(df / 2) +
    (df - 1) * log(v) - df * v^2 / 2
  stats::integrate(function(v) inner(v) * exp(logf(v)), 0, Inf, rel.tol = 1e-10, abs.tol = 0)$value
}
grid <- NULL
for (df in c(2, 3, 5, 8, 12, 20, 40, 100)) for (pat in list(c(3, 3), c(3, 3, 3, 3), c(2, 10, 10),
                                                          c(10, 2, 2, 2, 2, 2), c(4, 3, 3, 5, 3))) {
  lambda <- sqrt(pat[-1] / (pat[-1] + pat[1]))
  for (cc in c(1, 2, 2.5, 3.5, 5)) for (ts in c(TRUE, FALSE)) {
    a <- nested(cc, lambda, df, ts)
    b <- dunnett_prob(cc, lambda, df, ts)
    grid <- rbind(grid, data.frame(df = df, pattern = paste(pat, collapse = "-"), c = cc,
                                   two_sided = ts, nested = a, quadrature = b))
  }
}
grid$error <- grid$quadrature - grid$nested
mv <- NULL
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  for (pat in list(c(3, 3, 3, 3), c(4, 3, 3, 5, 3), c(10, 2, 2, 2, 2, 2))) {
    lambda <- sqrt(pat[-1] / (pat[-1] + pat[1]))
    r <- outer(lambda, lambda)
    diag(r) <- 1
    df <- sum(pat - 1)
    for (cc in c(1.5, 2.2, 2.8, 3.5)) {
      p <- withr::with_seed(1, mvtnorm::pmvt(lower = rep(-cc, length(lambda)),
                                             upper = rep(cc, length(lambda)), df = df, corr = r,
                                             algorithm = mvtnorm::GenzBretz(maxpts = 2e6, abseps = 1e-7)))
      mv <- rbind(mv, data.frame(pattern = paste(pat, collapse = "-"), df = df, c = cc,
                                 mvtnorm = as.numeric(p), mvtnorm_error = attr(p, "error"),
                                 quadrature = dunnett_prob(cc, lambda, df, TRUE)))
    }
  }
  mv$difference <- mv$quadrature - mv$mvtnorm
}
out <- c(out, "## 2. Dunnett probabilities", "",
         sprintf("Quadrature of `dunnett_prob()` against nested numerical integration (`stats::integrate()`, relative tolerance 1e-10): %d cases (df 2 to 100; 1 to 5 treatments; equal and unequal replication; c from 1 to 5; one- and two-sided). Largest absolute error: %.2g; with P > 0.9 (where small p-values are), %.2g.",
                 nrow(grid), max(abs(grid$error)), max(abs(grid$error[grid$nested > 0.9]))), "")
if (!is.null(mv)) {
  out <- c(out, sprintf("Against `mvtnorm::pmvt()` (two-sided, Genz-Bretz with a stated absolute error below 1e-6): %d cases, largest absolute difference %.2g.",
                        nrow(mv), max(abs(mv$difference))), "")
}
crit <- dunnett_critical(rep(sqrt(0.5), 4), 20, 0.05, TRUE)
crit1 <- dunnett_critical(rep(sqrt(0.5), 4), 20, 0.05, FALSE)
out <- c(out, sprintf("Critical values with 4 treatments, equal replication and 20 df: %.4f (two-sided) and %.4f (one-sided).", crit, crit1), "")
writeLines(out, file.path("dev", "validation", "validate_core.md"))
cat(out, sep = "\n")
