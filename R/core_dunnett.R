# Numerical core of Dunnett's test (internal) ----------------------------------
#
# Dunnett CW (1955) JASA 50:1096-1121. With treatment means, a control mean
# and a pooled variance with `df` degrees of freedom, the statistics
# T_i = (ybar_i - ybar_0) / (s sqrt(1/n_i + 1/n_0)) are multivariate t with
# correlations lambda_i lambda_j, lambda_i = sqrt(n_i / (n_i + n_0)).
# Conditional on the control mean (U_0 = z) and on V = s / sigma, they are
# independent, which turns the probability into a two-dimensional integral:
#   P(max_i |T_i| < c) = E_V E_Z prod_i [Phi((c V + lambda_i Z) / sqrt(1 - lambda_i^2))
#                                        - Phi((-c V + lambda_i Z) / sqrt(1 - lambda_i^2))]
# (one-sided: only the first term). The integral is computed with fixed
# nodes, so the result is deterministic: Gauss-Hermite in Z and the
# tanh-sinh (double exponential) rule on the probability scale of V,
# u = F_V(v), which copes with the singular behaviour of the quantile function
# at 0 and 1. Against nested stats::integrate() with relative tolerance 1e-10
# the absolute error was below 1e-10 for df from 2 to 100 and 1 to 5
# treatments (dev/validation/validate_core.R).

.quadrature <- new.env(parent = emptyenv())

# Tanh-sinh nodes and weights on (0, 1): x = (1 + tanh(pi/2 sinh(t))) / 2.
tanh_sinh_01 <- function(h = 1 / 12, t_max = 3.2) {
  t <- seq(-t_max, t_max, by = h)
  s <- (pi / 2) * sinh(t)
  x <- (1 + tanh(s)) / 2
  w <- h * (pi / 4) * cosh(t) / cosh(s)^2
  keep <- x > 0 & x < 1
  list(x = x[keep], w = w[keep])
}

# Probabilists' Gauss-Hermite: sum(w * f(x)) approximates E[f(Z)], Z ~ N(0, 1).
gauss_hermite_prob <- function(m) {
  i <- seq_len(m - 1L)
  b <- sqrt(i)
  jac <- matrix(0, m, m)
  jac[cbind(i, i + 1L)] <- b
  jac[cbind(i + 1L, i)] <- b
  e <- eigen(jac, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

quadrature_nodes <- function() {
  if (is.null(.quadrature$ts)) {
    assign("ts", tanh_sinh_01(), envir = .quadrature)
    assign("gh", gauss_hermite_prob(96L), envir = .quadrature)
  }
  list(ts = .quadrature$ts, gh = .quadrature$gh)
}

# P(max T_i < c) (one-sided) or P(max |T_i| < c) (two-sided) for each c.
dunnett_prob <- function(c, lambda, df, two_sided) {
  q <- quadrature_nodes()
  v <- sqrt(stats::qchisq(q$ts$x, df) / df)
  s <- sqrt(1 - lambda^2)
  vapply(c, function(ci) {
    cv <- matrix(ci * v, nrow = length(q$gh$x), ncol = length(v), byrow = TRUE)
    prod_mat <- matrix(1, nrow = length(q$gh$x), ncol = length(v))
    for (i in seq_along(lambda)) {
      lz <- lambda[i] * q$gh$x
      up <- stats::pnorm((cv + lz) / s[i])
      if (two_sided) {
        up <- up - stats::pnorm((-cv + lz) / s[i])
      }
      prod_mat <- prod_mat * up
    }
    sum(q$gh$w * (prod_mat %*% q$ts$w))
  }, numeric(1))
}

# Single-step adjusted p-values. `t` is oriented so that large values are
# evidence for the alternative (for two-sided tests its absolute value is used).
dunnett_p <- function(t, lambda, df, two_sided) {
  stat <- if (two_sided) abs(t) else t
  p <- 1 - dunnett_prob(stat, lambda, df, two_sided)
  pmin(1, pmax(0, p))
}

# Critical value c with P(max T < c) = 1 - alpha (or with |T|).
dunnett_critical <- function(lambda, df, alpha, two_sided, cache = NULL) {
  key <- paste(c(signif(lambda, 12), df, alpha, two_sided), collapse = "_")
  if (!is.null(cache)) {
    out <- cache[[key]]
    if (!is.null(out)) return(out)
  }
  f <- function(cc) dunnett_prob(cc, lambda, df, two_sided) - (1 - alpha)
  out <- stats::uniroot(f, lower = 0, upper = 10, extendInt = "upX", tol = 1e-9)$root
  if (!is.null(cache)) assign(key, out, envir = cache)
  out
}
