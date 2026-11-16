# Validation of the numerical core

Written by `dev/validation/validate_core.R` on 2026-11-16 with R 4.3.3.

## 1. Null distribution of Williams' statistic against the table of Williams (1972)

23 combinations of replication (n per group) and number of doses with df in the table; 1 000 000 draws each. Monte Carlo error of the 0.95 quantile: about 0.003.

| n | doses | df | simulated | table | difference |
|---|---|---|---|---|---|
| 2 | 4 | 5 | 2.2142 | 2.209 | +0.0052 |
| 2 | 5 | 6 | 2.1282 | 2.131 | -0.0028 |
| 2 | 6 | 7 | 2.0787 | 2.076 | +0.0027 |
| 2 | 8 | 9 | 2.0058 | 2.006 | -0.0002 |
| 3 | 2 | 6 | 2.0590 | 2.058 | +0.0010 |
| 3 | 3 | 8 | 1.9954 | 1.997 | -0.0016 |
| 3 | 4 | 10 | 1.9586 | 1.956 | +0.0026 |
| 3 | 5 | 12 | 1.9267 | 1.927 | -0.0003 |
| 3 | 6 | 14 | 1.9054 | 1.906 | -0.0006 |
| 3 | 8 | 18 | 1.8764 | 1.878 | -0.0016 |
| 4 | 2 | 9 | 1.9317 | 1.931 | +0.0007 |
| 4 | 3 | 12 | 1.9031 | 1.903 | +0.0001 |
| 4 | 4 | 15 | 1.8821 | 1.882 | +0.0001 |
| 4 | 5 | 18 | 1.8651 | 1.867 | -0.0019 |
| 5 | 2 | 12 | 1.8736 | 1.873 | +0.0006 |
| 5 | 3 | 16 | 1.8612 | 1.860 | +0.0012 |
| 5 | 4 | 20 | 1.8479 | 1.847 | +0.0009 |
| 5 | 5 | 24 | 1.8378 | 1.838 | -0.0002 |
| 5 | 6 | 28 | 1.8313 | 1.831 | +0.0003 |
| 6 | 2 | 15 | 1.8385 | 1.840 | -0.0015 |
| 6 | 3 | 20 | 1.8359 | 1.834 | +0.0019 |
| 6 | 5 | 30 | 1.8225 | 1.821 | +0.0015 |
| 6 | 6 | 35 | 1.8171 | 1.816 | +0.0011 |

Largest absolute difference: 0.0052; mean: +0.0004. The table has 3 decimals.

## 2. Dunnett probabilities

Quadrature of `dunnett_prob()` against nested numerical integration (`stats::integrate()`, relative tolerance 1e-10): 400 cases (df 2 to 100; 1 to 5 treatments; equal and unequal replication; c from 1 to 5; one- and two-sided). Largest absolute error: 5.9e-11; with P > 0.9 (where small p-values are), 3.6e-12.

Against `mvtnorm::pmvt()` (two-sided, Genz-Bretz with a stated absolute error below 1e-6): 12 cases, largest absolute difference 3e-07.

Critical values with 4 treatments, equal replication and 20 df: 2.6510 (two-sided) and 2.3044 (one-sided).

