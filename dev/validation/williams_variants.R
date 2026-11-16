# How much do NOECs change between the two versions of the step-down Williams
# test: the variance and isotonic means of all groups (Williams 1971, 1972;
# step_variance = "all", the default) or of the reduced design of each step
# (step_variance = "reduced", the original analysis script of the authors)?
# The script also checks step_variance = "all" against an independent
# implementation written here.
# Use, from the root of the package:
#   Rscript dev/validation/williams_variants.R
# Writes dev/validation/williams_variants.md.

pkgload::load_all(".", quiet = TRUE)

# Williams' null distribution for i doses with the df of the variance of all groups.
null_full <- function(n, df, n_sim, seed) {
  withr::with_seed(seed, {
    k <- length(n)
    z0 <- stats::rnorm(n_sim, sd = 1 / sqrt(n[1L]))
    ws <- numeric(n_sim); np <- 0; ymax <- rep(-Inf, n_sim)
    for (i in seq.int(k, 2L)) {
      ws <- ws + stats::rnorm(n_sim, sd = sqrt(n[i])); np <- np + n[i]
      ymax <- pmax(ymax, ws / np)
    }
    rs <- sqrt(stats::rchisq(n_sim, df = df) / df)
    sort((ymax - z0) / (rs * sqrt(1 / n[k] + 1 / n[1L])))
  })
}

# Isotonic (max-min) means of the treated groups, increasing.
iso_maxmin <- function(mu, n) {
  k <- length(mu)
  vapply(seq_len(k), function(i) {
    max(vapply(seq_len(i), function(u) {
      min(vapply(i:k, function(v) sum(n[u:v] * mu[u:v]) / sum(n[u:v]), numeric(1)))
    }, numeric(1)))
  }, numeric(1))
}

full_williams_noec <- function(d, alpha = 0.05, n_sim = 99999, seed = 123) {
  lev <- d$settings$treatments
  conc <- d$design$concentration[match(lev, d$design$treatment)]
  dat <- d$data
  cells <- unique(dat[c("taxon", "time")])
  cache <- new.env()
  res <- character(nrow(cells))
  for (ci in seq_len(nrow(cells))) {
    rows <- dat$taxon == cells$taxon[ci] & dat$time == cells$time[ci]
    g <- match(dat$treatment[rows], lev)
    prep <- prepare_cell(dat$response[rows], dat$value[rows], g, length(lev), FALSE)
    if (prep$status != "ok") { res[ci] <- NA_character_; next }
    k <- length(lev)
    df <- sum(prep$n - 1)
    s <- sqrt(sum(prep$ss) / df)
    found <- list()
    for (sign in c(-1, 1)) {
      m_iso <- iso_maxmin(sign * prep$mu[-1], prep$n[-1])
      rej <- logical(0)
      for (i in seq.int(k, 2L)) {
        stat <- (m_iso[i - 1] - sign * prep$mu[1]) / (s * sqrt(1 / prep$n[i] + 1 / prep$n[1]))
        if (i == 2L) {
          p <- stats::pt(stat, df, lower.tail = FALSE)
        } else {
          key <- paste(c(prep$n[seq_len(i)], df), collapse = "_")
          if (is.null(cache[[key]])) assign(key, null_full(prep$n[seq_len(i)], df, n_sim, seed), envir = cache)
          nl <- cache[[key]]
          p <- (1 + length(nl) - findInterval(stat, nl, left.open = TRUE)) / (length(nl) + 1)
        }
        r <- min(1, 2 * p) < alpha
        rej <- c(rej, stats::setNames(r, conc[i]))
        if (!r) break
      }
      found[[if (sign < 0) "decrease" else "increase"]] <- rej
    }
    top_rej <- vapply(found, function(r) isTRUE(r[1]), logical(1))
    if (sum(top_rej) == 0) {
      res[ci] <- paste0(">", format_conc(max(conc)))
    } else if (sum(top_rej) == 2) {
      res[ci] <- "ambiguous"
    } else {
      r <- found[[names(top_rej)[top_rej]]]
      loec <- min(as.numeric(names(r)[r]))
      res[ci] <- if (loec == min(conc[-1])) paste0("<", format_conc(loec)) else
        format_conc(conc[match(loec, conc) - 1L])
    }
  }
  cbind(cells, noec_full = res, stringsAsFactors = FALSE)
}

check_lines <- character()
sets <- list(
  list(community = "phytoplankton", A = 33.3, unit = "cells/mL"),
  list(community = "zooplankton", A = 10, unit = "individuals/L"),
  list(community = "macroinvertebrates", A = 2, unit = "individuals/sample")
)
out <- c("# Step-down Williams test: reduced design against the variance of all groups", "",
         sprintf("Written by `dev/validation/williams_variants.R` on %s. Data `azoxy`, ln(Ax + 1) with the A of Amador et al. (2024), alpha = 0.05, both directions (p x 2), 99 999 draws, seed 123.", format(Sys.Date())), "",
         "| Community | Habitat | Taxa x dates tested | Different NOECs | Detail (taxon, day: reduced -> all groups) |",
         "|---|---|---|---|---|")
for (st in sets) for (h in c("eutrophic_lake", "water_spring")) {
  x <- azoxy[azoxy$community == st$community & azoxy$habitat == h, ]
  d <- transform_response(meso_data(x, control = "0", zero_fill = FALSE, mesocosm_col = "tank",
                                    time_col = "day", treatment_col = "concentration",
                                    concentration_col = "concentration", value_col = "abundance",
                                    unit = st$unit), A = st$A)
  red <- derive_noec(test_control(d, method = "williams", step_variance = "reduced"))$noec
  pkg_all <- derive_noec(test_control(d, method = "williams", step_variance = "all"))$noec
  red$noec_red <- ifelse(red$status == "ambiguous_direction", "ambiguous", red$noec_label)
  full <- full_williams_noec(d)
  # Cross-check: step_variance = "all" against the independent code above
  pkg_lab <- ifelse(pkg_all$status == "ambiguous_direction", "ambiguous", pkg_all$noec_label)
  chk <- merge(data.frame(taxon = pkg_all$taxon, time = pkg_all$time, pkg = pkg_lab),
               full, by = c("taxon", "time"))
  n_check <- sum(!is.na(chk$noec_full))
  n_agree <- sum(!is.na(chk$noec_full) & chk$pkg %in% chk$noec_full & chk$pkg == chk$noec_full)
  check_lines <<- c(check_lines, sprintf("- %s, %s: %d of %d NOECs equal", st$community, h, n_agree, n_check))
  m <- merge(red[c("taxon", "time", "status", "noec_red")], full, by = c("taxon", "time"))
  m <- m[m$status %in% c("ok", "ambiguous_direction", "undefined_step"), ]
  diff <- m[!((!is.na(m$noec_red) & !is.na(m$noec_full) & m$noec_red == m$noec_full) |
                (is.na(m$noec_red) & is.na(m$noec_full))), ]
  detail <- if (nrow(diff) == 0) "-" else paste(sprintf("%s, %g: %s -> %s", diff$taxon, diff$time,
                                                         ifelse(is.na(diff$noec_red), "NA", diff$noec_red),
                                                         ifelse(is.na(diff$noec_full), "NA", diff$noec_full)),
                                                 collapse = "; ")
  out <- c(out, sprintf("| %s | %s | %d | %d | %s |", st$community, h, nrow(m), nrow(diff), detail))
}
out <- c(out, "", "`test_control(step_variance = \"all\")` against the independent implementation of this script (equal NOECs / taxa x dates tested):", "", check_lines)
writeLines(out, file.path("dev", "validation", "williams_variants.md"))
cat(out, sep = "\n")
