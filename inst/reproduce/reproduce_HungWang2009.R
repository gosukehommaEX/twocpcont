# =============================================================================
# reproduce_HungWang2009.R
#
# Reproduction of the sample size ratio of Hung and Wang (2009, Section 3
# and Figure 1) for two co-primary endpoints: the ratio n / m1 of the sample
# size for two independent endpoints with equal standardized effect sizes to
# the single-endpoint sample size,
#   n / m1 = {(z_alpha + z_p) / (z_alpha + z_beta)}^2,  z_p = qnorm(sqrt(1 - beta)),
# at alpha = 0.025 (one-sided) and 1 - beta = 0.80, 0.85 and 0.90. In the
# twocpcont article this ratio equals 1 / {1 - R_max(1, beta)}, which is
# computed here with r_max() by both of its routes (exact = TRUE, the default,
# and the closed-form approximation exact = FALSE, which coincide at
# kappa = 1).
#
# Verdicts (one row per target power and route)
#   PASS  |package - published| < 1e-10
#   FAIL  otherwise
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file:
#     devtools::load_all()
#     source("inst/reproduce/reproduce_HungWang2009.R")
#
# Output  inst/reproduce/results/reproduce_HungWang2009.csv
# The script stops with an error if any verdict is FAIL.
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "twocpcont")) {
  stop("Run this script from the root of the twocpcont package ",
       "(open twocpcont.Rproj first).", call. = FALSE)
}
if (!requireNamespace("twocpcont", quietly = TRUE)) {
  stop("The twocpcont package is neither loaded nor installed.", call. = FALSE)
}
if (!"package:twocpcont" %in% search()) {
  library(twocpcont)
}

results_dir <- file.path("inst", "reproduce", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

alpha     <- 0.025
power_set <- c(0.80, 0.85, 0.90)
z_alpha   <- qnorm(1 - alpha)

ratio_hw <- ((z_alpha + qnorm(sqrt(power_set))) /
               (z_alpha + qnorm(power_set))) ^ 2

res <- do.call(rbind, lapply(c(TRUE, FALSE), function(ex) {
  pkg <- sapply(power_set, function(pw) {
    1 / (1 - r_max(kappa = 1, alpha = alpha, beta = 1 - pw, exact = ex))
  })
  data.frame(item = sprintf("1 - beta = %.2f", power_set),
             route = if (ex) "r_max(exact = TRUE)" else "r_max(exact = FALSE)",
             published = ratio_hw, package = pkg,
             difference = pkg - ratio_hw,
             verdict = ifelse(abs(pkg - ratio_hw) < 1e-10, "PASS", "FAIL"),
             stringsAsFactors = FALSE)
}))

out_file <- file.path(results_dir, "reproduce_HungWang2009.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf("reproduce_HungWang2009.R: PASS %d, FAIL %d; written to %s",
                sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"),
                out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("reproduce_HungWang2009: %d FAIL verdict(s); see %s.",
               sum(res$verdict == "FAIL"), out_file), call. = FALSE)
}
