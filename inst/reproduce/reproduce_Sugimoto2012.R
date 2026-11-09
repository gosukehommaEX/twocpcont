# =============================================================================
# reproduce_Sugimoto2012.R
#
# Reproduction of the published values of Sugimoto, Sozu and Hamasaki
# (2012), Pharmaceutical Statistics 11(2): 118-128 (doi:10.1002/pst.505):
#   (a) Supporting Information, Appendix A.2.1: the constant
#       C = 0.9988124 computed by their Newton-Raphson algorithm for
#       K = 2, gamma = 8/7, rho = 0.5, alpha = 0.025 and 1 - beta = 0.8,
#       and the sample size n = 2 (C + z_0.025)^2 / 0.35^2 = 142.93 -> 143
#       for delta = (0.4, 0.35) and p = 1, reproduced with newton_ck() and
#       twocpcont_ss(method = "newton");
#   (b) Table I: the sequential search for delta = (0.55, 0.50),
#       rho = 0.5, alpha = 0.025 and p = 1 at 1 - beta = 0.8 and 0.9, i.e.
#       the starting points n0 of their (4), the co-primary power at each
#       candidate n and the resulting n, reproduced with twocpcont_power()
#       and twocpcont_ss(method = "bivariate").
# The exact root of their equation (6) is listed for information (INFO):
# it differs from the published C in the seventh decimal, because the
# Newton-Raphson algorithm stops once the power is within 1e-5 above the
# target.
#
# Verdicts (one row per published value)
#   PASS  the value of the package, rounded to the decimals of the
#         publication, equals the published value
#   INFO  shown for comparison only
#   FAIL  otherwise
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file:
#     devtools::load_all()
#     source("inst/reproduce/reproduce_Sugimoto2012.R")
#
# Output  inst/reproduce/results/reproduce_Sugimoto2012.csv
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

alpha   <- 0.025
z_alpha <- qnorm(1 - alpha)

rows <- list()
add_row <- function(item, published, digits, package, route,
                    info = FALSE) {
  shown <- round(package, digits)
  verdict <- if (info) "INFO" else
    if (abs(shown - published) < 1e-9) "PASS" else "FAIL"
  rows[[length(rows) + 1L]] <<- data.frame(
    item = item, route = route, published = published, digits = digits,
    package = package, difference = shown - published, verdict = verdict,
    stringsAsFactors = FALSE)
}

# =============================================================================
# (a) Appendix A.2.1: gamma = 8/7, rho = 0.5, alpha = 0.025, 1 - beta = 0.8
# =============================================================================
ck <- newton_ck(kappa = 8 / 7, rho = 0.5, alpha = alpha, beta = 0.2)
add_row("Appendix A.2.1: C (gamma = 8/7, rho = 0.5, 1 - beta = 0.80)",
        0.9988124, 7, ck$C, "newton_ck()")
C_exact <- uniroot(function(cc) {
  pbivnorm::pbivnorm(8 / 7 * cc + z_alpha * (8 / 7 - 1), cc, rho = 0.5) - 0.8
}, lower = -2, upper = 5, tol = 1e-13)$root
add_row("Appendix A.2.1: C (gamma = 8/7, rho = 0.5, 1 - beta = 0.80)",
        0.9988124, 7, C_exact, "exact root of (6) by uniroot()", info = TRUE)

ss <- twocpcont_ss(delta1 = 0.4, delta2 = 0.35, sd1 = 1, sd2 = 1, rho = 0.5,
                   r = 1, alpha = alpha, beta = 0.2, method = "newton")
add_row("Appendix A.2.1: n = 2 (C + z_0.025)^2 / 0.35^2, delta = (0.4, 0.35)",
        142.93, 2, ss$n2_cont, "twocpcont_ss(method = \"newton\")$n2_cont")
add_row("Appendix A.2.1: n rounded up", 143, 0, ss$n2,
        "twocpcont_ss(method = \"newton\")$n2")

# =============================================================================
# (b) Table I: delta = (0.55, 0.50), rho = 0.5, alpha = 0.025, p = 1
# =============================================================================
table1 <- list(
  list(beta = 0.2, n0 = 62.8, n = 63:72,
       power = c(0.734, 0.742, 0.750, 0.758, 0.765, 0.773, 0.780, 0.787,
                 0.794, 0.8001),
       digits = c(rep(3, 9), 4)),
  list(beta = 0.1, n0 = 84.1, n = 85:93,
       power = c(0.871, 0.875, 0.879, 0.883, 0.888, 0.891, 0.895, 0.899,
                 0.902),
       digits = rep(3, 9))
)
for (tb in table1) {
  # Starting point (4): the single-endpoint sample size with the smaller
  # effect size, n0 = (1 / kappa_p) (z_beta + z_alpha)^2 / 0.50^2 with
  # kappa_p = p / (1 + p) = 1 / 2
  add_row(sprintf("Table I: starting point n0, 1 - beta = %.2f", 1 - tb$beta),
          tb$n0, 1, 2 * (qnorm(1 - tb$beta) + z_alpha) ^ 2 / 0.50 ^ 2,
          "(4) of Sugimoto et al.")
  for (i in seq_along(tb$n)) {
    pw <- twocpcont_power(n1 = tb$n[i], n2 = tb$n[i], delta1 = 0.55,
                          delta2 = 0.50, sd1 = 1, sd2 = 1, rho = 0.5,
                          alpha = alpha)$powerCoprimary
    add_row(sprintf("Table I: power at n = %d, 1 - beta = %.2f", tb$n[i],
                    1 - tb$beta),
            tb$power[i], tb$digits[i], pw, "twocpcont_power()")
  }
  ss_b <- twocpcont_ss(delta1 = 0.55, delta2 = 0.50, sd1 = 1, sd2 = 1,
                       rho = 0.5, r = 1, alpha = alpha, beta = tb$beta,
                       method = "bivariate")
  add_row(sprintf("Table I: required n, 1 - beta = %.2f", 1 - tb$beta),
          max(tb$n), 0, ss_b$n2, "twocpcont_ss(method = \"bivariate\")$n2")
}

res <- do.call(rbind, rows)
out_file <- file.path(results_dir, "reproduce_Sugimoto2012.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf(
  "reproduce_Sugimoto2012.R: PASS %d, INFO %d, FAIL %d; written to %s",
  sum(res$verdict == "PASS"), sum(res$verdict == "INFO"),
  sum(res$verdict == "FAIL"), out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("reproduce_Sugimoto2012: %d FAIL verdict(s); see %s.",
               sum(res$verdict == "FAIL"), out_file), call. = FALSE)
}
