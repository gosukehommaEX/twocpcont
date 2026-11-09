# =============================================================================
# reproduce_Sozu2015.R
#
# Reproduction of the constant C_2 of the convenient formula of Sugimoto et
# al. (2012), as tabulated by Sozu et al. (2015, Tables 4.3 and 4.4) for
# K = 2, alpha = 0.025 (one-sided), 1 - beta = 0.80 and 0.90, 15 effect size
# ratios gamma1 and 7 correlations (210 values). In the notation of the
# twocpcont article, C_2 = lambda_w - z_alpha, where lambda_w is the
# noncentrality parameter of the weaker endpoint at which the target power is
# attained.
#
# Each published value is reproduced by two routes of the package:
#   exact        the root of the bivariate normal power equation
#                Phi2(lambda - z_alpha, gamma1 lambda - z_alpha; rho) = 1 - beta,
#                with pbivnorm::pbivnorm() and uniroot()
#   closed_form  lambda_w* of twocpcont_ss(method = "univariate"), through
#                n2_cont with delta_w = sd_w = 1 and r = 1 (c = 2)
#
# Verdicts (one row per published value and route)
#   PASS       the value rounded to three decimals equals the published value
#   EXPLAINED  closed form only: it differs by 0.001 for rho <= 0.8 while
#              the exact route reproduces the value and the closed form is
#              within 4e-4 of the exact root, or it differs by at most 0.002
#              at rho = 0.95 (the reasons are registered in
#              inst/validation/explained.csv)
#   FAIL       any other difference
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file:
#     devtools::load_all()
#     source("inst/reproduce/reproduce_Sozu2015.R")
#
# Input   inst/extdata/sozu2015_C2_table.csv (the published values)
# Output  inst/reproduce/results/reproduce_Sozu2015.csv
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

alpha    <- 0.025
gl_nodes <- 5L
z_alpha  <- qnorm(1 - alpha)

sozu_tab <- utils::read.csv(file.path("inst", "extdata",
                                      "sozu2015_C2_table.csv"))

C2_exact <- mapply(function(g, rh, pw) {
  uniroot(function(l) {
    pbivnorm::pbivnorm(l - z_alpha, g * l - z_alpha, rho = rh) - pw
  }, lower = z_alpha, upper = z_alpha + 5, tol = 1e-12)$root - z_alpha
}, sozu_tab$gamma1, sozu_tab$rho, sozu_tab$power)
C2_cf <- mapply(function(g, rh, pw) {
  ss <- twocpcont_ss(delta1 = g, delta2 = 1, sd1 = 1, sd2 = 1,
                     rho = rh, r = 1, alpha = alpha, beta = 1 - pw,
                     method = "univariate", gl_nodes = gl_nodes)
  sqrt(ss$n2_cont / 2) - z_alpha
}, sozu_tab$gamma1, sozu_tab$rho, sozu_tab$power)

dev_exact <- round(round(C2_exact, 3) - sozu_tab$C2, 3)
dev_cf    <- round(round(C2_cf, 3) - sozu_tab$C2, 3)
err_cf    <- C2_cf - C2_exact

verdict_exact <- ifelse(abs(dev_exact) < 1e-9, "PASS", "FAIL")
verdict_cf <- ifelse(
  abs(dev_cf) < 1e-9, "PASS",
  ifelse(sozu_tab$rho <= 0.8 & abs(abs(dev_cf) - 0.001) < 1e-9 &
           verdict_exact == "PASS" & abs(err_cf) < 4e-4, "EXPLAINED",
         ifelse(sozu_tab$rho > 0.8 & abs(dev_cf) <= 0.002 + 1e-9,
                "EXPLAINED", "FAIL")))
reason_cf <- ifelse(verdict_cf != "EXPLAINED", "",
                    ifelse(sozu_tab$rho <= 0.8, "closed_form_rho_le_0.8",
                           "closed_form_rho_0.95"))

key <- sprintf("Table %s, 1 - beta = %.2f, gamma1 = %.2f, rho = %.2f",
               sozu_tab$table, sozu_tab$power, sozu_tab$gamma1, sozu_tab$rho)
res <- rbind(
  data.frame(item = key, route = "exact", published = sozu_tab$C2,
             package = C2_exact, difference = dev_exact,
             verdict = verdict_exact, explained = "",
             stringsAsFactors = FALSE),
  data.frame(item = key, route = "closed_form", published = sozu_tab$C2,
             package = C2_cf, difference = dev_cf,
             verdict = verdict_cf, explained = reason_cf,
             stringsAsFactors = FALSE)
)

# Every EXPLAINED verdict must have its reason in explained.csv
explained <- utils::read.csv(file.path("inst", "validation", "explained.csv"),
                             stringsAsFactors = FALSE)
reg <- explained$item[explained$script == "reproduce_Sozu2015.R"]
unreg <- setdiff(unique(res$explained[res$verdict == "EXPLAINED"]), reg)
if (length(unreg) > 0L) {
  res$verdict[res$explained %in% unreg] <- "FAIL"
}

out_file <- file.path(results_dir, "reproduce_Sozu2015.csv")
utils::write.csv(res, out_file, row.names = FALSE)

message("reproduce_Sozu2015.R: ", nrow(sozu_tab), " published values")
for (rt in c("exact", "closed_form")) {
  tb <- table(factor(res$verdict[res$route == rt],
                     levels = c("PASS", "EXPLAINED", "FAIL")))
  message(sprintf("  %-12s PASS %d, EXPLAINED %d, FAIL %d", rt,
                  tb[["PASS"]], tb[["EXPLAINED"]], tb[["FAIL"]]))
}
message("Results written to ", out_file)
if (length(unreg) > 0L) {
  message("  EXPLAINED without a reason in explained.csv: ",
          paste(unreg, collapse = ", "))
}
if (any(res$verdict == "FAIL")) {
  stop(sprintf("reproduce_Sozu2015: %d FAIL verdict(s); see %s.",
               sum(res$verdict == "FAIL"), out_file), call. = FALSE)
}
