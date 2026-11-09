# =============================================================================
# validate_claims_example.R
#
# Checks of the statements of Section 4 (real trial example), Section 5
# (discussion) and Section S5 of the Supporting Information (software) of
# the twocpcont article that are not single numbers: monotonicity and
# agreement of the two methods in Table 4, the position of the effect size
# ratio relative to the thresholds, the protocol's table of marginal powers,
# and the power checks described in the discussion and in Section S5. The
# numbers cited in the text are checked against the ledger by
# inst/scripts/number_check.R.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file after data_generate_main.R:
#     devtools::load_all()
#     source("inst/validation/validate_claims_example.R")
#   For a test of the script itself, SCRIPTS_DATA_DIR (a copy of
#   inst/scripts/data) and VALIDATION_OUT_DIR can be assigned before
#   sourcing; both are removed when the script starts.
#
# Output  inst/validation/results/validate_claims_example.csv
#         (item, claim, verdict PASS or FAIL, detail)
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

scripts_data_dir <- file.path("inst", "scripts", "data")
out_dir <- file.path("inst", "validation", "results")
if (exists("SCRIPTS_DATA_DIR", envir = globalenv())) {
  scripts_data_dir <- get("SCRIPTS_DATA_DIR", envir = globalenv())
  rm("SCRIPTS_DATA_DIR", envir = globalenv())
}
if (exists("VALIDATION_OUT_DIR", envir = globalenv())) {
  out_dir <- get("VALIDATION_OUT_DIR", envir = globalenv())
  rm("VALIDATION_OUT_DIR", envir = globalenv())
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

read_scripts_data <- function(name) {
  path <- file.path(scripts_data_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run inst/scripts/data_generate_main.R ",
                 path), "first.", call. = FALSE)
  }
  readRDS(path)
}

checks <- list()
check <- function(item, claim, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    item = item, claim = claim,
    verdict = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = detail, stringsAsFactors = FALSE)
}

# =============================================================================
# Section 4: EXPEDITION 1 (Table 4)
# =============================================================================
tab4 <- read_scripts_data("tab4_real_example")
t4   <- tab4$data[order(tab4$data$rho), ]

check("trial_N_nonincreasing_in_rho",
      "N_prop decreases monotonically as rho increases from 0 to 0.8",
      all(diff(t4$N_prop) <= 0),
      sprintf("N_prop = %s", paste(t4$N_prop, collapse = ", ")))
check("trial_methods_agree",
      "The two methods agree at the four correlations",
      nrow(t4) == 4L && all(t4$N_prop == t4$N_conv),
      sprintf("N_conv = %s; N_prop = %s", paste(t4$N_conv, collapse = ", "),
              paste(t4$N_prop, collapse = ", ")))
check("trial_N_rho0_equals_protocol",
      "At rho = 0 the formula gives the 2 x 500 patients of the protocol",
      t4$N_prop[t4$rho == 0] == 2 * 500,
      sprintf("N_prop(rho = 0) = %d", t4$N_prop[t4$rho == 0]))

check("trial_rmax_below_5pct",
      "R_max at kappa = 1.292 and 1 - beta = 0.87 is below 5%",
      tab4$R_max_pct < 5, sprintf("R_max = %.4f%%", tab4$R_max_pct))

red_int  <- t4$reduction_pct[t4$rho == 0.8]
red_cont <- tab4$lambda$reduction_cont_pct[tab4$lambda$rho == 0.8]
check("trial_continuous_reduction_below_rmax_and_integer",
      paste0("The continuous reduction at rho = 0.8 is below R_max and below ",
             "the integer reduction"),
      red_cont < tab4$R_max_pct && red_cont < red_int,
      sprintf("continuous = %.4f%%, integer = %.2f%%, R_max = %.4f%%",
              red_cont, red_int, tab4$R_max_pct))

check("trial_kstar05_lt_kappa_lt_kstar01",
      "kappa_star(0.05, 0.13) < kappa = 1.292 < kappa_star(0.01, 0.13)",
      tab4$kappa_star_eps05 < tab4$kappa && tab4$kappa < tab4$kappa_star_eps01,
      sprintf("%.4f < %.4f < %.4f", tab4$kappa_star_eps05, tab4$kappa,
              tab4$kappa_star_eps01))
check("trial_kappa_round_below_kstar05",
      paste0("With the rounded effect size 0.25, kappa = 1.25 is below ",
             "kappa_star(0.05, 0.13) and R_max is above 5%"),
      tab4$kappa_round < tab4$kappa_star_eps05 && tab4$R_max_round_pct > 5,
      sprintf("kappa = %.4f, kappa_star(0.05) = %.4f, R_max = %.4f%%",
              tab4$kappa_round, tab4$kappa_star_eps05, tab4$R_max_round_pct))

# The protocol's table of marginal powers (Table LZAM.4: 81/85/89% for
# effect size 0.20 and 94/96/98% for 0.25 at 400/450/500 patients per group,
# two-sided 0.05) agrees with 0.25 and not with 3.1 / 12
prot_pow <- function(es) {
  round(100 * pnorm(es * sqrt(c(400, 450, 500) / 2) - qnorm(0.975)))
}
check("protocol_power_table_matches_0.25",
      paste0("The protocol's table of marginal powers agrees with the ",
             "effect size 0.25 and not with 3.1 / 12"),
      all(prot_pow(0.20) == c(81, 85, 89)) &&
        all(prot_pow(0.25) == c(94, 96, 98)) &&
        any(prot_pow(3.1 / 12) != c(94, 96, 98)),
      sprintf("0.20: %s; 0.25: %s; 3.1/12: %s",
              paste(prot_pow(0.20), collapse = "/"),
              paste(prot_pow(0.25), collapse = "/"),
              paste(prot_pow(3.1 / 12), collapse = "/")))

# =============================================================================
# Section 5: discussion
# =============================================================================
incr_t <- qnorm(1 - 0.025) ^ 2 / 2
check("schouten_increment_about_two",
      "z_alpha^2 / 2 at alpha = 0.025 is about two patients",
      round(incr_t) == 2, sprintf("z_alpha^2 / 2 = %.4f", incr_t))

# One power check with n2 increased by one gives the conventional sample
# size in the shortfall of Section 3.2
grid <- read_scripts_data("tab2_sample_size")$data
mis  <- grid[grid$N_prop != grid$N_conv, , drop = FALSE]
fix_ok <- FALSE
fix_detail <- sprintf("%d mismatch(es)", nrow(mis))
if (nrow(mis) == 1L) {
  n2_fix <- mis$n2_prop + 1
  n1_fix <- ceiling(mis$r * n2_fix)
  pw_fix <- twocpcont_power(n1 = n1_fix, n2 = n2_fix,
                            delta1 = mis$delta1, delta2 = mis$delta2,
                            sd1 = mis$sd1, sd2 = mis$sd2, rho = mis$rho,
                            alpha = 0.025)$powerCoprimary
  fix_ok <- pw_fix >= mis$power && n1_fix + n2_fix == mis$N_conv
  fix_detail <- sprintf("n2 = %d, N = %d, power = %.7f", as.integer(n2_fix),
                        as.integer(n1_fix + n2_fix), pw_fix)
}
check("shortfall_fixed_by_one_power_check",
      paste0("Increasing n2 by one when the target is missed gives the ",
             "conventional sample size in the shortfall of Section 3.2"),
      fix_ok, fix_detail)

# =============================================================================
# Section S5: software example (EXPEDITION 1 with rho = 0.5)
# =============================================================================
ss_app <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                       rho = 0.5, r = 1, alpha = 0.025, beta = 0.13,
                       method = "univariate")
pw_app <- twocpcont_power(n1 = ss_app$n1, n2 = ss_app$n2,
                          delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                          rho = 0.5, alpha = 0.025)$powerCoprimary
check("software_example_matches_table4",
      paste0("In Section S5 the total sample size equals that of Table 4 at ",
             "rho = 0.5 and the co-primary power is at least 0.87"),
      ss_app$N == t4$N_prop[t4$rho == 0.5] && pw_app >= 0.87,
      sprintf("N = %d (Table 4: %d), power = %.6f", ss_app$N,
              t4$N_prop[t4$rho == 0.5], pw_app))

# =============================================================================
# Results
# =============================================================================
res <- do.call(rbind, checks)
out_file <- file.path(out_dir, "validate_claims_example.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf("validate_claims_example.R: PASS %d, FAIL %d; written to %s",
                sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"),
                out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("validate_claims_example: FAIL in %s.",
               paste(res$item[res$verdict == "FAIL"], collapse = ", ")),
       call. = FALSE)
}
