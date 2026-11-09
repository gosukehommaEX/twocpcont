# =============================================================================
# validate_claims_supplement.R
#
# Checks of the statements of Section 3 and of Sections S5 to S8 of the
# Supporting Information of the twocpcont article that are not single
# numbers, on the data written by inst/sup_info/data_generate_supplement.R:
#   - Table S1: the closed form evaluates no bivariate normal distribution
#     function, it is not faster than Brent's method, and the computing
#     times of the Newton-Raphson algorithm, Brent's method and the closed
#     form are of the same order;
#   - Figure S1: the anchor is exact at kappa = 1 and rho = 0, the error of
#     the quadrature alone is much smaller than that of the Taylor
#     approximation alone for rho <= 0.8, the error of the proposed formula
#     is close to the latter, and the two errors are of the same order at
#     rho = 0.95;
#   - Table S2 and Figure S2: the error is largest for kappa close to one and
#     small for kappa >= 2, it changes little with alpha, it is of the same
#     order for negative rho as at rho = 0, and the power at the continuous
#     n2 differs from the target by less than 1e-3 for rho <= 0.9;
#   - Table S3: the closed form never gives a smaller total than the
#     sequential search, the differences occur mostly with r < 1, and the
#     power at the total of the closed form is at least the target;
#   - Table S4 and Figure S3: a planning correlation above the true one gives
#     a power below the target, planning at rho = 0 attains the target at
#     every nonnegative true correlation but not at a negative one, and the
#     effect size of ADCS-ADL at kappa_star(0.05, 0.13) lies between the
#     value of the protocol and that of Table 4.
# The numbers cited in the text are checked against the ledger by
# inst/scripts/number_check.R.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file after
#   inst/sup_info/data_generate_supplement.R:
#     devtools::load_all()
#     source("inst/validation/validate_claims_supplement.R")
#   For a test of the script itself, SUPINFO_DATA_DIR (a copy of
#   inst/sup_info/data) and VALIDATION_OUT_DIR can be assigned before
#   sourcing; both are removed when the script starts.
#
# Output  inst/validation/results/validate_claims_supplement.csv
#         (item, claim, verdict PASS or FAIL, detail)
# The script stops with an error if any verdict is FAIL.
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "twocpcont")) {
  stop("Run this script from the root of the twocpcont package ",
       "(open twocpcont.Rproj first).", call. = FALSE)
}

supinfo_data_dir <- file.path("inst", "sup_info", "data")
out_dir <- file.path("inst", "validation", "results")
if (exists("SUPINFO_DATA_DIR", envir = globalenv())) {
  supinfo_data_dir <- get("SUPINFO_DATA_DIR", envir = globalenv())
  rm("SUPINFO_DATA_DIR", envir = globalenv())
}
if (exists("VALIDATION_OUT_DIR", envir = globalenv())) {
  out_dir <- get("VALIDATION_OUT_DIR", envir = globalenv())
  rm("VALIDATION_OUT_DIR", envir = globalenv())
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

read_supinfo_data <- function(name) {
  path <- file.path(supinfo_data_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run ", path),
         "inst/sup_info/data_generate_supplement.R first.", call. = FALSE)
  }
  readRDS(path)
}

checks <- list()
check <- function(item, claim, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    item = item, claim = claim, verdict = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = detail, stringsAsFactors = FALSE)
}

# =============================================================================
# Table S1: computational cost (Section 3.4 and Section S5)
# =============================================================================
cost_obj <- read_supinfo_data("cost_comparison")
cd <- cost_obj$data
cf <- cd[cd$method == "univariate", ]
check("closed_form_no_phi2_evaluation",
      "The closed form evaluates Phi2 at no point",
      nrow(cf) == 96L && all(cf$n_phi2 == 0L),
      sprintf("cells = %d, max count = %d", nrow(cf), max(cf$n_phi2)))

ms <- stats::setNames(cost_obj$timing$ms_per_design, cost_obj$timing$method)
check("closed_form_not_faster_than_brent",
      "The closed form is not faster than Brent's method (uniroot)",
      ms[["univariate"]] >= ms[["uniroot"]],
      sprintf("ms per design: closed form %.3f, uniroot %.3f",
              ms[["univariate"]], ms[["uniroot"]]))

ms3 <- ms[c("newton", "uniroot", "univariate")]
check("iterative_and_closed_form_same_order",
      paste0("The computing times of the Newton-Raphson algorithm, Brent's ",
             "method and the closed form are of the same order (ratio < 10)"),
      max(ms3) / min(ms3) < 10,
      sprintf("ms per design: %s; ratio = %.2f",
              paste(sprintf("%s %.3f", names(ms3), ms3), collapse = ", "),
              max(ms3) / min(ms3)))

# =============================================================================
# Figure S1: error by approximation step (Section 3.3 and Section S6)
# =============================================================================
dec <- read_supinfo_data("error_decomposition")$data
abs_rel <- function(d, v) abs(d[[paste0("rel_n2_", v)]])

at_k1 <- abs(dec$kappa - 1) < 1e-9 & abs(dec$rho) < 1e-9
check("anchor_exact_at_kappa1_rho0",
      "At kappa = 1 and rho = 0 the anchor is the exact root",
      sum(at_k1) > 0L && max(abs_rel(dec[at_k1, ], "anchor")) < 1e-12,
      sprintf("designs = %d, max |relative error| = %.1e", sum(at_k1),
              max(abs_rel(dec[at_k1, ], "anchor"))))

lo <- dec[dec$rho <= 0.8 + 1e-9, ]
gl_lo <- max(abs_rel(lo, "gl_root"))
tq_lo <- max(abs_rel(lo, "exact_quad"))
check("quadrature_error_small_rho_le_0.8",
      paste0("For rho <= 0.8 the error of the quadrature alone is much ",
             "smaller (less than a tenth) than that of the Taylor ",
             "approximation alone"),
      gl_lo < 0.1 * tq_lo,
      sprintf("max |relative error|: quadrature %.1e, Taylor %.1e",
              gl_lo, tq_lo))

d_pt <- max(abs(lo$rel_n2_proposed - lo$rel_n2_exact_quad))
check("proposed_close_to_taylor_rho_le_0.8",
      paste0("For rho <= 0.8 the error of the proposed formula is close to ",
             "that of the Taylor approximation alone (difference less than ",
             "a tenth of the largest Taylor error)"),
      d_pt < 0.1 * tq_lo,
      sprintf("max |difference| = %.1e, max Taylor error = %.1e", d_pt,
              tq_lo))

hi <- dec[abs(dec$rho - 0.95) < 1e-9, ]
r95 <- max(abs_rel(hi, "gl_root")) / max(abs_rel(hi, "exact_quad"))
check("quadrature_and_taylor_same_order_rho_0.95",
      paste0("At rho = 0.95 the errors of the quadrature alone and of the ",
             "Taylor approximation alone are of the same order (ratio ",
             "between 0.1 and 10)"),
      nrow(hi) > 0L && r95 > 0.1 && r95 < 10,
      sprintf("ratio of the largest errors = %.2f", r95))

# =============================================================================
# Table S2 and Figure S2: accuracy over a wide grid (Section 3.3 and S6)
# =============================================================================
acc <- read_supinfo_data("accuracy_grid")$data
acc$abs_rel <- abs(acc$rel_n2)

big_k <- acc$kappa >= 2 - 1e-9
check("error_small_for_kappa_ge_2",
      "For kappa >= 2 the relative error is below 1e-4 at every correlation",
      sum(big_k) > 0L && max(acc$abs_rel[big_k]) < 1e-4,
      sprintf("designs = %d, max |relative error| = %.1e", sum(big_k),
              max(acc$abs_rel[big_k])))

rho_set <- sort(unique(acc$rho))
k_at_max <- sapply(rho_set, function(rr) {
  d <- acc[abs(acc$rho - rr) < 1e-9, ]
  d$kappa[which.max(d$abs_rel)]
})
check("error_largest_near_kappa_1",
      paste0("For every correlation the largest relative error occurs at ",
             "kappa close to one (at most 1.5)"),
      all(k_at_max <= 1.5 + 1e-9),
      sprintf("kappa at the largest error: %s",
              paste(sprintf("%g: %.2f", rho_set, k_at_max), collapse = ", ")))

alpha_ratio <- sapply(rho_set, function(rr) {
  d <- acc[abs(acc$rho - rr) < 1e-9, ]
  m <- tapply(d$abs_rel, d$alpha, max)
  max(m) / min(m)
})
check("error_changes_little_with_alpha",
      paste0("For every correlation the largest relative errors at the four ",
             "significance levels differ by less than a factor of 2"),
      all(alpha_ratio < 2),
      sprintf("largest ratio = %.2f", max(alpha_ratio)))

neg_ratio <- max(acc$abs_rel[acc$rho < 0]) /
  max(acc$abs_rel[abs(acc$rho) < 1e-9])
check("negative_rho_same_order_as_zero",
      paste0("For negative rho the largest relative error is of the same ",
             "order as at rho = 0 (ratio between 0.1 and 10)"),
      any(acc$rho < 0) && neg_ratio > 0.1 && neg_ratio < 10,
      sprintf("ratio = %.2f", neg_ratio))

pow_lo <- max(abs(acc$power_diff[acc$rho <= 0.9 + 1e-9]))
check("power_difference_below_1e-3",
      paste0("For rho <= 0.9 the power at the continuous closed-form n2 ",
             "differs from the target by less than 1e-3"),
      pow_lo < 1e-3, sprintf("max |difference| = %.1e", pow_lo))

# =============================================================================
# Table S3: allocation ratios (Section 3.2 and Section S7)
# =============================================================================
al <- read_supinfo_data("allocation")$data
d_N <- al$N_prop - al$N_conv
check("closed_form_never_smaller",
      "The closed form never gives a smaller total than the sequential search",
      nrow(al) == 120L && all(d_N >= 0),
      sprintf("combinations = %d, smallest N_prop - N_conv = %d", nrow(al),
              as.integer(min(d_N))))

ne <- d_N != 0
check("differences_mostly_r_below_1",
      "The combinations with N_prop != N_conv mostly have r < 1",
      sum(ne) > 0L && sum(ne & al$r < 1) > sum(ne) / 2,
      sprintf("differences = %d, of which r < 1: %d", sum(ne),
              sum(ne & al$r < 1)))

check("closed_form_power_at_least_target",
      "The power at the total of the closed form is at least the target",
      all(al$power_prop >= al$power),
      sprintf("smallest excess = %.5f", min(al$power_prop - al$power)))

# =============================================================================
# Table S4 and Figure S3: planning correlation and effect size (Sections 3.5,
# 4 and S8)
# =============================================================================
sens <- read_supinfo_data("expedition_sensitivity")
pl <- sens$planning
over <- pl$rho_plan > pl$rho_true
check("planning_rho_above_true_loses_power",
      paste0("When the planning correlation exceeds the true one, the power ",
             "falls below the target"),
      sum(over) > 0L && all(pl$power_true[over] < pl$power[over]),
      sprintf("rows = %d, largest power - target = %.4f", sum(over),
              max(pl$power_true[over] - pl$power[over])))

p0_nonneg <- pl$rho_plan == 0 & pl$rho_true >= 0
p0_neg    <- pl$rho_plan == 0 & pl$rho_true < 0
check("planning_at_zero_conservative_only_for_nonnegative_rho",
      paste0("Planning at rho = 0 attains the target at every nonnegative ",
             "true correlation but not at a negative one, in both designs"),
      sum(p0_nonneg) > 0L && sum(p0_neg) > 0L &&
        all(pl$power_true[p0_nonneg] >= pl$power[p0_nonneg]) &&
        all(pl$power_true[p0_neg] < pl$power[p0_neg]),
      sprintf(paste0("smallest excess (rho >= 0) = %.4f, largest excess ",
                     "(rho < 0) = %.4f"),
              min(pl$power_true[p0_nonneg] - pl$power[p0_nonneg]),
              max(pl$power_true[p0_neg] - pl$power[p0_neg])))

thr <- sens$thresholds
es2_k05 <- sens$es1 * thr$kappa_star[abs(thr$epsilon - 0.05) < 1e-9]
check("kstar_effect_size_between_protocol_and_table4",
      paste0("The effect size of ADCS-ADL at kappa_star(0.05, 0.13) lies ",
             "between the value of the protocol and that of Table 4"),
      length(es2_k05) == 1L && es2_k05 > sens$es2_protocol &&
        es2_k05 < sens$es2_design,
      sprintf("protocol %.3f < %.4f < Table 4 %.4f", sens$es2_protocol,
              es2_k05, sens$es2_design))

# =============================================================================
# Results
# =============================================================================
res <- do.call(rbind, checks)
out_file <- file.path(out_dir, "validate_claims_supplement.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf(
  "validate_claims_supplement.R: PASS %d, FAIL %d; written to %s",
  sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"), out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("validate_claims_supplement: FAIL in %s.",
               paste(res$item[res$verdict == "FAIL"], collapse = ", ")),
       call. = FALSE)
}
