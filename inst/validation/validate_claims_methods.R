# =============================================================================
# validate_claims_methods.R
#
# Checks of the statements of Section 2 of the twocpcont article and of
# Section S2 of its Supporting Information that are not single numbers
# (bounds, monotonicity, inequalities). The numbers cited in the text are
# checked against the ledger by inst/scripts/number_check.R; this script
# checks the claims made about them.
#
# Claims about the computed results of the article are checked on the data
# written by inst/scripts/data_generate_main.R; claims about the functions
# of the package are checked by computing them here, against an independent
# route where one exists (adaptive quadrature for the Gauss-Legendre rule,
# uniroot() for the closed-form R_max).
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file after data_generate_main.R:
#     devtools::load_all()
#     source("inst/validation/validate_claims_methods.R")
#   For a test of the script itself, SCRIPTS_DATA_DIR (a copy of
#   inst/scripts/data) and VALIDATION_OUT_DIR can be assigned before
#   sourcing; both are removed when the script starts.
#
# Output  inst/validation/results/validate_claims_methods.csv
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

alpha     <- 0.025
power_set <- c(0.80, 0.85, 0.90)
z_alpha   <- qnorm(1 - alpha)

# =============================================================================
# Section 2.3 and S2: accuracy of the 5-point Gauss-Legendre rule
# =============================================================================
# The exact Plackett integral is computed by adaptive quadrature over the
# grid of (a, b, rho) of Section S2, which covers the arguments arising in
# the article for |rho| <= 0.8.
plackett_exact <- function(a, b, rho) {
  integrate(function(t) {
    exp(-(a ^ 2 + b ^ 2 - 2 * t * a * b) / (2 * (1 - t ^ 2))) /
      (2 * pi * sqrt(1 - t ^ 2))
  }, lower = 0, upper = rho, rel.tol = 1e-12, abs.tol = 1e-15)$value
}
gl_grid <- expand.grid(a = seq(-1, 4, by = 0.25), b = seq(-1, 8, by = 0.25),
                       rho = seq(-0.8, 0.8, by = 0.1))
gl_err <- mapply(function(a, b, rho) {
  abs(plackett_exact(a, b, rho) -
        plackett_gl_full(a, b, rho, kappa = 1, gl_nodes = 5L)$I_GL)
}, gl_grid$a, gl_grid$b, gl_grid$rho)
check("gl_error_below_1e-5",
      "|I - I_GL| < 1e-5 for |rho| <= 0.8 on the 13209-point grid of S2",
      max(gl_err) < 1e-5 && nrow(gl_grid) == 13209L,
      sprintf("grid points = %d, max error = %.2e", nrow(gl_grid),
              max(gl_err)))

# =============================================================================
# Section 2.2.1: the anchor weight w* exceeds the exact root of (S2)
# =============================================================================
fig1 <- read_scripts_data("fig1_anchor_weight")
gap <- fig1$gap
check("w_star_at_least_exact_root",
      "w* >= exact root of F(w; kappa) = 0 for 1 < kappa <= 3",
      min(gap$gap) > -1e-8 && min(gap$kappa) > 1 && max(gap$kappa) >= 3 &&
        setequal(unique(gap$power), power_set),
      sprintf("kappa in [%.3f, %.3f], min gap = %.2e, max gap = %.4f",
              min(gap$kappa), max(gap$kappa), min(gap$gap), max(gap$gap)))

# =============================================================================
# Section 2.4: maximum reduction rate and threshold effect size ratio
# =============================================================================
# The root of the rho = 0 equation (eq:rho0_root) lies in
# [z_alpha + z_beta, z_alpha + u_p]
in_bracket <- TRUE
for (b in 1 - power_set) {
  z_beta <- qnorm(1 - b)
  u_p    <- qnorm(sqrt(1 - b))
  for (k in seq(1, 3, by = 0.05)) {
    lam0 <- (z_alpha + z_beta) / sqrt(1 - r_max(kappa = k, alpha = alpha,
                                                 beta = b))
    in_bracket <- in_bracket && lam0 >= z_alpha + z_beta - 1e-10 &&
      lam0 <= z_alpha + u_p + 1e-10
  }
}
check("rmax_root_in_bracket",
      "The root of (eq:rho0_root) lies in [z_alpha + z_beta, z_alpha + u_p]",
      in_bracket, "kappa in [1, 3] by 0.05, 1 - beta in {0.80, 0.85, 0.90}")

mono_k <- all(sapply(1 - power_set, function(b) {
  v <- sapply(seq(1, 3, by = 0.01), r_max, alpha = alpha, beta = b)
  all(diff(v) <= 1e-12)
}))
check("rmax_nonincreasing_in_kappa",
      "R_max(kappa, beta) is non-increasing in kappa", mono_k,
      "kappa in [1, 3] by 0.01")

rmax_diff <- max(sapply(c(0.30, 0.20, 0.15, 0.13, 0.10), function(b) {
  k_seq <- seq(1, 3, by = 0.001)
  max(abs(sapply(k_seq, r_max, alpha = alpha, beta = b, exact = TRUE) -
            sapply(k_seq, r_max, alpha = alpha, beta = b, exact = FALSE)))
}))
check("rmax_closed_form_vs_exact_below_1e-4",
      paste0("Closed-form and exact R_max differ by less than 1e-4 for ",
             "kappa in [1, 3], 1 - beta in {0.70, 0.80, 0.85, 0.87, 0.90}"),
      rmax_diff < 1e-4, sprintf("max |difference| = %.2e", rmax_diff))

inv_err <- max(sapply(1 - power_set, function(b) {
  sapply(c(0.05, 0.10, 0.15), function(e) {
    abs(r_max(kappa = kappa_star(epsilon = e, alpha = alpha, beta = b),
              alpha = alpha, beta = b) - e)
  })
}))
check("kstar_inverse_of_rmax", "R_max(kappa_star(epsilon), beta) = epsilon",
      inv_err < 1e-10, sprintf("max |R_max - epsilon| = %.2e", inv_err))

ks_bound <- sapply(1 - power_set, function(b) {
  r1 <- r_max(kappa = 1, alpha = alpha, beta = b)
  c(kappa_star(r1, alpha = alpha, beta = b),
    kappa_star(min(0.99, r1 + 0.05), alpha = alpha, beta = b),
    kappa_star(r1 - 0.01, alpha = alpha, beta = b))
})
check("kstar_one_at_and_above_rmax1",
      paste0("kappa_star = 1 for epsilon >= R_max(1, beta) and > 1 for ",
             "epsilon < R_max(1, beta)"),
      all(abs(ks_bound[1:2, ] - 1) < 1e-8) && all(ks_bound[3, ] > 1),
      sprintf("at R_max(1): %s; just below: %s",
              paste(sprintf("%.6f", ks_bound[1, ]), collapse = ", "),
              paste(sprintf("%.4f", ks_bound[3, ]), collapse = ", ")))

arg_ok <- all(sapply(1 - power_set, function(b) {
  z_b <- qnorm(1 - b)
  sapply(seq(0.001, 0.999, by = 0.001), function(e) {
    nu  <- (z_alpha + z_b) * (1 / sqrt(1 - e) - 1)
    arg <- (1 - b) / pnorm(z_b + nu)
    arg > 0 && arg < 1
  })
}))
check("kstar_quantile_argument_in_unit_interval",
      "The argument (1 - beta) / Phi(z_beta + nu) of (eq:kappa_star) lies in (0, 1)",
      arg_ok, "epsilon in [0.001, 0.999] by 0.001")

# =============================================================================
# Results
# =============================================================================
res <- do.call(rbind, checks)
out_file <- file.path(out_dir, "validate_claims_methods.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf("validate_claims_methods.R: PASS %d, FAIL %d; written to %s",
                sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"),
                out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("validate_claims_methods: FAIL in %s.",
               paste(res$item[res$verdict == "FAIL"], collapse = ", ")),
       call. = FALSE)
}
