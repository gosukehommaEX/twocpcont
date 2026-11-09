# =============================================================================
# validate_error_decomposition.R
#
# Checks of the computations behind Figure S1 (error of the continuous n2 by
# approximation step) and Table S2 and Figure S2 (accuracy over a wide grid)
# of the Supporting Information of the twocpcont article, on the data
# written by inst/sup_info/data_generate_supplement.R:
#   - the analytic P, P_1 and P_2 at the anchor, for the Gauss-Legendre
#     approximation and for the exact bivariate normal distribution
#     function, agree with numerical derivatives;
#   - the steps (first-order, quadratic) are computed from these values, and
#     the quadratic correction of the script equals lambda_w* of the package;
#   - the root of the Gauss-Legendre power equation and the exact root solve
#     their equations;
#   - at rho = 0, where the Gauss-Legendre term vanishes, the quadrature-only
#     and Taylor-only routes coincide with the exact root and with the
#     proposed formula;
#   - the discriminant of the correction term is nonnegative in every design
#     of the two parts.
# The numbers cited in the text are checked against the ledger by
# inst/scripts/number_check.R.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file after
#   inst/sup_info/data_generate_supplement.R:
#     devtools::load_all()
#     source("inst/validation/validate_error_decomposition.R")
#   For a test of the script itself, SUPINFO_DATA_DIR (a copy of
#   inst/sup_info/data) and VALIDATION_OUT_DIR can be assigned before
#   sourcing; both are removed when the script starts.
#
# Output  inst/validation/results/validate_error_decomposition.csv
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

decomp_obj <- read_supinfo_data("error_decomposition")
dec   <- decomp_obj$data
alpha <- decomp_obj$alpha
z_alpha <- qnorm(1 - alpha)

# Power as a function of lambda_w: Gauss-Legendre approximation (5 nodes)
# and the exact bivariate normal distribution function
power_gl <- function(lam, kappa, rho) {
  u_w <- lam - z_alpha
  u_s <- kappa * lam - z_alpha
  pnorm(u_w) * pnorm(u_s) +
    plackett_gl_full(u_w, u_s, rho, kappa, gl_nodes = 5L)$I_GL
}
power_exact <- function(lam, kappa, rho) {
  pbivnorm::pbivnorm(lam - z_alpha, kappa * lam - z_alpha, rho = rho)
}

# Central differences: step H1 for the first derivative and H2 for the
# second derivative (truncation and rounding errors below 1e-10 and 1e-6)
H1 <- 1e-5
H2 <- 1e-3
num_derivs <- function(f, lam, kappa, rho) {
  f0 <- f(lam, kappa, rho)
  c(P  = f0,
    P1 = (f(lam + H1, kappa, rho) - f(lam - H1, kappa, rho)) / (2 * H1),
    P2 = (f(lam + H2, kappa, rho) - 2 * f0 + f(lam - H2, kappa, rho)) / H2 ^ 2)
}
quad_step <- function(C0, P1, P2) -2 * C0 / (P1 + sqrt(P1 ^ 2 - 2 * P2 * C0))

# =============================================================================
# Analytic P, P_1 and P_2 at the anchor against numerical derivatives
# =============================================================================
TOL_P  <- 1e-14
TOL_P1 <- 1e-8
TOL_P2 <- 1e-5
for (route in c("gl", "ex")) {
  f <- if (route == "gl") power_gl else power_exact
  nd <- t(mapply(function(l, k, rr) num_derivs(f, l, k, rr),
                 dec$lam_anchor, dec$kappa, dec$rho))
  d_P  <- max(abs(nd[, "P"] - dec[[paste0("P_", route)]]))
  d_P1 <- max(abs(nd[, "P1"] - dec[[paste0("P1_", route)]]))
  d_P2 <- max(abs(nd[, "P2"] - dec[[paste0("P2_", route)]]))
  label <- if (route == "gl") "Gauss-Legendre approximation" else
    "exact bivariate normal distribution function"
  check(sprintf("derivatives_%s", route),
        sprintf(paste0("P, P_1 and P_2 at the anchor (%s) agree with the ",
                       "function and its central differences"), label),
        d_P < TOL_P && d_P1 < TOL_P1 && d_P2 < TOL_P2,
        sprintf(paste0("designs = %d, max |P| diff = %.1e, max |P_1| diff = ",
                       "%.1e, max |P_2| diff = %.1e (tolerances %.0e, %.0e, ",
                       "%.0e)"), nrow(dec), d_P, d_P1, d_P2, TOL_P, TOL_P1,
                TOL_P2))
}

# =============================================================================
# The steps from the anchor
# =============================================================================
TOL_LAM <- 1e-12
tp <- dec$power
d_first <- max(abs(dec$lam_anchor - (dec$P_gl - tp) / dec$P1_gl -
                     dec$lam_first))
check("first_order_step",
      paste0("The first-order correction is one Newton step on the ",
             "Gauss-Legendre power function from the anchor"),
      d_first < TOL_LAM, sprintf("max |difference| = %.1e", d_first))

d_prop <- max(abs(dec$lam_anchor +
                    quad_step(dec$P_gl - tp, dec$P1_gl, dec$P2_gl) -
                    dec$lam_proposed))
d_taylor <- max(abs(dec$lam_anchor +
                      quad_step(dec$P_ex - tp, dec$P1_ex, dec$P2_ex) -
                      dec$lam_exact_quad))
check("quadratic_steps",
      paste0("The proposed and the Taylor-only routes are the quadratic ",
             "correction with the Gauss-Legendre and the exact P, P_1, P_2"),
      d_prop < TOL_LAM && d_taylor < TOL_LAM,
      sprintf("max |difference| = %.1e (proposed), %.1e (Taylor only)",
              d_prop, d_taylor))

d_pkg <- max(abs(dec$lam_proposed - dec$lam_package))
check("proposed_equals_package",
      paste0("The quadratic correction of the script equals lambda_w* of ",
             "twocpcont_ss(method = \"univariate\")"),
      d_pkg < TOL_LAM, sprintf("max |difference| = %.1e", d_pkg))

# =============================================================================
# The roots solve their equations
# =============================================================================
TOL_RES <- 1e-12
res_gl <- max(abs(mapply(power_gl, dec$lam_gl_root, dec$kappa, dec$rho) - tp))
check("gl_root_solves",
      paste0("The quadrature-only route is a root of the Gauss-Legendre ",
             "power equation"),
      res_gl < TOL_RES, sprintf("max |P(lambda) - (1 - beta)| = %.1e", res_gl))

res_ex <- max(abs(mapply(power_exact, dec$lam_exact, dec$kappa, dec$rho) - tp))
grid_obj <- read_supinfo_data("accuracy_grid")
acc <- grid_obj$data
res_grid <- max(abs(mapply(function(l, k, rr, a) {
  pbivnorm::pbivnorm(l - qnorm(1 - a), k * l - qnorm(1 - a), rho = rr)
}, acc$lam_exact, acc$kappa, acc$rho, acc$alpha) - acc$power))
check("exact_root_solves",
      paste0("The exact lambda_w of the error decomposition and of the ",
             "accuracy grid solves the power equation with the bivariate ",
             "normal distribution function"),
      res_ex < TOL_RES && res_grid < TOL_RES,
      sprintf(paste0("max |Phi2 - (1 - beta)| = %.1e (decomposition, %d ",
                     "designs), %.1e (grid, %d designs)"), res_ex, nrow(dec),
              res_grid, nrow(acc)))

# =============================================================================
# rho = 0: the Gauss-Legendre term vanishes
# =============================================================================
d0 <- dec[dec$rho == 0, ]
d0_root <- max(abs(d0$lam_gl_root - d0$lam_exact))
d0_quad <- max(abs(d0$lam_proposed - d0$lam_exact_quad))
check("rho0_routes_coincide",
      paste0("At rho = 0 the quadrature-only route equals the exact root and ",
             "the proposed route equals the Taylor-only route"),
      nrow(d0) > 0L && d0_root < 1e-10 && d0_quad < TOL_LAM,
      sprintf("designs = %d, max |difference| = %.1e (roots), %.1e (steps)",
              nrow(d0), d0_root, d0_quad))

# =============================================================================
# Relative errors and the accuracy grid
# =============================================================================
d_rel <- max(vapply(c("anchor", "first", "proposed", "gl_root", "exact_quad"),
                    function(v) {
  max(abs((dec[[paste0("lam_", v)]] / dec$lam_exact) ^ 2 - 1 -
            dec[[paste0("rel_n2_", v)]]))
}, numeric(1L)))
d_rel_grid <- max(abs((acc$lam_cf / acc$lam_exact) ^ 2 - 1 - acc$rel_n2))
d_pow_grid <- max(abs(mapply(function(l, k, rr, a) {
  pbivnorm::pbivnorm(l - qnorm(1 - a), k * l - qnorm(1 - a), rho = rr)
}, acc$lam_cf, acc$kappa, acc$rho, acc$alpha) - acc$power - acc$power_diff))
check("relative_errors",
      paste0("The relative errors of the continuous n2 and the power ",
             "differences are those of the stored noncentrality parameters"),
      d_rel < 1e-14 && d_rel_grid < 1e-14 && d_pow_grid < 1e-14,
      sprintf(paste0("max |difference| = %.1e (decomposition), %.1e ",
                     "(grid), %.1e (power, grid)"), d_rel, d_rel_grid,
              d_pow_grid))

check("discriminant_nonnegative_si",
      paste0("The discriminant of the correction term is nonnegative in ",
             "every design of the error decomposition and of the accuracy ",
             "grid"),
      all(dec$disc_gl >= 0) && all(acc$disc >= 0),
      sprintf("designs = %d, min discriminant = %.4f", nrow(dec) + nrow(acc),
              min(c(dec$disc_gl, acc$disc))))

# =============================================================================
# Results
# =============================================================================
res <- do.call(rbind, checks)
out_file <- file.path(out_dir, "validate_error_decomposition.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf(
  "validate_error_decomposition.R: PASS %d, FAIL %d; written to %s",
  sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"), out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("validate_error_decomposition: FAIL in %s.",
               paste(res$item[res$verdict == "FAIL"], collapse = ", ")),
       call. = FALSE)
}
