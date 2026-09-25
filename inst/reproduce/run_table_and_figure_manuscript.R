# Heavy computation script for manuscript tables and figures.
#
# Performs all numerical computations needed to produce the manuscript
# tables and figures under the known-variance assumption and saves the
# results, together with timing and session information, to .rds files
# in the subdirectory "table_and_figure_manuscript" (created
# automatically if missing).
#
# Output files (under table_and_figure_manuscript/):
#   gl_nodes.rds                : Table 1 raw data
#   anchor_weight.rds           : Figure 1 (anchor weight w*) data
#   table2_sample_size_comparison.rds   : Table 2 raw data
#   reduction_curve.rds         : Figure 2 (R_max curve) data
#   kappa_star_table.rds        : Table 3 (kappa_star) raw data
#   table4_real_example.rds     : Table 4 (real trial example) raw data
#
# Companion script create_table_and_figure_manuscript.R reads these
# .rds files and writes LaTeX tables and PDF figures.
#
# Before running this script, set the working directory to a writable
# folder that contains this script.  If the R source files of
# twocpcont are also in that folder, they are used; otherwise the
# installed twocpcont package is used.

# ------------------------------------------------------------
# Load the functions.  When the R source files of twocpcont are in
# the working directory (layout of the flat code supplement), they
# are sourced; otherwise the installed twocpcont package is used.
# ------------------------------------------------------------
src_files <- c("GL_nodes_and_weights.R",
               "plackett_gl_full.R",
               "twocpcont_power.R",
               "twocpcont_ss.R",
               "print_twocpcont_ss.R",
               "print_twocpcont_power.R",
               "r_max.R",
               "kappa_star.R")
if (all(file.exists(src_files))) {
  for (f in src_files) source(f)
} else {
  library(twocpcont)
}

library(pbivnorm)

# ------------------------------------------------------------
# Output directory for all cached .rds files
# ------------------------------------------------------------
output_dir <- "table_and_figure_manuscript"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
out_path <- function(name) file.path(output_dir, name)

# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------
alpha     <- 0.025
power_set <- c(0.80, 0.85, 0.90)  # Target overall powers 1 - beta
gl_nodes  <- 5L   # Number of Gauss-Legendre quadrature points

# Scenario design for Table 2: 8 patterns (2 delta x 2 sigma x 2 r)
# x 4 rho x 3 target powers = 96 cells
# - delta:  (0.3, 0.3) for delta1 = delta2; (0.7, 0.5) for delta1 > delta2
# - sigma:  (1, 1)     for sigma1 = sigma2; (1, 1.5)   for sigma1 < sigma2
# - r:      1 or 2.  Because the standard deviations are common to the
#           two groups, equal allocation minimizes N for any (sigma1,
#           sigma2); r = 2 represents designs that allocate more patients
#           to one group for reasons other than efficiency.
scenarios <- list(
  list(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1.0, r = 1),
  list(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1.0, r = 2),
  list(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1.5, r = 1),
  list(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1.5, r = 2),
  list(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.0, r = 1),
  list(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.0, r = 2),
  list(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5, r = 1),
  list(delta1 = 0.7, delta2 = 0.5, sd1 = 1, sd2 = 1.5, r = 2)
)
rho_vals <- c(0.0, 0.3, 0.5, 0.8)

# ------------------------------------------------------------
# Table 1: Gauss-Legendre nodes and weights
# ------------------------------------------------------------
cat(sprintf("--- Computing Table 1 (%d-point Gauss-Legendre) ---\n", gl_nodes))
t1_start <- Sys.time()
gl <- GL_nodes_and_weights(gl_nodes = gl_nodes)
t1_elapsed <- as.numeric(difftime(Sys.time(), t1_start, units = "secs"))

gl_obj <- list(
  x          = gl$x,
  w          = gl$w,
  gl_nodes   = gl_nodes,
  elapsed    = t1_elapsed,
  computed   = Sys.time(),
  R_version  = R.version.string
)
saveRDS(gl_obj, out_path("gl_nodes.rds"))
cat(sprintf("  Saved: gl_nodes.rds (elapsed %.3f s)\n", t1_elapsed))

# ------------------------------------------------------------
# Figure 1 data: anchor weight w* as a function of kappa
# Implements eq. (eq:w_star) directly:
#   w_star = 1/2 + (kappa - 1) (u_p + z_alpha) phi(u_p) /
#                  [(1 + kappa) |2 log p| p]
# ------------------------------------------------------------
cat("\n--- Computing anchor weight curve (Figure 1) ---\n")
t3_start <- Sys.time()

w_star_fun <- function(kappa, beta, alpha) {
  z_alpha <- qnorm(1 - alpha)
  p       <- sqrt(1 - beta)
  u_p     <- qnorm(p)
  phi_u_p <- dnorm(u_p)
  L_abs   <- abs(2 * log(p))
  w <- 0.5 +
    (kappa - 1) * (u_p + z_alpha) * phi_u_p /
    ((1 + kappa) * L_abs * p)
  pmin(1, pmax(0, w))
}

kappa_seq_w <- seq(1, 4, length.out = 301)
power_set_w <- power_set

df_w <- do.call(rbind, lapply(power_set_w, function(pw) {
  data.frame(
    kappa  = kappa_seq_w,
    power  = pw,
    w_star = w_star_fun(kappa_seq_w, beta = 1 - pw, alpha = alpha)
  )
}))
t3_elapsed <- as.numeric(difftime(Sys.time(), t3_start, units = "secs"))

anchor_obj <- list(
  data       = df_w,
  alpha      = alpha,
  power_set  = power_set_w,
  kappa_seq  = kappa_seq_w,
  elapsed    = t3_elapsed,
  computed   = Sys.time(),
  R_version  = R.version.string
)
saveRDS(anchor_obj, out_path("anchor_weight.rds"))
cat(sprintf("  Saved: anchor_weight.rds (elapsed %.3f s)\n", t3_elapsed))

# ------------------------------------------------------------
# Helper: one cell of the Table 2 grid
# ------------------------------------------------------------
compute_cell <- function(sc, rho, beta) {

  ss_conv <- twocpcont_ss(
    delta1 = sc$delta1, delta2 = sc$delta2,
    sd1 = sc$sd1, sd2 = sc$sd2, rho = rho, r = sc$r,
    alpha = alpha, beta = beta,
    method = "bivariate"
  )

  ss_prop <- twocpcont_ss(
    delta1 = sc$delta1, delta2 = sc$delta2,
    sd1 = sc$sd1, sd2 = sc$sd2, rho = rho, r = sc$r,
    alpha = alpha, beta = beta,
    method = "univariate", gl_nodes = gl_nodes
  )

  pwr_prop <- twocpcont_power(
    n1 = ss_prop$n1, n2 = ss_prop$n2,
    delta1 = sc$delta1, delta2 = sc$delta2,
    sd1 = sc$sd1, sd2 = sc$sd2, rho = rho,
    alpha = alpha, method = "bivariate"
  )

  list(
    delta1 = sc$delta1,
    delta2 = sc$delta2,
    sd1    = sc$sd1,
    sd2    = sc$sd2,
    r      = sc$r,
    rho    = rho,
    # Rounded so that the value matches power_set exactly
    power  = round(1 - beta, 10),
    N_conv = ss_conv$N,
    N_prop = ss_prop$N,
    diff   = ss_prop$N - ss_conv$N,
    pow    = pwr_prop$powerCoprimary
  )
}

# ------------------------------------------------------------
# Table 2: Sample size comparison
# ------------------------------------------------------------
cat(sprintf("\n--- Computing Table 2 (%d cells) ---\n",
            length(scenarios) * length(rho_vals) * length(power_set)))
t2_start <- Sys.time()
rows_known <- list()
for (pw in power_set) {
  for (sc in scenarios) {
    for (rho in rho_vals) {
      rows_known[[length(rows_known) + 1L]] <- compute_cell(sc, rho, 1 - pw)
    }
  }
}
t2_elapsed <- as.numeric(difftime(Sys.time(), t2_start, units = "secs"))

df_known <- do.call(rbind, lapply(rows_known, as.data.frame))

table2_obj <- list(
  data       = df_known,
  alpha      = alpha,
  power_set  = power_set,
  gl_nodes   = gl_nodes,
  elapsed    = t2_elapsed,
  computed   = Sys.time(),
  R_version  = R.version.string
)
saveRDS(table2_obj, out_path("table2_sample_size_comparison.rds"))
cat(sprintf("  Saved: table2_sample_size_comparison.rds (elapsed %.2f s)\n",
            t2_elapsed))

# ------------------------------------------------------------
# Figure 2 data: R_max(kappa) curve at three power levels
# ------------------------------------------------------------
cat("\n--- Computing R_max(kappa) curve data ---\n")
t4_start <- Sys.time()

kappa_seq <- seq(1.0, 2.0, length.out = 201)
beta_set  <- 1 - power_set

rows_rmax <- list()
for (k in seq_along(beta_set)) {
  for (kk in kappa_seq) {
    rows_rmax[[length(rows_rmax) + 1L]] <- list(
      kappa  = kk,
      beta   = beta_set[k],
      power  = power_set[k],
      R_max  = r_max(kappa = kk, alpha = alpha, beta = beta_set[k])
    )
  }
}
df_rmax <- do.call(rbind, lapply(rows_rmax, as.data.frame))
t4_elapsed <- as.numeric(difftime(Sys.time(), t4_start, units = "secs"))

reduction_obj <- list(
  data       = df_rmax,
  alpha      = alpha,
  power_set  = power_set,
  kappa_seq  = kappa_seq,
  elapsed    = t4_elapsed,
  computed   = Sys.time(),
  R_version  = R.version.string
)
saveRDS(reduction_obj, out_path("reduction_curve.rds"))
cat(sprintf("  Saved: reduction_curve.rds (elapsed %.3f s)\n", t4_elapsed))

# ------------------------------------------------------------
# Table 3 data: kappa_star(epsilon, beta) for selected tolerances
# ------------------------------------------------------------
cat("\n--- Computing kappa_star table ---\n")
t5_start <- Sys.time()

epsilon_set <- c(0.05, 0.10, 0.15)
rows_kstar <- list()
for (p in power_set) {
  for (eps in epsilon_set) {
    rows_kstar[[length(rows_kstar) + 1L]] <- list(
      power      = p,
      beta       = 1 - p,
      epsilon    = eps,
      kappa_star = kappa_star(epsilon = eps, alpha = alpha,
                              beta = 1 - p)
    )
  }
}
df_kstar <- do.call(rbind, lapply(rows_kstar, as.data.frame))
t5_elapsed <- as.numeric(difftime(Sys.time(), t5_start, units = "secs"))

kappa_star_obj <- list(
  data        = df_kstar,
  alpha       = alpha,
  power_set   = power_set,
  epsilon_set = epsilon_set,
  elapsed     = t5_elapsed,
  computed    = Sys.time(),
  R_version   = R.version.string
)
saveRDS(kappa_star_obj, out_path("kappa_star_table.rds"))
cat(sprintf("  Saved: kappa_star_table.rds (elapsed %.3f s)\n", t5_elapsed))

# ------------------------------------------------------------
# Table 4 data: real trial example based on Doody et al. (2014),
# EXPEDITION 1 of solanezumab (NEJM 2014; 370: 311-321), Protocol
# H8A-MC-LZAM Amendment (b), Section 12.1.1.
#
# The protocol prescribed 500 patients per group with co-primary
# endpoints ADAS-Cog11 and ADCS-ADL, assuming independence between
# the two endpoints.  We re-derive the planned sample size with our
# closed-form formula at rho = 0 to confirm the match, then vary
# rho to quantify the design-time benefit of accounting for the
# between-endpoint correlation.
# ------------------------------------------------------------
cat("\n--- Computing Table 4 (real trial example) ---\n")
t6_start <- Sys.time()

trial_delta1 <- 1.8   # ADAS-Cog11 treatment difference at 18 months
trial_delta2 <- 3.1   # ADCS-ADL  treatment difference at 18 months
trial_sd1    <- 9     # standard deviation of ADAS-Cog11 change
trial_sd2    <- 12    # standard deviation of ADCS-ADL  change
trial_alpha  <- 0.025 # one-sided significance level
trial_beta   <- 0.13  # overall power 0.87 = 0.89 * 0.98
trial_r      <- 1     # allocation ratio (1:1)
trial_rhos   <- c(0.0, 0.3, 0.5, 0.8)

rows_t4 <- list()
for (rho in trial_rhos) {
  ss_conv <- twocpcont_ss(
    delta1 = trial_delta1, delta2 = trial_delta2,
    sd1 = trial_sd1, sd2 = trial_sd2, rho = rho, r = trial_r,
    alpha = trial_alpha, beta = trial_beta,
    method = "bivariate"
  )
  ss_prop <- twocpcont_ss(
    delta1 = trial_delta1, delta2 = trial_delta2,
    sd1 = trial_sd1, sd2 = trial_sd2, rho = rho, r = trial_r,
    alpha = trial_alpha, beta = trial_beta,
    method = "univariate", gl_nodes = gl_nodes
  )
  rows_t4[[length(rows_t4) + 1L]] <- list(
    rho    = rho,
    N_conv = ss_conv$N,
    N_prop = ss_prop$N,
    diff   = ss_prop$N - ss_conv$N
  )
}
df_t4 <- do.call(rbind, lapply(rows_t4, as.data.frame))

N_baseline <- df_t4$N_prop[df_t4$rho == 0]
df_t4$reduction_pct <- (N_baseline - df_t4$N_prop) / N_baseline * 100

kappa_trial <- max(trial_delta1 / trial_sd1, trial_delta2 / trial_sd2) /
  min(trial_delta1 / trial_sd1, trial_delta2 / trial_sd2)
R_max_trial <- r_max(kappa = kappa_trial, alpha = trial_alpha,
                     beta = trial_beta) * 100
kappa_star_trial_eps01 <- kappa_star(epsilon = 0.01, alpha = trial_alpha,
                                     beta = trial_beta)
kappa_star_trial_eps05 <- kappa_star(epsilon = 0.05, alpha = trial_alpha,
                                     beta = trial_beta)

t6_elapsed <- as.numeric(difftime(Sys.time(), t6_start, units = "secs"))

table4_obj <- list(
  data         = df_t4,
  delta1       = trial_delta1,
  delta2       = trial_delta2,
  sd1          = trial_sd1,
  sd2          = trial_sd2,
  alpha        = trial_alpha,
  beta         = trial_beta,
  r            = trial_r,
  rho_set      = trial_rhos,
  kappa        = kappa_trial,
  R_max_pct    = R_max_trial,
  kappa_star_eps01 = kappa_star_trial_eps01,
  kappa_star_eps05 = kappa_star_trial_eps05,
  gl_nodes     = gl_nodes,
  elapsed      = t6_elapsed,
  computed     = Sys.time(),
  R_version    = R.version.string
)
saveRDS(table4_obj, out_path("table4_real_example.rds"))
cat(sprintf("  Saved: table4_real_example.rds (elapsed %.3f s)\n",
            t6_elapsed))

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------
cat("\n----- Summary -----\n")
for (pw in power_set) {
  d <- df_known[df_known$power == pw, ]
  cat(sprintf(paste0("Table 2 (1-beta = %.2f): match (diff = 0) %d / %d, ",
                     "achieved power range [%.7f, %.7f]\n"),
              pw, sum(d$diff == 0L), nrow(d), min(d$pow), max(d$pow)))
}
d_short <- df_known[df_known$pow < df_known$power, ]
if (nrow(d_short) > 0L) {
  cat("Table 2: cells in which N_prop falls short of the target power:\n")
  print(d_short, row.names = FALSE, digits = 8)
}
cat(sprintf("R_max range: kappa = 1, 1-beta = 0.80 -> %.2f%%; ",
            df_rmax$R_max[df_rmax$kappa == 1.0 &
                            df_rmax$power == 0.80] * 100))
df_rmax_80 <- df_rmax[df_rmax$power == 0.80, ]
cat(sprintf("kappa = 1.5 -> %.2f%%\n",
            df_rmax_80$R_max[which.min(abs(df_rmax_80$kappa - 1.5))] * 100))
cat(sprintf("kappa_star(eps = 0.05, 1-beta = 0.80) = %.4f\n",
            df_kstar$kappa_star[df_kstar$power == 0.80 &
                                  df_kstar$epsilon == 0.05]))
cat(sprintf("Table 4: rho = 0 reproduces N = %d (Doody 2014 planned 1000)\n",
            df_t4$N_prop[df_t4$rho == 0]))
cat(sprintf("Table 4: kappa = %.3f, R_max upper bound = %.2f%%\n",
            kappa_trial, R_max_trial))
cat(sprintf("Table 4: kappa_star(eps = 0.01, 1-beta = 0.87) = %.4f\n",
            kappa_star_trial_eps01))
cat(sprintf("Table 4: kappa_star(eps = 0.05, 1-beta = 0.87) = %.4f\n",
            kappa_star_trial_eps05))

cat(sprintf("\nAll .rds files saved to: %s/\n", output_dir))
cat("Run create_table_and_figure_manuscript.R next.\n")
