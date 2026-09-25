# ============================================================
# In-text numerical-claim verification for Sections 3 and 4
# ------------------------------------------------------------
# Verifies every concrete number that appears in the body text
# of Section 3 (Numerical investigation) and Section 4 (Real
# trial example) of the manuscript.  Table and figure contents
# themselves are produced by run_table_and_figure_manuscript.R
# and are not re-checked here; only the numbers cited inside
# the prose of those two sections are tested.
#
# Each claim is a separate PASS / FAIL line so that any drift
# between manuscript prose and code output is immediately
# visible.
#
# How to use:
#   Set the working directory to a writable folder that contains
#   this script.  If the R source files of twocpcont are also in
#   that folder, they are used; otherwise the installed twocpcont
#   package is used.
#
# Output:
#   In addition to printing to the console, a copy of the full
#   PASS / FAIL log is written to
#   "verify_intext_numbers_manuscript/verify_intext_numbers_log.txt"
#   (the folder is created automatically if it does not exist).
# ============================================================

# ------------------------------------------------------------
# Output directory for the log file
# ------------------------------------------------------------
output_dir <- "verify_intext_numbers_manuscript"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
out_path <- function(name) file.path(output_dir, name)

# Open a fresh log file (every cat() below also writes to it)
log_path <- out_path("verify_intext_numbers_log.txt")
cat("", file = log_path, append = FALSE)
log_msg <- function(msg) cat(msg, file = log_path, append = TRUE)

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

# Track results
results <- list()
record  <- function(tag, ok, detail = "") {
  results[[length(results) + 1L]] <<- list(tag = tag, ok = ok, detail = detail)
  msg <- sprintf("  [%s] %s  %s\n",
                 if (ok) "PASS" else "FAIL", tag, detail)
  cat(msg)
  log_msg(msg)
}

section <- function(s) {
  msg <- sprintf("\n==============================================================\n%s\n==============================================================\n", s)
  cat(msg)
  log_msg(msg)
}

# Settings shared across Section 3
alpha <- 0.025
beta  <- 0.20

# ============================================================
# SECTION 3.1: Setting
# ============================================================
section("Section 3.1 (Setting): effect-size ratios kappa per pattern")

# (3.1-1) kappa = 1.0 for (0.3, 0.3, 1, 1)
k_a <- max(0.3 / 1, 0.3 / 1) / min(0.3 / 1, 0.3 / 1)
record("(3.1-1) kappa = 1.0 for (delta1, delta2, sd1, sd2) = (0.3, 0.3, 1, 1)",
       abs(k_a - 1.0) < 0.05,
       sprintf("kappa = %.4f", k_a))

# (3.1-2) kappa = 1.4 for (0.7, 0.5, 1, 1)
k_b <- max(0.7 / 1, 0.5 / 1) / min(0.7 / 1, 0.5 / 1)
record("(3.1-2) kappa = 1.4 for (delta1, delta2, sd1, sd2) = (0.7, 0.5, 1, 1)",
       abs(k_b - 1.4) < 0.05,
       sprintf("kappa = %.4f", k_b))

# (3.1-3) kappa = 1.5 for (0.3, 0.3, 1, 1.5)
k_c <- max(0.3 / 1, 0.3 / 1.5) / min(0.3 / 1, 0.3 / 1.5)
record("(3.1-3) kappa = 1.5 for (delta1, delta2, sd1, sd2) = (0.3, 0.3, 1, 1.5)",
       abs(k_c - 1.5) < 0.05,
       sprintf("kappa = %.4f", k_c))

# (3.1-4) kappa = 2.1 for (0.7, 0.5, 1, 1.5)
k_d <- max(0.7 / 1, 0.5 / 1.5) / min(0.7 / 1, 0.5 / 1.5)
record("(3.1-4) kappa = 2.1 for (delta1, delta2, sd1, sd2) = (0.7, 0.5, 1, 1.5)",
       abs(k_d - 2.1) < 0.05,
       sprintf("kappa = %.4f", k_d))

# ============================================================
# Helper: compute Table 2 grid for downstream checks
# ============================================================
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

rows <- list()
for (sc in scenarios) {
  for (rho in rho_vals) {
    ss_p <- twocpcont_ss(delta1 = sc$delta1, delta2 = sc$delta2,
                         sd1 = sc$sd1, sd2 = sc$sd2, rho = rho,
                         r = sc$r, alpha = alpha, beta = beta,
                         method = "univariate")
    pw_p <- twocpcont_power(n1 = ss_p$n1, n2 = ss_p$n2,
                            delta1 = sc$delta1, delta2 = sc$delta2,
                            sd1 = sc$sd1, sd2 = sc$sd2, rho = rho,
                            alpha = alpha, method = "bivariate")
    rows[[length(rows) + 1L]] <- list(
      delta1 = sc$delta1, delta2 = sc$delta2,
      sd1 = sc$sd1, sd2 = sc$sd2, r = sc$r, rho = rho,
      n1 = ss_p$n1, n2 = ss_p$n2, N = ss_p$N,
      power = pw_p$powerCoprimary
    )
  }
}
df32 <- do.call(rbind, lapply(rows, as.data.frame))

# Helper: locate one row of df32
get_row <- function(d1, d2, s1, s2, rr, rho_val) {
  idx <- df32$delta1 == d1 & df32$delta2 == d2 &
         df32$sd1 == s1 & df32$sd2 == s2 &
         df32$r == rr & abs(df32$rho - rho_val) < 1e-9
  df32[idx, , drop = FALSE]
}

# ============================================================
# SECTION 3.2: Sample size comparison and achieved power
# ============================================================
section("Section 3.2 (Sample size comparison and achieved power)")

# (3.2-1) Power range approximately [0.8000, 0.8078]
pw_min <- min(df32$power)
pw_max <- max(df32$power)
record("(3.2-1) Achieved power range [0.8000, 0.8078]",
       abs(pw_min - 0.8000) < 0.0001 && abs(pw_max - 0.8078) < 0.0001,
       sprintf("[%.4f, %.4f]", pw_min, pw_max))

# (3.2-2) Maximum power 0.8078 occurs at
#         (delta1, delta2, sd1, sd2, r, rho) = (0.7, 0.5, 1, 1, 2, 0.5)
row_argmax <- df32[which.max(df32$power), ]
ok_max <- (row_argmax$delta1 == 0.7 && row_argmax$delta2 == 0.5 &&
           row_argmax$sd1    == 1   && row_argmax$sd2    == 1   &&
           row_argmax$r      == 2   && abs(row_argmax$rho - 0.5) < 1e-9 &&
           abs(row_argmax$power - 0.8078) < 0.0001)
record("(3.2-2) Maximum power 0.8078 at (0.7, 0.5, 1, 1, 2, 0.5)",
       ok_max,
       sprintf("(d1, d2, s1, s2, r, rho) = (%.1f, %.1f, %.0f, %.0f, %d, %.1f), power = %.4f",
               row_argmax$delta1, row_argmax$delta2,
               row_argmax$sd1, row_argmax$sd2,
               row_argmax$r, row_argmax$rho, row_argmax$power))

# (3.2-3) At the argmax, n2 = 49
record("(3.2-3) At the max-power cell, n2 = 49",
       row_argmax$n2 == 49L,
       sprintf("n2 = %d", row_argmax$n2))

# (3.2-4) Pattern (0.3, 0.3, 1, 1, 1): rho = 0 -> N = 460
row_460 <- get_row(0.3, 0.3, 1, 1, 1, 0.0)
record("(3.2-4) Pattern (0.3, 0.3, 1, 1, 1) at rho = 0 gives N = 460",
       nrow(row_460) == 1L && row_460$N == 460L,
       sprintf("N = %d", row_460$N))

# (3.2-5) Same pattern: rho = 0.8 -> N = 408
row_408 <- get_row(0.3, 0.3, 1, 1, 1, 0.8)
record("(3.2-5) Pattern (0.3, 0.3, 1, 1, 1) at rho = 0.8 gives N = 408",
       nrow(row_408) == 1L && row_408$N == 408L,
       sprintf("N = %d", row_408$N))

# (3.2-6) Reduction from 460 to 408 ~ 11%
red_pct_a <- (460 - 408) / 460 * 100
record("(3.2-6) Reduction from 460 to 408 is approximately 11%",
       abs(red_pct_a - 11) < 1,
       sprintf("reduction = %.2f%%", red_pct_a))

# (3.2-7) Pattern (0.7, 0.5, 1, 1.5, 1) at all 4 rho values -> N = 284
rows_kp <- df32[df32$delta1 == 0.7 & df32$delta2 == 0.5 &
                df32$sd1 == 1 & df32$sd2 == 1.5 & df32$r == 1, ]
record("(3.2-7) Pattern (0.7, 0.5, 1, 1.5, 1) yields N = 284 at all four rho",
       nrow(rows_kp) == 4L && all(rows_kp$N == 284L),
       sprintf("N values: %s",
               paste(rows_kp$N, collapse = ", ")))

# (3.2-8) Intermediate patterns: kappa = 1.4 -> reduction ~ 3-4%
get_red <- function(d1, d2, s1, s2, rr) {
  r0  <- get_row(d1, d2, s1, s2, rr, 0.0)$N
  r08 <- get_row(d1, d2, s1, s2, rr, 0.8)$N
  (r0 - r08) / r0 * 100
}
red_k14_a <- get_red(0.7, 0.5, 1, 1, 1)
red_k14_b <- get_red(0.7, 0.5, 1, 1, 2)
record("(3.2-8) kappa = 1.4 patterns: reduction in 3-4% range",
       red_k14_a >= 2.5 && red_k14_a <= 4.5 &&
       red_k14_b >= 2.5 && red_k14_b <= 4.5,
       sprintf("(0.7, 0.5, 1, 1, 1): %.2f%%, (0.7, 0.5, 1, 1, 2): %.2f%%",
               red_k14_a, red_k14_b))

# (3.2-9) kappa = 1.5 patterns: reduction ~ 2%
red_k15_a <- get_red(0.3, 0.3, 1, 1.5, 1)
red_k15_b <- get_red(0.3, 0.3, 1, 1.5, 2)
record("(3.2-9) kappa = 1.5 patterns: reduction approximately 2%",
       red_k15_a >= 1.0 && red_k15_a <= 3.0 &&
       red_k15_b >= 1.0 && red_k15_b <= 3.0,
       sprintf("(0.3, 0.3, 1, 1.5, 1): %.2f%%, (0.3, 0.3, 1, 1.5, 2): %.2f%%",
               red_k15_a, red_k15_b))

# (3.2-10) Monotonically non-increasing in rho for every pattern
all_mono <- all(sapply(scenarios, function(sc) {
  Ns <- sapply(rho_vals, function(rr) {
    get_row(sc$delta1, sc$delta2, sc$sd1, sc$sd2, sc$r, rr)$N
  })
  all(diff(Ns) <= 0L)
}))
record("(3.2-10) N is monotonically non-increasing in rho for all 8 patterns",
       all_mono,
       "")

# ============================================================
# SECTION 3.3: R_max and kappa_star
# ============================================================
section("Section 3.3 (R_max and kappa_star)")

# (3.3-1) R_max(kappa = 1, 1 - beta = 0.70) = 28.63%
Rmax_1_pwr70 <- r_max(kappa = 1, alpha = 0.025, beta = 0.30) * 100
record("(3.3-1) R_max(1, 1-beta = 0.70) approximately 28.63%",
       abs(Rmax_1_pwr70 - 28.63) < 0.01,
       sprintf("R_max = %.4f%%", Rmax_1_pwr70))

# (3.3-2) R_max(kappa = 1, 1 - beta = 0.80) = 23.85%
Rmax_1_pwr80 <- r_max(kappa = 1, alpha = 0.025, beta = 0.20) * 100
record("(3.3-2) R_max(1, 1-beta = 0.80) approximately 23.85%",
       abs(Rmax_1_pwr80 - 23.85) < 0.01,
       sprintf("R_max = %.4f%%", Rmax_1_pwr80))

# (3.3-3) R_max(kappa = 1, 1 - beta = 0.90) = 18.57%
Rmax_1_pwr90 <- r_max(kappa = 1, alpha = 0.025, beta = 0.10) * 100
record("(3.3-3) R_max(1, 1-beta = 0.90) approximately 18.57%",
       abs(Rmax_1_pwr90 - 18.57) < 0.01,
       sprintf("R_max = %.4f%%", Rmax_1_pwr90))

# (3.3-4) R_max crosses 5% at kappa approximately 1.50 (beta = 0.30)
ks_5pct_b30 <- kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.30)
record("(3.3-4) R_max(kappa, 1-beta = 0.70) = 5% at kappa approximately 1.50",
       abs(ks_5pct_b30 - 1.50) < 0.01,
       sprintf("kappa_star = %.4f", ks_5pct_b30))

# (3.3-5) R_max crosses 5% at kappa approximately 1.37 (beta = 0.20)
ks_5pct_b20 <- kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.20)
record("(3.3-5) R_max(kappa, 1-beta = 0.80) = 5% at kappa approximately 1.37",
       abs(ks_5pct_b20 - 1.37) < 0.01,
       sprintf("kappa_star = %.4f", ks_5pct_b20))

# (3.3-6) R_max crosses 5% at kappa approximately 1.24 (beta = 0.10)
ks_5pct_b10 <- kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.10)
record("(3.3-6) R_max(kappa, 1-beta = 0.90) = 5% at kappa approximately 1.24",
       abs(ks_5pct_b10 - 1.24) < 0.01,
       sprintf("kappa_star = %.4f", ks_5pct_b10))

# (3.3-7) R_max crosses 1% at kappa approximately 1.79 (beta = 0.30)
ks_1pct_b30 <- kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.30)
record("(3.3-7) R_max(kappa, 1-beta = 0.70) = 1% at kappa approximately 1.79",
       abs(ks_1pct_b30 - 1.79) < 0.01,
       sprintf("kappa_star = %.4f", ks_1pct_b30))

# (3.3-8) R_max crosses 1% at kappa approximately 1.61 (beta = 0.20)
ks_1pct_b20 <- kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.20)
record("(3.3-8) R_max(kappa, 1-beta = 0.80) = 1% at kappa approximately 1.61",
       abs(ks_1pct_b20 - 1.61) < 0.01,
       sprintf("kappa_star = %.4f", ks_1pct_b20))

# (3.3-9) R_max crosses 1% at kappa approximately 1.44 (beta = 0.10)
ks_1pct_b10 <- kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.10)
record("(3.3-9) R_max(kappa, 1-beta = 0.90) = 1% at kappa approximately 1.44",
       abs(ks_1pct_b10 - 1.44) < 0.01,
       sprintf("kappa_star = %.4f", ks_1pct_b10))

# (3.3-10) Table 3 reference: kappa_star(eps = 0.005, 1-beta = 0.80) = 1.700
ks_t3_a <- kappa_star(epsilon = 0.005, alpha = 0.025, beta = 0.20)
record("(3.3-10) kappa_star(eps = 0.005, 1-beta = 0.80) approximately 1.700",
       abs(ks_t3_a - 1.700) < 0.001,
       sprintf("kappa_star = %.4f", ks_t3_a))

# (3.3-11) kappa_star(eps = 0.050, 1-beta = 0.80) = 1.369
ks_t3_b <- kappa_star(epsilon = 0.050, alpha = 0.025, beta = 0.20)
record("(3.3-11) kappa_star(eps = 0.050, 1-beta = 0.80) approximately 1.369",
       abs(ks_t3_b - 1.369) < 0.001,
       sprintf("kappa_star = %.4f", ks_t3_b))

# (3.3-12) kappa_star(eps = 0.005, 1-beta = 0.90) = 1.512
ks_t3_c <- kappa_star(epsilon = 0.005, alpha = 0.025, beta = 0.10)
record("(3.3-12) kappa_star(eps = 0.005, 1-beta = 0.90) approximately 1.512",
       abs(ks_t3_c - 1.512) < 0.001,
       sprintf("kappa_star = %.4f", ks_t3_c))

# (3.3-13) kappa_star(eps = 0.050, 1-beta = 0.90) = 1.239
ks_t3_d <- kappa_star(epsilon = 0.050, alpha = 0.025, beta = 0.10)
record("(3.3-13) kappa_star(eps = 0.050, 1-beta = 0.90) approximately 1.239",
       abs(ks_t3_d - 1.239) < 0.001,
       sprintf("kappa_star = %.4f", ks_t3_d))

# (3.3-14) All 12 (epsilon, 1 - beta) combinations: kappa_star in [1.24, 1.89]
eps_set <- c(0.005, 0.010, 0.020, 0.050)
pwr_set <- c(0.70, 0.80, 0.90)
all_ks <- numeric(0)
for (p in pwr_set) {
  for (e in eps_set) {
    all_ks <- c(all_ks, kappa_star(epsilon = e, alpha = 0.025, beta = 1 - p))
  }
}
record("(3.3-14) Across all 12 (eps, 1-beta), kappa_star in [1.24, 1.89]",
       min(all_ks) >= 1.235 && max(all_ks) <= 1.895,
       sprintf("range [%.4f, %.4f]", min(all_ks), max(all_ks)))

# (3.3-15) Universal upper bound R_max(1, 1-beta = 0.80) approximately 23.85%
record("(3.3-15) Universal upper bound R_max(1, 1-beta = 0.80) = 23.85% (same as 3.3-2)",
       abs(Rmax_1_pwr80 - 23.85) < 0.01,
       sprintf("R_max = %.4f%%", Rmax_1_pwr80))

# ============================================================
# SECTION 4: Real trial example (Doody / EXPEDITION 1)
# ============================================================
section("Section 4 (Real trial example)")

# (4-1) delta1 / sd1 = 0.200
record("(4-1) delta1 / sd1 = 1.8 / 9 = 0.200",
       abs(1.8 / 9 - 0.200) < 1e-4,
       sprintf("ratio = %.4f", 1.8 / 9))

# (4-2) delta2 / sd2 = 0.258 (with the manuscript rounding 0.258)
record("(4-2) delta2 / sd2 = 3.1 / 12 approximately 0.258",
       abs(3.1 / 12 - 0.258) < 1e-3,
       sprintf("ratio = %.6f", 3.1 / 12))

# (4-3) kappa = 1.29
kappa_doody <- (3.1 / 12) / (1.8 / 9)
record("(4-3) Doody example: kappa = 1.29",
       abs(kappa_doody - 1.29) < 0.005,
       sprintf("kappa = %.4f", kappa_doody))

# (4-4) Marginal-power product 0.89 * 0.98 = 0.87 (within rounding)
record("(4-4) 0.89 * 0.98 approximately 0.87 (within 0.01)",
       abs(0.89 * 0.98 - 0.87) < 0.01,
       sprintf("product = %.4f", 0.89 * 0.98))

# (4-5) Doody example: rho = 0 yields N = 1000
ss_doody_0 <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1,
                           sd1 = 9, sd2 = 12, rho = 0, r = 1,
                           alpha = 0.025, beta = 0.13,
                           method = "univariate")
record("(4-5) Doody example: rho = 0 yields N = 1000",
       ss_doody_0$N == 1000L,
       sprintf("N = %d", ss_doody_0$N))

# (4-6) Doody example: rho = 0.8 yields N = 960
ss_doody_8 <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1,
                           sd1 = 9, sd2 = 12, rho = 0.8, r = 1,
                           alpha = 0.025, beta = 0.13,
                           method = "univariate")
record("(4-6) Doody example: rho = 0.8 yields N = 960",
       ss_doody_8$N == 960L,
       sprintf("N = %d", ss_doody_8$N))

# (4-7) Reduction of 40 patients (1000 - 960)
record("(4-7) Doody example: rho = 0 -> 0.8 reduction = 40 patients",
       (ss_doody_0$N - ss_doody_8$N) == 40L,
       sprintf("reduction = %d patients", ss_doody_0$N - ss_doody_8$N))

# (4-8) Approximately 4% reduction (40 / 1000)
red_doody <- (ss_doody_0$N - ss_doody_8$N) / ss_doody_0$N * 100
record("(4-8) Doody example: rho = 0 -> 0.8 reduction approximately 4%",
       abs(red_doody - 4) < 0.5,
       sprintf("reduction = %.2f%%", red_doody))

# (4-9) Monotonic decrease in N across rho = 0, 0.3, 0.5, 0.8
N_doody_seq <- sapply(c(0, 0.3, 0.5, 0.8), function(rr) {
  twocpcont_ss(delta1 = 1.8, delta2 = 3.1,
               sd1 = 9, sd2 = 12, rho = rr, r = 1,
               alpha = 0.025, beta = 0.13,
               method = "univariate")$N
})
record("(4-9) Doody example: N monotonically non-increasing in rho",
       all(diff(N_doody_seq) <= 0L),
       sprintf("N values: %s", paste(N_doody_seq, collapse = ", ")))

# (4-10) The two methods agree at all four rho values
N_doody_conv <- sapply(c(0, 0.3, 0.5, 0.8), function(rr) {
  twocpcont_ss(delta1 = 1.8, delta2 = 3.1,
               sd1 = 9, sd2 = 12, rho = rr, r = 1,
               alpha = 0.025, beta = 0.13,
               method = "bivariate")$N
})
record("(4-10) Doody example: bivariate and univariate agree at all four rho",
       all(N_doody_conv == N_doody_seq),
       sprintf("N_conv: %s; N_prop: %s",
               paste(N_doody_conv, collapse = ", "),
               paste(N_doody_seq, collapse = ", ")))

# (4-11) R_max(kappa = 1.29, 1-beta = 0.87) approximately 4.6%
Rmax_doody <- r_max(kappa = kappa_doody, alpha = 0.025,
                    beta = 0.13) * 100
record("(4-11) R_max(kappa = 1.29, 1-beta = 0.87) approximately 4.6%",
       abs(Rmax_doody - 4.6) < 0.1,
       sprintf("R_max = %.4f%%", Rmax_doody))

# (4-12) Integer reduction of 4.0% recovers roughly 87% of the analytic ceiling
#        i.e. 4.0 / R_max ~ 0.87
ratio_recover <- red_doody / Rmax_doody
record("(4-12) Integer reduction 4.0% recovers roughly 87% of analytic ceiling",
       abs(ratio_recover - 0.87) < 0.03,
       sprintf("4.0 / R_max = %.4f (target 0.87)", ratio_recover))

# (4-13) Remaining gap is approximately 0.6 percentage points
#        (R_max - integer reduction = 4.6 - 4.0 = 0.6)
gap_pp <- Rmax_doody - red_doody
record("(4-13) Remaining gap = R_max - 4.0% approximately 0.6 pp",
       abs(gap_pp - 0.6) < 0.1,
       sprintf("gap = %.4f pp (target 0.6)", gap_pp))

# (4-14) kappa_star(eps = 0.05, 1-beta = 0.87) approximately 1.28
ks_doody_05 <- kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.13)
record("(4-14) kappa_star(eps = 0.05, 1-beta = 0.87) approximately 1.28",
       abs(ks_doody_05 - 1.28) < 0.01,
       sprintf("kappa_star = %.4f", ks_doody_05))

# (4-15) kappa_star(eps = 0.01, 1-beta = 0.87) approximately 1.49
ks_doody_01 <- kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.13)
record("(4-15) kappa_star(eps = 0.01, 1-beta = 0.87) approximately 1.49",
       abs(ks_doody_01 - 1.49) < 0.01,
       sprintf("kappa_star = %.4f", ks_doody_01))

# (4-16) kappa_star(0.05) < kappa < kappa_star(0.01)
record("(4-16) Doody example: kappa_star(0.05) < kappa < kappa_star(0.01)",
       ks_doody_05 < kappa_doody && kappa_doody < ks_doody_01,
       sprintf("%.4f < %.4f < %.4f",
               ks_doody_05, kappa_doody, ks_doody_01))

# (4-17) Bound on gain: approximately 46 patients out of 1000
#        Manuscript states "approximately 46 patients out of 1000"
#        i.e. R_max * 1000 rounds to ~46
bound_patients <- Rmax_doody / 100 * 1000
record("(4-17) Gain bound: R_max * 1000 approximately 46 patients",
       abs(bound_patients - 46) < 1.0,
       sprintf("R_max * 1000 = %.2f patients (target 46)", bound_patients))

# (4-18) Manuscript claims "deviates ... by less than 5%"
#        Check that R_max < 5%
record("(4-18) R_max < 5% (so 'less than 5%' is correct)",
       Rmax_doody < 5.0,
       sprintf("R_max = %.4f%% (less than 5%%)", Rmax_doody))

# ============================================================
# Final summary
# ============================================================
section("Summary")

n_total <- length(results)
n_pass  <- sum(sapply(results, function(z) z$ok))
n_fail  <- n_total - n_pass

msg <- sprintf("\n  TOTAL  : %d\n  PASSED : %d\n  FAILED : %d\n",
               n_total, n_pass, n_fail)
cat(msg); log_msg(msg)

if (n_fail > 0L) {
  msg <- "\n  Failing tests:\n"
  cat(msg); log_msg(msg)
  for (z in results) {
    if (!z$ok) {
      msg <- sprintf("    %s  (%s)\n", z$tag, z$detail)
      cat(msg); log_msg(msg)
    }
  }
} else {
  msg <- "\n  All in-text numerical claims match the code output.\n"
  cat(msg); log_msg(msg)
}

cat(sprintf("\n  Log written to %s\n\n", log_path))
