# ============================================================
# In-text numerical-claim verification
# ------------------------------------------------------------
# Verifies every concrete number that appears in the body text
# of Section 2.4 (reduction in sample size due to correlation),
# Section 3 (numerical investigation) and Section 4 (real trial
# example) of the manuscript.  The Discussion mostly repeats numbers
# from Sections 3 and 4, which are covered by those checks; its one
# new number is checked in the Section 5 block.
# Table and figure contents themselves are produced by
# run_table_and_figure_manuscript.R and are not re-checked here.
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

# Settings shared across Sections 2 and 3
alpha     <- 0.025
power_set <- c(0.80, 0.85, 0.90)

# ============================================================
# SECTION 1: Introduction
# ============================================================
section("Section 1 (Introduction)")

record("(1-1) Product of two marginal powers 0.8^2 = 0.64",
       abs(0.8 ^ 2 - 0.64) < 1e-12,
       sprintf("0.8^2 = %.4f", 0.8 ^ 2))

# ============================================================
# SECTION 2.3 and Appendix B: accuracy of the 5-point rule
# ============================================================
section("Section 2.3 and Appendix B (Gauss-Legendre approximation)")

# (2.3-1) |I - I_GL| < 1e-5 for |rho| <= 0.8.  The exact Plackett integral
#         is computed by adaptive quadrature over a grid of (a, b, rho)
#         that covers the arguments arising in this article.
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
record("(2.3-1) 5-point rule: |I - I_GL| < 1e-5 on the 13209-point grid of Appendix B",
       max(gl_err) < 1e-5 && nrow(gl_grid) == 13209L,
       sprintf("grid points = %d, max error = %.2e", nrow(gl_grid),
               max(gl_err)))

# (2.3-2), (2.3-3) Anchor weight w* of eq. (8) against the exact root of
#         F(w; kappa) = 0 in eq. (A2).  w* exceeds the exact root for
#         1 < kappa <= 3 by at most about 0.09, and w* reaches one at
#         kappa of about 1.41, 1.36 and 1.30 for 1-beta = 0.80, 0.85, 0.90.
z_a_w <- qnorm(1 - alpha)
w_star_fun <- function(k, tp) {
  p   <- sqrt(tp)
  u_p <- qnorm(p)
  0.5 + (k - 1) * (u_p + z_a_w) * dnorm(u_p) /
    ((1 + k) * abs(2 * log(p)) * p)
}
w_root_fun <- function(k, tp) {
  if (k == 1) return(0.5)
  f_w <- function(w) {
    qnorm(tp ^ (1 - w)) - k * qnorm(tp ^ w) - (k - 1) * z_a_w
  }
  # F is increasing in w and tends to +Inf as w -> 1, but in double
  # precision F(1 - 1e-12) can still be negative for large kappa; the
  # root is then within 1e-12 of one.
  if (f_w(1 - 1e-12) <= 0) return(1)
  uniroot(f_w, lower = 0.5, upper = 1 - 1e-12, tol = 1e-13)$root
}
k_one <- sapply(power_set, function(tp) {
  uniroot(function(k) w_star_fun(k, tp) - 1, lower = 1, upper = 3,
          tol = 1e-10)$root
})
record("(2.3-2) w* reaches one at kappa = 1.41, 1.36, 1.30 (0.80, 0.85, 0.90)",
       all(abs(round(k_one, 2) - c(1.41, 1.36, 1.30)) < 1e-9),
       sprintf("kappa = %s", paste(sprintf("%.4f", k_one), collapse = ", ")))
k_seq_w <- seq(1.001, 3, by = 0.001)
w_gap <- sapply(power_set, function(tp) {
  g <- sapply(k_seq_w, function(k) {
    min(1, w_star_fun(k, tp)) - w_root_fun(k, tp)
  })
  c(min(g), max(g))
})
record("(2.3-3) 1 < kappa <= 3: w* exceeds the exact root, by at most about 0.09",
       all(w_gap[1, ] > -1e-8) && abs(round(max(w_gap[2, ]), 2) - 0.09) < 1e-9,
       sprintf("min gap = %.2e, max gap = %s", min(w_gap[1, ]),
               paste(sprintf("%.4f", w_gap[2, ]), collapse = ", ")))

# ============================================================
# SECTION 2.4: Maximum reduction rate and threshold
# ============================================================
section("Section 2.4 (Maximum reduction rate and threshold effect size ratio)")

# (2.4-1) Rmax(1, 0.20) approximately 23.8%
Rmax_1_80 <- r_max(kappa = 1, alpha = alpha, beta = 0.20) * 100
record("(2.4-1) Rmax(1, 1-beta = 0.80) approximately 23.8%",
       abs(Rmax_1_80 - 23.8) < 0.05,
       sprintf("Rmax = %.4f%%", Rmax_1_80))

# (2.4-2) The root of the rho = 0 equation lies in
#         [z_alpha + z_beta, z_alpha + u_p]
z_alpha <- qnorm(1 - alpha)
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
record("(2.4-2) Root of the rho = 0 equation lies in [z_a + z_b, z_a + u_p]",
       in_bracket, "kappa in [1, 3], 1-beta in {0.80, 0.85, 0.90}")

# (2.4-3) Rmax is non-increasing in kappa
mono_k <- all(sapply(1 - power_set, function(b) {
  v <- sapply(seq(1, 3, by = 0.01), r_max, alpha = alpha, beta = b)
  all(diff(v) <= 1e-12)
}))
record("(2.4-3) Rmax(kappa, beta) is non-increasing in kappa", mono_k,
       "kappa in [1, 3] by 0.01")

# (2.4-4) Closed-form and exact Rmax differ by less than 1e-4
diff_max <- max(sapply(c(0.30, 0.20, 0.15, 0.13, 0.10), function(b) {
  k_seq <- seq(1, 3, by = 0.001)
  max(abs(sapply(k_seq, r_max, alpha = alpha, beta = b, exact = TRUE) -
            sapply(k_seq, r_max, alpha = alpha, beta = b, exact = FALSE)))
}))
record("(2.4-4) Closed-form and exact Rmax differ by less than 1e-4",
       diff_max < 1e-4,
       sprintf("max |difference| = %.2e", diff_max))

# (2.4-5) kappa_star(0.05, 1-beta = 0.80) approximately 1.37
ks_05_80 <- kappa_star(epsilon = 0.05, alpha = alpha, beta = 0.20)
record("(2.4-5) kappa_star(eps = 0.05, 1-beta = 0.80) approximately 1.37",
       abs(ks_05_80 - 1.37) < 0.005,
       sprintf("kappa_star = %.4f", ks_05_80))

# (2.4-6) kappa_star is the exact inverse of Rmax
inv_err <- max(sapply(1 - power_set, function(b) {
  sapply(c(0.05, 0.10, 0.15), function(e) {
    abs(r_max(kappa = kappa_star(epsilon = e, alpha = alpha, beta = b),
              alpha = alpha, beta = b) - e)
  })
}))
record("(2.4-6) Rmax(kappa_star(eps)) = eps", inv_err < 1e-10,
       sprintf("max |Rmax - eps| = %.2e", inv_err))

# (2.4-7) kappa_star = 1 when eps >= Rmax(1, beta) and > 1 otherwise
ks_bound <- sapply(1 - power_set, function(b) {
  r1 <- r_max(kappa = 1, alpha = alpha, beta = b)
  c(kappa_star(r1, alpha = alpha, beta = b),
    kappa_star(min(0.99, r1 + 0.05), alpha = alpha, beta = b),
    kappa_star(r1 - 0.01, alpha = alpha, beta = b))
})
record("(2.4-7) kappa_star = 1 for eps >= Rmax(1, beta), > 1 for eps < Rmax(1, beta)",
       all(abs(ks_bound[1:2, ] - 1) < 1e-8) && all(ks_bound[3, ] > 1),
       sprintf("at Rmax(1): %s; just below: %s",
               paste(sprintf("%.6f", ks_bound[1, ]), collapse = ", "),
               paste(sprintf("%.4f", ks_bound[3, ]), collapse = ", ")))

# (2.4-8) The argument (1 - beta) / Phi(z_beta + nu) of eq. (20) lies in
#         (0, 1) for every eps in (0, 1)
arg_ok <- all(sapply(1 - power_set, function(b) {
  z_b <- qnorm(1 - b)
  sapply(seq(0.001, 0.999, by = 0.001), function(e) {
    nu  <- (z_alpha + z_b) * (1 / sqrt(1 - e) - 1)
    arg <- (1 - b) / pnorm(z_b + nu)
    arg > 0 && arg < 1
  })
}))
record("(2.4-8) Argument of the quantile function in eq. (20) lies in (0, 1)",
       arg_ok, "eps in [0.001, 0.999] by 0.001")

# ============================================================
# SECTION 3.1: Setting
# ============================================================
section("Section 3.1 (Setting)")

kappa_of <- function(d1, d2, s1, s2) {
  max(d1 / s1, d2 / s2) / min(d1 / s1, d2 / s2)
}
k_vals <- c(kappa_of(0.3, 0.3, 1, 1), kappa_of(0.7, 0.5, 1, 1),
            kappa_of(0.3, 0.3, 1, 1.5), kappa_of(0.7, 0.5, 1, 1.5))
record("(3.1-1) kappa = 1.0, 1.4, 1.5, 2.1 for the four (delta, sigma) pairs",
       all(abs(k_vals - c(1.0, 1.4, 1.5, 2.1)) < 0.05),
       sprintf("kappa = %s", paste(sprintf("%.4f", k_vals), collapse = ", ")))

# (3.1-2) Equal allocation minimizes N (1/n1 + 1/n2 at fixed N)
N_fix  <- 600
n1_seq <- 1:(N_fix - 1)
record("(3.1-2) 1/n1 + 1/n2 at fixed N is minimized at n1 = n2",
       n1_seq[which.min(1 / n1_seq + 1 / (N_fix - n1_seq))] == N_fix / 2,
       sprintf("N = %d", N_fix))

# ------------------------------------------------------------
# Helper: compute the 96-cell grid of Table 2
# ------------------------------------------------------------
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
for (pw in power_set) {
  for (sc in scenarios) {
    for (rho in rho_vals) {
      ss_p <- twocpcont_ss(delta1 = sc$delta1, delta2 = sc$delta2,
                           sd1 = sc$sd1, sd2 = sc$sd2, rho = rho,
                           r = sc$r, alpha = alpha, beta = 1 - pw,
                           method = "univariate")
      ss_c <- twocpcont_ss(delta1 = sc$delta1, delta2 = sc$delta2,
                           sd1 = sc$sd1, sd2 = sc$sd2, rho = rho,
                           r = sc$r, alpha = alpha, beta = 1 - pw,
                           method = "bivariate")
      pw_p <- twocpcont_power(n1 = ss_p$n1, n2 = ss_p$n2,
                              delta1 = sc$delta1, delta2 = sc$delta2,
                              sd1 = sc$sd1, sd2 = sc$sd2, rho = rho,
                              alpha = alpha, method = "bivariate")
      rows[[length(rows) + 1L]] <- list(
        power = pw, delta1 = sc$delta1, delta2 = sc$delta2,
        sd1 = sc$sd1, sd2 = sc$sd2, r = sc$r, rho = rho,
        n2 = ss_p$n2, N = ss_p$N, N_conv = ss_c$N,
        pow = pw_p$powerCoprimary
      )
    }
  }
}
grid <- do.call(rbind, lapply(rows, as.data.frame))

# Helper: N_prop for one pattern over rho
N_over_rho <- function(d1, d2, s1, s2, rr, pw) {
  d <- grid[grid$delta1 == d1 & grid$delta2 == d2 & grid$sd1 == s1 &
              grid$sd2 == s2 & grid$r == rr & grid$power == pw, ]
  d$N[order(d$rho)]
}
red_pct <- function(N_vec) (N_vec[1] - N_vec[length(N_vec)]) / N_vec[1] * 100

record("(3.1-3) The grid has 96 combinations", nrow(grid) == 96L,
       sprintf("rows = %d", nrow(grid)))

record("(3.1-4) Computations performed in R version 4.6.0",
       paste(R.version$major, R.version$minor, sep = ".") == "4.6.0",
       R.version.string)

# ============================================================
# SECTION 3.2: Sample size comparison and achieved power
# ============================================================
section("Section 3.2 (Sample size comparison and achieved power)")

n_match <- sapply(power_set, function(pw) {
  sum(grid$N[grid$power == pw] == grid$N_conv[grid$power == pw])
})
record("(3.2-1) Matches: 32 / 32 at 0.80 and 0.85, 31 / 32 at 0.90",
       all(n_match == c(32L, 32L, 31L)),
       sprintf("matches = %s", paste(n_match, collapse = ", ")))

mis <- grid[grid$N != grid$N_conv, ]
ok_mis <- nrow(mis) == 1L && mis$delta1 == 0.3 && mis$delta2 == 0.3 &&
  mis$sd1 == 1 && mis$sd2 == 1 && mis$r == 1 && mis$rho == 0.5 &&
  mis$power == 0.90 && mis$N == 556 && mis$N_conv == 558
record("(3.2-2) Mismatch at (0.3, 0.3, 1, 1, 1, 0.5), 0.90: N = 556 vs 558",
       ok_mis,
       sprintf("N_prop = %s, N_conv = %s",
               paste(mis$N, collapse = ", "),
               paste(mis$N_conv, collapse = ", ")))

# (3.2-3) Continuous n2: exact 278.0009, closed form 277.9970
lam_exact <- uniroot(function(lam) {
  pbivnorm::pbivnorm(lam - z_alpha, lam - z_alpha, rho = 0.5) - 0.90
}, lower = 2, upper = 5, tol = 1e-14)$root
n2_exact <- 2 / 0.3 ^ 2 * lam_exact ^ 2
# Closed-form continuous n2 (before the ceiling) returned by twocpcont_ss()
n2_cf <- twocpcont_ss(delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1,
                      rho = 0.5, r = 1, alpha = 0.025, beta = 0.10,
                      method = "univariate")$n2_cont
record("(3.2-3) Continuous n2: exact 278.0009, closed form 277.9970",
       abs(n2_exact - 278.0009) < 5e-5 && abs(n2_cf - 277.9970) < 5e-5,
       sprintf("exact = %.4f, closed form = %.4f", n2_exact, n2_cf))
record("(3.2-4) Approximation error in n2 about 0.004",
       abs((n2_exact - n2_cf) - 0.004) < 0.0005,
       sprintf("difference = %.5f", n2_exact - n2_cf))

record("(3.2-5) Achieved power at N = 556 is 0.8999989 (short by about 1e-6)",
       abs(mis$pow - 0.8999989) < 5e-8 && abs(0.90 - mis$pow - 1e-6) < 1e-7,
       sprintf("power = %.7f, shortfall = %.2e", mis$pow, 0.90 - mis$pow))

ok_rows <- grid[grid$N == grid$N_conv, ]
rng <- sapply(power_set, function(pw) {
  range(ok_rows$pow[ok_rows$power == pw])
})
record("(3.2-6) Other 95: achieved power >= target",
       all(ok_rows$pow >= ok_rows$power),
       sprintf("min(power - target) = %.2e", min(ok_rows$pow - ok_rows$power)))
record(paste0("(3.2-7) Ranges [0.8000, 0.8078], [0.8501, 0.8544], ",
              "[0.9000, 0.9040]"),
       all(abs(round(rng, 4) - rbind(c(0.8000, 0.8501, 0.9000),
                                     c(0.8078, 0.8544, 0.9040))) < 1e-9),
       sprintf("%s", paste(sprintf("[%.4f, %.4f]", rng[1, ], rng[2, ]),
                           collapse = ", ")))

N_k1_80 <- N_over_rho(0.3, 0.3, 1, 1, 1, 0.80)
N_k1_90 <- N_over_rho(0.3, 0.3, 1, 1, 1, 0.90)
record("(3.2-8) kappa = 1, 1-beta = 0.80: 460 -> 408 (about 11%)",
       N_k1_80[1] == 460 && N_k1_80[4] == 408 &&
         abs(red_pct(N_k1_80) - 11) < 0.5,
       sprintf("%d -> %d, %.2f%%", N_k1_80[1], N_k1_80[4], red_pct(N_k1_80)))
record("(3.2-9) kappa = 1, 1-beta = 0.90: 574 -> 530 (about 8%)",
       N_k1_90[1] == 574 && N_k1_90[4] == 530 &&
         abs(red_pct(N_k1_90) - 8) < 0.5,
       sprintf("%d -> %d, %.2f%%", N_k1_90[1], N_k1_90[4], red_pct(N_k1_90)))

const_k21 <- all(sapply(power_set, function(pw) {
  length(unique(N_over_rho(0.7, 0.5, 1, 1.5, 1, pw))) == 1L
}))
record("(3.2-10) kappa = 2.1, r = 1: N constant in rho at all three powers",
       const_k21, "")

red_mid <- unlist(lapply(power_set, function(pw) {
  c(red_pct(N_over_rho(0.7, 0.5, 1, 1, 1, pw)),
    red_pct(N_over_rho(0.7, 0.5, 1, 1, 2, pw)),
    red_pct(N_over_rho(0.3, 0.3, 1, 1.5, 1, pw)),
    red_pct(N_over_rho(0.3, 0.3, 1, 1.5, 2, pw)))
}))
record("(3.2-11) kappa = 1.4 and 1.5: reduction at most 4%",
       max(red_mid) <= 4 + 1e-9,
       sprintf("max reduction = %.2f%%", max(red_mid)))

mono_all <- all(sapply(power_set, function(pw) {
  all(sapply(scenarios, function(sc) {
    all(diff(N_over_rho(sc$delta1, sc$delta2, sc$sd1, sc$sd2, sc$r,
                        pw)) <= 0)
  }))
}))
record("(3.2-12) N_prop non-increasing in rho in all 24 series", mono_all, "")

r2 <- grid[grid$r == 2, ]
record("(3.2-13) All 48 combinations with r = 2 give identical N",
       nrow(r2) == 48L && all(r2$N == r2$N_conv),
       sprintf("rows = %d, matches = %d", nrow(r2), sum(r2$N == r2$N_conv)))

# ============================================================
# SECTION 3.3: Comparison with published tables and computational cost
# (reads the rds files written by run_table_and_figure_manuscript.R)
# ============================================================
section("Section 3.3 (Comparison with published tables and computational cost)")

pub <- readRDS(file.path("table_and_figure_manuscript",
                         "published_table_comparison.rds"))$sozu
lo  <- pub$rho <= 0.8
hi  <- !lo

record("(3.3-1) Exact root reproduces all 210 tabulated values",
       nrow(pub) == 210L && all(abs(pub$dev_exact) < 1e-9),
       sprintf("rows = %d, matches = %d", nrow(pub),
               sum(abs(pub$dev_exact) < 1e-9)))

record("(3.3-2) rho <= 0.8: closed form agrees with 177 of 180, others differ by 0.001",
       sum(lo) == 180L && sum(abs(pub$dev_cf[lo]) < 1e-9) == 177L &&
         all(abs(abs(pub$dev_cf[lo & abs(pub$dev_cf) > 1e-9]) - 0.001) < 1e-9),
       sprintf("rows = %d, matches = %d", sum(lo),
               sum(abs(pub$dev_cf[lo]) < 1e-9)))

max_err <- max(abs(pub$err_cf[lo]))
record("(3.3-3) rho <= 0.8: largest |lambda_w* - exact root| below 4e-4",
       max_err < 4e-4,
       sprintf("max error = %.3e", max_err))

rel_n2 <- max(abs(((pub$C2_cf[lo] + z_alpha) /
                     (pub$C2_exact[lo] + z_alpha)) ^ 2 - 1)) * 100
record("(3.3-4) rho <= 0.8: c (lambda_w*)^2 changes by less than 0.03%",
       rel_n2 < 0.03,
       sprintf("max relative change = %.4f%%", rel_n2))

record("(3.3-5) rho = 0.95: 20 of 30 agree, largest difference 0.002",
       sum(hi) == 30L && sum(abs(pub$dev_cf[hi]) < 1e-9) == 20L &&
         abs(max(abs(pub$dev_cf[hi])) - 0.002) < 1e-9,
       sprintf("rows = %d, matches = %d, max |dev| = %.3f", sum(hi),
               sum(abs(pub$dev_cf[hi]) < 1e-9), max(abs(pub$dev_cf[hi]))))

ev <- readRDS(file.path("table_and_figure_manuscript",
                        "phi2_evaluations.rds"))$data
record("(3.3-6) Phi2 evaluations of the conventional search: median 31.5, range 2 to 125",
       nrow(ev) == 96L && abs(stats::median(ev$n_phi2_conv) - 31.5) < 1e-9 &&
         min(ev$n_phi2_conv) == 2L && max(ev$n_phi2_conv) == 125L,
       sprintf("cells = %d, median = %g, range = [%d, %d]", nrow(ev),
               stats::median(ev$n_phi2_conv), min(ev$n_phi2_conv),
               max(ev$n_phi2_conv)))

ev_kappa <- mapply(kappa_of, ev$delta1, ev$delta2, ev$sd1, ev$sd2)
top <- ev$n_phi2_conv == max(ev$n_phi2_conv)
record("(3.3-7) Largest Phi2 counts occur at kappa = 1.5",
       all(abs(ev_kappa[top] - 1.5) < 1e-9) &&
         all(abs(ev_kappa[order(-ev$n_phi2_conv)[1:5]] - 1.5) < 1e-9),
       sprintf("kappa of the five largest counts: %s",
               paste(sprintf("%.2f", ev_kappa[order(-ev$n_phi2_conv)[1:5]]),
                     collapse = ", ")))

record("(3.3-8) The closed form evaluates Phi2 at no point",
       all(ev$n_phi2_prop == 0L),
       sprintf("max count = %d", max(ev$n_phi2_prop)))

record("(3.3-9) Tables: 2 powers, 15 ratios in [1.00, 2.00], 7 correlations",
       setequal(unique(pub$power), c(0.80, 0.90)) &&
         length(unique(pub$gamma1)) == 15L &&
         abs(min(pub$gamma1) - 1) < 1e-9 && abs(max(pub$gamma1) - 2) < 1e-9 &&
         setequal(unique(pub$rho), c(0, 0.2, 0.3, 0.5, 0.7, 0.8, 0.95)),
       sprintf("powers = %s; ratios = %d; correlations = %s",
               paste(unique(pub$power), collapse = ", "),
               length(unique(pub$gamma1)),
               paste(unique(pub$rho), collapse = ", ")))

k_trial <- (3.1 / 12) / (1.8 / 9)
record("(3.3-10) EXPEDITION 1 (0.87, 1.292) is not covered by the tables",
       all(abs(unique(pub$power) - 0.87) > 1e-9) &&
         all(abs(unique(pub$gamma1) - k_trial) > 1e-3),
       sprintf("kappa = %.4f between tabulated %.2f and %.2f", k_trial,
               max(pub$gamma1[pub$gamma1 < k_trial]),
               min(pub$gamma1[pub$gamma1 > k_trial])))

# (3.3-11) For rho >= 0 the starting value (6) already attains the target
#          power, so the search of Algorithm 1 proceeds downward
start_ok <- sapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  z_b0 <- qnorm(sqrt(g$power))
  n2_0 <- ceiling(max(g$sd1 ^ 2 * (1 + 1 / g$r) / g$delta1 ^ 2,
                      g$sd2 ^ 2 * (1 + 1 / g$r) / g$delta2 ^ 2) *
                    (z_alpha + z_b0) ^ 2)
  twocpcont_power(n1 = ceiling(g$r * n2_0), n2 = n2_0,
                  delta1 = g$delta1, delta2 = g$delta2,
                  sd1 = g$sd1, sd2 = g$sd2, rho = g$rho,
                  alpha = alpha)$powerCoprimary >= g$power
})
record("(3.3-11) Starting value attains the target in all 96 cells (downward search)",
       length(start_ok) == 96L && all(start_ok),
       sprintf("cells = %d, downward = %d", length(start_ok), sum(start_ok)))

# (3.3-12) The discriminant of eq. (13) is nonnegative in every design of
#          Sections 3 and 4: the 96 cells of Table 2, the 210 tabulated
#          designs of Sozu et al. (2015), the four EXPEDITION 1 designs,
#          and the rho = 0 designs behind the closed-form R_max (kappa in
#          [1, 3] by 0.001, 1-beta in {0.70, 0.80, 0.85, 0.87, 0.90}).
disc_fun <- function(k, rr, tp) {
  p     <- sqrt(tp)
  u_p   <- qnorm(p)
  w     <- 0.5 + (k - 1) * (u_p + z_alpha) * dnorm(u_p) /
    ((1 + k) * abs(2 * log(p)) * p)
  w     <- min(1, max(0, w))
  lam0  <- z_alpha + qnorm(tp ^ w)
  u_w   <- lam0 - z_alpha
  u_s   <- k * lam0 - z_alpha
  pk    <- plackett_gl_full(u_w, u_s, rr, k, gl_nodes = 5L)
  P_at  <- pnorm(u_w) * pnorm(u_s) + pk$I_GL
  P_1   <- dnorm(u_w) * pnorm(u_s) + k * dnorm(u_s) * pnorm(u_w) + pk$I_GL_1
  P_2   <- -u_w * dnorm(u_w) * pnorm(u_s) +
    2 * k * dnorm(u_w) * dnorm(u_s) -
    k ^ 2 * u_s * dnorm(u_s) * pnorm(u_w) + pk$I_GL_2
  P_1 ^ 2 - 2 * P_2 * (P_at - tp)
}
disc_grid <- mapply(function(d1, d2, s1, s2, rr, tp) {
  disc_fun(kappa_of(d1, d2, s1, s2), rr, tp)
}, grid$delta1, grid$delta2, grid$sd1, grid$sd2, grid$rho, grid$power)
disc_pub <- mapply(disc_fun, pub$gamma1, pub$rho, pub$power)
disc_trial <- sapply(c(0, 0.3, 0.5, 0.8), function(rr) {
  disc_fun(k_trial, rr, 0.87)
})
disc_rho0 <- unlist(lapply(c(0.70, 0.80, 0.85, 0.87, 0.90), function(tp) {
  sapply(seq(1, 3, by = 0.001), disc_fun, rr = 0, tp = tp)
}))
disc_all <- c(disc_grid, disc_pub, disc_trial, disc_rho0)
record("(3.3-12) Discriminant of eq. (13) is nonnegative in all designs of Sections 3 and 4",
       length(disc_grid) == 96L && length(disc_pub) == 210L &&
         all(disc_all >= 0),
       sprintf("designs = %d, min discriminant = %.4f", length(disc_all),
               min(disc_all)))

# ============================================================
# SECTION 3.4: Maximum reduction rate and threshold
# ============================================================
section("Section 3.4 (Maximum reduction rate and threshold effect size ratio)")

Rmax_1 <- sapply(1 - power_set, function(b) {
  r_max(kappa = 1, alpha = alpha, beta = b) * 100
})
record("(3.4-1) Rmax(1) = 23.8%, 21.3%, 18.6% at 0.80, 0.85, 0.90",
       all(abs(Rmax_1 - c(23.8, 21.3, 18.6)) < 0.05),
       sprintf("%s", paste(sprintf("%.4f%%", Rmax_1), collapse = ", ")))

inc_cs <- (1 / (1 - r_max(kappa = 1, alpha = alpha, beta = 0.10)) - 1) * 100
record("(3.4-2) 1 / (1 - Rmax(1, 0.10)) - 1 approximately 22.8%",
       abs(inc_cs - 22.8) < 0.05,
       sprintf("increase = %.4f%%", inc_cs))

ks <- sapply(c(0.05, 0.15), function(e) {
  sapply(c(0.20, 0.10), function(b) kappa_star(e, alpha = alpha, beta = b))
})
record(paste0("(3.4-3) kappa_star: 1.37 -> 1.13 at 0.80, ",
              "1.24 -> 1.05 at 0.90"),
       all(abs(ks - rbind(c(1.37, 1.13), c(1.24, 1.05))) < 0.005),
       sprintf("0.80: %.4f -> %.4f; 0.90: %.4f -> %.4f",
               ks[1, 1], ks[1, 2], ks[2, 1], ks[2, 2]))

ks_05_max <- max(sapply(1 - power_set, function(b) {
  kappa_star(0.05, alpha = alpha, beta = b)
}))
record("(3.4-4) kappa_star(0.05) <= 1.4 at all three target powers",
       ks_05_max <= 1.4,
       sprintf("max kappa_star(0.05) = %.4f", ks_05_max))

record("(3.4-5) Rmax(1, 0.20) about 24% (upper bound at 1-beta = 0.80)",
       abs(Rmax_1[1] - 24) < 0.5,
       sprintf("Rmax = %.4f%%", Rmax_1[1]))

hw <- sapply(c(0.20, 0.15, 0.10), function(b) {
  1 / (1 - r_max(kappa = 1, alpha = alpha, beta = b))
})
hw_formula <- sapply(c(0.20, 0.15, 0.10), function(b) {
  ((z_alpha + qnorm(sqrt(1 - b))) / (z_alpha + qnorm(1 - b))) ^ 2
})
record("(3.4-6) Hung and Wang ratio 1 / (1 - Rmax(1)) = 1.31, 1.27, 1.23",
       all(abs(hw - c(1.31, 1.27, 1.23)) < 0.005) &&
         all(abs(hw - hw_formula) < 1e-10),
       sprintf("%s", paste(sprintf("%.4f", hw), collapse = ", ")))

Rmax_varga <- r_max(kappa = 0.5 / 0.4, alpha = alpha, beta = 0.20) * 100
record("(3.4-7) Varga et al. example: Rmax(1.25, 0.20) approximately 9.1%",
       abs(Rmax_varga - 9.1) < 0.05,
       sprintf("Rmax = %.4f%%", Rmax_varga))

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
record("(4-3) Doody example: kappa = 1.292",
       abs(kappa_doody - 1.292) < 0.0005,
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

# (4-13) On the continuous scale of eq. (14), the reduction at rho = 0.8
#        is about 3.9%; the gap to R_max is due to rho = 0.8 < 1, and the
#        ceiling operation increases the realized (integer) reduction
lam_doody <- function(rr) {
  uniroot(function(lam) {
    pbivnorm::pbivnorm(lam - z_alpha, kappa_doody * lam - z_alpha,
                       rho = rr) - 0.87
  }, lower = 2, upper = 6, tol = 1e-14)$root
}
red_cont <- (1 - lam_doody(0.8) ^ 2 / lam_doody(0) ^ 2) * 100
record("(4-13) Continuous reduction at rho = 0.8 about 3.9%, below R_max and 4.0%",
       abs(round(red_cont, 1) - 3.9) < 1e-9 && red_cont < Rmax_doody &&
         red_cont < red_doody,
       sprintf("continuous = %.4f%%, integer = %.2f%%, R_max = %.4f%%",
               red_cont, red_doody, Rmax_doody))

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
# SECTION 5: Discussion
# ============================================================
section("Section 5 (Discussion)")

# (5-1) Schouten (1999): increment z_alpha^2 / 2 is about two patients
incr_t <- qnorm(1 - 0.025) ^ 2 / 2
record("(5-1) z_alpha^2 / 2 at alpha = 0.025 is about two patients",
       round(incr_t) == 2,
       sprintf("z_alpha^2 / 2 = %.4f", incr_t))

# (5-2) One power check with n2 increased by one gives the conventional
#       sample size in the single shortfall of Section 3.2 (r = 1)
n2_fix <- mis$n2 + 1
N_plus <- ceiling(1 * n2_fix) + n2_fix
pw_fix <- twocpcont_power(n1 = ceiling(1 * n2_fix), n2 = n2_fix,
                          delta1 = 0.3, delta2 = 0.3, sd1 = 1, sd2 = 1,
                          rho = 0.5, alpha = alpha)$powerCoprimary
record("(5-2) Shortfall cell: n2 + 1 meets the target and gives N_conv = 558",
       pw_fix >= 0.90 && N_plus == mis$N_conv,
       sprintf("n2 = %d, N = %d, power = %.7f", as.integer(n2_fix),
               as.integer(N_plus), pw_fix))

# ============================================================
# APPENDIX E: Software implementation (claims in the prose around
# the transcript written by appendix_software_example.R)
# ============================================================
section("Appendix E (Software implementation)")

ss_app <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                       rho = 0.5, r = 1, alpha = 0.025, beta = 0.13,
                       method = "univariate")
pw_app <- twocpcont_power(n1 = ss_app$n1, n2 = ss_app$n2,
                          delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,
                          rho = 0.5, alpha = 0.025)$powerCoprimary
record("(E-1) rho = 0.5: N = 980 (Table 4) and co-primary power >= 0.87",
       ss_app$N == 980 && pw_app >= 0.87,
       sprintf("N = %d, power = %.6f", ss_app$N, pw_app))

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
