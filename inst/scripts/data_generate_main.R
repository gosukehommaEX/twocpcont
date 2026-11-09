# =============================================================================
# data_generate_main.R
#
# Computation stage for the manuscript figures and tables of the twocpcont
# article (a noniterative closed-form sample size formula for clinical
# trials with two co-primary continuous endpoints).
#
# This script performs all numerical computation and writes the results to
# inst/scripts/data/. The companion script table_and_figure_manuscript.R
# reads these files and produces the figures, the LaTeX tables and the
# ledger of the numbers cited in the text without any computation.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file:
#     devtools::load_all()
#     source("inst/scripts/data_generate_main.R")
#   Parts of the computation can be skipped with the logical vector RUN below;
#   the data files of the skipped parts are left unchanged.
#
# Required packages: twocpcont, pbivnorm
#
# Running time (rough guide; the time of each step is written to the log)
#   each part : a few seconds or less
#
# The computations are deterministic: no random numbers are used.
#
# Outputs (relative to the package root)
#   inst/scripts/data/tab1_gl_nodes.rds          Table 1 (Gauss-Legendre rule)
#   inst/scripts/data/fig1_anchor_weight.rds     Figure 1 (anchor weight w*)
#                                                and the exact root of (S2)
#   inst/scripts/data/tab2_sample_size.rds       Table 2 (96 design cells)
#   inst/scripts/data/phi2_evaluations.rds       Number of bivariate normal
#                                                CDF evaluations per cell
#   inst/scripts/data/fig2_reduction_curve.rds   Figure 2 (R_max curve)
#   inst/scripts/data/tab3_kappa_star.rds        Table 3 (kappa_star)
#   inst/scripts/data/tab4_real_example.rds      Table 4 (EXPEDITION 1)
#   inst/scripts/data/published_tables.rds       Sozu et al. (2015) and
#                                                Hung and Wang (2009)
#   inst/scripts/data/fig3_reduction_rho.rds     Figure 3 (reduction of the
#                                                sample size over rho)
#   inst/scripts/data/csv/*.csv        Data frames of the .rds files as CSV
#   inst/scripts/data/data_generate_main.log
#   inst/scripts/data/sessionInfo_data_generate_main.txt
# =============================================================================

# --------------------------------------------------------------------------- #
# Working directory and packages
# --------------------------------------------------------------------------- #
if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "twocpcont")) {
  stop("Run this script from the root of the twocpcont package ",
       "(open twocpcont.Rproj first).", call. = FALSE)
}
# A version loaded by devtools::load_all() is used if present; otherwise the
# installed package is loaded.
if (!requireNamespace("twocpcont", quietly = TRUE)) {
  stop("The twocpcont package is neither loaded nor installed.", call. = FALSE)
}
if (!"package:twocpcont" %in% search()) {
  library(twocpcont)
}
library(pbivnorm)

# Parts of the computation to run. RUN can also be assigned before this
# file is sourced to rerun only some parts, e.g.
#   RUN <- c(tab1 = FALSE, fig1 = TRUE, tab2 = FALSE, phi2 = FALSE,
#            fig2 = FALSE, tab3 = FALSE, tab4 = FALSE, pub = FALSE,
#            fig3 = FALSE)
RUN_PARTS <- c("tab1", "fig1", "tab2", "phi2", "fig2", "tab3", "tab4", "pub",
               "fig3")
if (!exists("RUN", inherits = FALSE)) {
  RUN <- c(tab1 = TRUE,   # Table 1
           fig1 = TRUE,   # Figure 1
           tab2 = TRUE,   # Table 2
           phi2 = TRUE,   # bivariate normal CDF evaluations (Section 3.3)
           fig2 = TRUE,   # Figure 2
           tab3 = TRUE,   # Table 3
           tab4 = TRUE,   # Table 4
           pub  = TRUE,   # published tables (Section 3.3)
           fig3 = TRUE)   # Figure 3
}
if (!is.logical(RUN) || !setequal(names(RUN), RUN_PARTS)) {
  stop("RUN must be a logical vector with the elements ",
       paste(RUN_PARTS, collapse = ", "), ".", call. = FALSE)
}

# --------------------------------------------------------------------------- #
# Output location and logging
# --------------------------------------------------------------------------- #
data_dir <- file.path("inst", "scripts", "data")
csv_dir  <- file.path(data_dir, "csv")
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_main.log")
cat(sprintf("data_generate_main.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}

# Evaluate an expression, log its running time, and write its warnings to
# the log instead of the console
run_step <- function(label, expr) {
  t0  <- proc.time()[["elapsed"]]
  val <- withCallingHandlers(expr, warning = function(w) {
    log_msg("  [warning] %s: %s", label, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  log_msg("%-44s %8.1f sec", label, proc.time()[["elapsed"]] - t0)
  val
}

# Save a list as .rds and each of its data frame elements as .csv
save_data <- function(obj, name) {
  obj$settings <- SETTINGS
  obj$meta     <- META
  saveRDS(obj, file.path(data_dir, paste0(name, ".rds")))
  for (el in names(obj)) {
    if (is.data.frame(obj[[el]])) {
      utils::write.csv(obj[[el]],
                       file.path(csv_dir, paste0(name, "_", el, ".csv")),
                       row.names = FALSE)
    }
  }
  log_msg("Saved %s", file.path(data_dir, paste0(name, ".rds")))
}

# Short commit hash and modified tracked files of the package repository
git_info <- function() {
  run_git <- function(args) {
    out <- tryCatch(suppressWarnings(system2("git", args, stdout = TRUE,
                                             stderr = FALSE)),
                    error = function(e) NULL)
    if (is.null(out) || !is.null(attr(out, "status"))) NULL else out
  }
  hash  <- run_git(c("rev-parse", "--short", "HEAD"))
  dirty <- run_git(c("status", "--porcelain", "--untracked-files=no"))
  list(commit   = if (length(hash) == 1L) hash else NA_character_,
       modified = if (is.null(dirty)) character(0) else trimws(dirty))
}

GIT  <- git_info()
META <- list(date              = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
             twocpcont_version = as.character(getNamespaceVersion("twocpcont")),
             twocpcont_path    = getNamespaceInfo("twocpcont", "path"),
             git_commit        = GIT$commit,
             git_modified      = GIT$modified,
             R_version         = R.version.string)

log_msg("twocpcont %s loaded from %s", META$twocpcont_version,
        META$twocpcont_path)
log_msg("git commit %s; modified tracked files: %s", META$git_commit,
        if (length(META$git_modified) == 0L) "none" else
          paste(META$git_modified, collapse = "; "))

# =============================================================================
# Common settings
# =============================================================================
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

# Real trial example based on Doody et al. (2014), EXPEDITION 1 of
# solanezumab (NEJM 2014; 370: 311-321), Protocol H8A-MC-LZAM Amendment (b),
# Section 12.1.1
trial_delta1 <- 1.8   # ADAS-Cog11 treatment difference at 18 months
trial_delta2 <- 3.1   # ADCS-ADL  treatment difference at 18 months
trial_sd1    <- 9     # standard deviation of ADAS-Cog11 change
trial_sd2    <- 12    # standard deviation of ADCS-ADL  change
trial_alpha  <- 0.025 # one-sided significance level
trial_beta   <- 0.13  # overall power 0.87 = 0.89 * 0.98
trial_r      <- 1     # allocation ratio (1:1)
trial_rhos   <- c(0.0, 0.3, 0.5, 0.8)
trial_es2_protocol <- 0.25   # effect size of ADCS-ADL stated in the protocol

SETTINGS <- list(
  alpha       = alpha,
  power_set   = power_set,
  gl_nodes    = gl_nodes,
  scenarios   = scenarios,
  rho_vals    = rho_vals,
  kappa_fig1  = seq(1, 4, length.out = 301),
  kappa_gap   = seq(1.001, 3, by = 0.001),
  kappa_fig2  = seq(1.0, 2.0, length.out = 201),
  epsilon_set = c(0.05, 0.10, 0.15),
  kappa_fig3  = c(1, 1.1, 1.2, 1.3, 1.5, 2),
  rho_fig3    = seq(0, 0.99, by = 0.01),
  power_fig3  = c(0.80, 0.90),
  varga       = c(es_w = 0.4, es_s = 0.5, power = 0.80),
  trial       = list(delta1 = trial_delta1, delta2 = trial_delta2,
                     sd1 = trial_sd1, sd2 = trial_sd2, alpha = trial_alpha,
                     beta = trial_beta, r = trial_r, rhos = trial_rhos,
                     es2_protocol = trial_es2_protocol)
)

z_alpha <- qnorm(1 - alpha)

# Anchor weight w* of eq. (eq:w_star), before the clipping to [0, 1]
w_star_unclipped <- function(kappa, beta, alpha) {
  z_alpha <- qnorm(1 - alpha)
  p       <- sqrt(1 - beta)
  u_p     <- qnorm(p)
  phi_u_p <- dnorm(u_p)
  L_abs   <- abs(2 * log(p))
  0.5 + (kappa - 1) * (u_p + z_alpha) * phi_u_p / ((1 + kappa) * L_abs * p)
}

# Exact noncentrality parameter of the weaker endpoint: the root in lambda
# of Phi2(lambda - z_alpha, kappa lambda - z_alpha; rho) = 1 - beta, for
# rho >= 0. The root lies between z_alpha + z_beta and z_alpha + u_p.
lambda_exact <- function(kappa, rho, beta, alpha) {
  z_alpha <- qnorm(1 - alpha)
  uniroot(function(lam) {
    pbivnorm(lam - z_alpha, kappa * lam - z_alpha, rho = rho) - (1 - beta)
  }, lower = z_alpha + qnorm(1 - beta) - 0.5,
  upper = z_alpha + qnorm(sqrt(1 - beta)) + 0.5, tol = 1e-14)$root
}

# =============================================================================
# Table 1: Gauss-Legendre nodes and weights
# =============================================================================
if (RUN[["tab1"]]) {
  gl_obj <- run_step("Table 1 (Gauss-Legendre rule)", {
    gl <- GL_nodes_and_weights(gl_nodes = gl_nodes)
    list(x = gl$x, w = gl$w, gl_nodes = gl_nodes)
  })
  save_data(gl_obj, "tab1_gl_nodes")
}

# =============================================================================
# Figure 1: anchor weight w* as a function of kappa, eq. (eq:w_star):
#   w_star = 1/2 + (kappa - 1) (u_p + z_alpha) phi(u_p) /
#                  [(1 + kappa) |2 log p| p]
# together with the exact root of F(w; kappa) = 0 in (S2) of the
# Supporting Information, and the kappa at which w* reaches one
# =============================================================================
if (RUN[["fig1"]]) {
  anchor_obj <- run_step("Figure 1 (anchor weight)", {
    # Exact root of F(w; kappa) = 0, eq. (S2). F is increasing in w and
    # tends to +Inf as w -> 1, but in double precision F(1 - 1e-12) can
    # still be negative for large kappa; the root is then within 1e-12 of
    # one.
    w_root_fun <- function(k, tp) {
      if (k == 1) return(0.5)
      f_w <- function(w) {
        qnorm(tp ^ (1 - w)) - k * qnorm(tp ^ w) - (k - 1) * z_alpha
      }
      if (f_w(1 - 1e-12) <= 0) return(1)
      uniroot(f_w, lower = 0.5, upper = 1 - 1e-12, tol = 1e-13)$root
    }

    # w* (clipped to [0, 1]) and the exact root on the grid of Figure 1
    kappa_seq_w <- SETTINGS$kappa_fig1
    df_w <- do.call(rbind, lapply(power_set, function(pw) {
      data.frame(
        kappa  = kappa_seq_w,
        power  = pw,
        w_star = pmin(1, pmax(0, w_star_unclipped(kappa_seq_w, beta = 1 - pw,
                                                  alpha = alpha))),
        w_root = sapply(kappa_seq_w, w_root_fun, tp = pw)
      )
    }))
    df_gap <- do.call(rbind, lapply(power_set, function(pw) {
      k_seq  <- SETTINGS$kappa_gap
      w_cf   <- pmin(1, w_star_unclipped(k_seq, beta = 1 - pw, alpha = alpha))
      w_ex   <- sapply(k_seq, w_root_fun, tp = pw)
      data.frame(kappa = k_seq, power = pw, w_star = w_cf, w_root = w_ex,
                 gap = w_cf - w_ex)
    }))

    # kappa at which the unclipped w* equals one
    df_one <- data.frame(
      power = power_set,
      kappa = sapply(power_set, function(pw) {
        uniroot(function(k) w_star_unclipped(k, beta = 1 - pw,
                                             alpha = alpha) - 1,
                lower = 1, upper = 3, tol = 1e-10)$root
      })
    )

    list(data = df_w, gap = df_gap, w_one = df_one, alpha = alpha,
         power_set = power_set, kappa_seq = kappa_seq_w)
  })
  save_data(anchor_obj, "fig1_anchor_weight")
}

# =============================================================================
# Table 2: sample size comparison over the 96 design cells
# =============================================================================
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

  # Exact continuous n2: c (lambda_w)^2 at the exact root of the power
  # equation (eq:power_biv)
  g1 <- sc$delta1 / sc$sd1
  g2 <- sc$delta2 / sc$sd2
  gamma_w <- min(g1, g2)
  kappa   <- max(g1, g2) / gamma_w
  c_scale <- (1 + 1 / sc$r) / gamma_w ^ 2
  n2_exact <- c_scale * lambda_exact(kappa, rho, beta, alpha) ^ 2

  list(
    delta1   = sc$delta1,
    delta2   = sc$delta2,
    sd1      = sc$sd1,
    sd2      = sc$sd2,
    r        = sc$r,
    rho      = rho,
    # Rounded so that the value matches power_set exactly
    power    = round(1 - beta, 10),
    kappa    = kappa,
    N_conv   = ss_conv$N,
    N_prop   = ss_prop$N,
    diff     = ss_prop$N - ss_conv$N,
    pow      = pwr_prop$powerCoprimary,
    n2_conv  = ss_conv$n2,
    n2_prop  = ss_prop$n2,
    n2_cont  = ss_prop$n2_cont,
    n2_exact = n2_exact
  )
}

if (RUN[["tab2"]]) {
  table2_obj <- run_step("Table 2 (96 cells)", {
    rows_known <- list()
    for (pw in power_set) {
      for (sc in scenarios) {
        for (rho in rho_vals) {
          rows_known[[length(rows_known) + 1L]] <- compute_cell(sc, rho,
                                                                1 - pw)
        }
      }
    }
    df_known <- do.call(rbind, lapply(rows_known, as.data.frame))
    list(data = df_known, alpha = alpha, power_set = power_set,
         gl_nodes = gl_nodes)
  })
  save_data(table2_obj, "tab2_sample_size")
}

# =============================================================================
# Number of bivariate normal CDF evaluations per Table 2 cell
#
# Every co-primary power evaluation of the conventional method calls
# pbivnorm::pbivnorm() exactly once.  The calls are counted with
# trace(), so that the count is taken from the implementation itself
# rather than from a re-implementation of the search.  The same count
# is taken for the closed form, which should make no call.  The count
# does not depend on the computer.
# =============================================================================
if (RUN[["phi2"]]) {
  phi2_obj <- run_step("Bivariate normal CDF evaluations", {
    phi2_counter <- new.env()
    phi2_counter$n <- 0L
    suppressMessages(trace("pbivnorm", where = asNamespace("pbivnorm"),
                           tracer = quote(phi2_counter$n <-
                                            phi2_counter$n + 1L),
                           print = FALSE))
    count_calls <- function(expr) {
      phi2_counter$n <- 0L
      force(expr)
      phi2_counter$n
    }

    rows_eval <- list()
    for (pw in power_set) {
      for (sc in scenarios) {
        for (rho in rho_vals) {
          n_conv <- count_calls(twocpcont_ss(
            delta1 = sc$delta1, delta2 = sc$delta2,
            sd1 = sc$sd1, sd2 = sc$sd2, rho = rho, r = sc$r,
            alpha = alpha, beta = 1 - pw, method = "bivariate"))
          n_prop <- count_calls(twocpcont_ss(
            delta1 = sc$delta1, delta2 = sc$delta2,
            sd1 = sc$sd1, sd2 = sc$sd2, rho = rho, r = sc$r,
            alpha = alpha, beta = 1 - pw, method = "univariate",
            gl_nodes = gl_nodes))
          rows_eval[[length(rows_eval) + 1L]] <- list(
            power  = pw,
            delta1 = sc$delta1, delta2 = sc$delta2,
            sd1    = sc$sd1,    sd2    = sc$sd2,
            r      = sc$r,      rho    = rho,
            kappa  = max(sc$delta1 / sc$sd1, sc$delta2 / sc$sd2) /
              min(sc$delta1 / sc$sd1, sc$delta2 / sc$sd2),
            n_phi2_conv = n_conv,
            n_phi2_prop = n_prop
          )
        }
      }
    }
    suppressMessages(untrace("pbivnorm", where = asNamespace("pbivnorm")))
    df_eval <- do.call(rbind, lapply(rows_eval, as.data.frame))
    list(data = df_eval, alpha = alpha, power_set = power_set,
         gl_nodes = gl_nodes)
  })
  save_data(phi2_obj, "phi2_evaluations")
}

# =============================================================================
# Figure 2: R_max(kappa) curve at three power levels, and the values of
# R_max cited in Section 3.4
# =============================================================================
if (RUN[["fig2"]]) {
  reduction_obj <- run_step("Figure 2 (reduction curve)", {
    kappa_seq <- SETTINGS$kappa_fig2
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

    # Example of Varga et al. (2017): standardized effect sizes 0.4 and 0.5
    # at 80% power
    vg <- SETTINGS$varga
    kappa_varga <- vg[["es_s"]] / vg[["es_w"]]
    df_varga <- data.frame(
      es_w  = vg[["es_w"]],
      es_s  = vg[["es_s"]],
      kappa = kappa_varga,
      power = vg[["power"]],
      R_max = r_max(kappa = kappa_varga, alpha = alpha,
                    beta = 1 - vg[["power"]])
    )

    list(data = df_rmax, varga = df_varga, alpha = alpha,
         power_set = power_set, kappa_seq = kappa_seq)
  })
  save_data(reduction_obj, "fig2_reduction_curve")
}

# =============================================================================
# Table 3: kappa_star(epsilon, beta) for selected tolerances
# =============================================================================
if (RUN[["tab3"]]) {
  kappa_star_obj <- run_step("Table 3 (kappa_star)", {
    epsilon_set <- SETTINGS$epsilon_set
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
    list(data = df_kstar, alpha = alpha, power_set = power_set,
         epsilon_set = epsilon_set)
  })
  save_data(kappa_star_obj, "tab3_kappa_star")
}

# =============================================================================
# Table 4: real trial example based on Doody et al. (2014), EXPEDITION 1.
#
# The protocol prescribed 500 patients per group with co-primary
# endpoints ADAS-Cog11 and ADCS-ADL, assuming independence between
# the two endpoints.  The planned sample size is re-derived with the
# closed-form formula at rho = 0, and rho is varied to quantify the
# design-time benefit of accounting for the between-endpoint correlation.
# =============================================================================
if (RUN[["tab4"]]) {
  table4_obj <- run_step("Table 4 (EXPEDITION 1)", {
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

    # Exact noncentrality parameter at each rho; the reduction on the
    # continuous scale of eq. (eq:R_def) follows from its square
    df_lam <- data.frame(
      rho    = trial_rhos,
      lambda = sapply(trial_rhos, function(rr) {
        lambda_exact(kappa_trial, rr, trial_beta, trial_alpha)
      })
    )
    df_lam$reduction_cont_pct <-
      (1 - df_lam$lambda ^ 2 / df_lam$lambda[df_lam$rho == 0] ^ 2) * 100

    # Effect size ratio with the rounded effect size 0.25 of the protocol
    kappa_round <- trial_es2_protocol / (trial_delta1 / trial_sd1)
    R_max_round <- r_max(kappa = kappa_round, alpha = trial_alpha,
                         beta = trial_beta) * 100

    list(
      data         = df_t4,
      lambda       = df_lam,
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
      kappa_round  = kappa_round,
      R_max_round_pct = R_max_round,
      gl_nodes     = gl_nodes
    )
  })
  save_data(table4_obj, "tab4_real_example")
}

# =============================================================================
# Comparison with published values (Section 3.3)
#
# (a) Sozu et al. (2015, Tables 4.3 and 4.4) tabulate the constant
#     C_2 = lambda_w - z_alpha for alpha = 0.025, 1 - beta = 0.80 and
#     0.90, 15 effect size ratios gamma1 and 7 correlations.  C_2 is
#     computed here by (i) the exact root of the bivariate normal power
#     equation and (ii) the closed form, through n2_cont of
#     twocpcont_ss() with delta_w = sd_w = 1 and r = 1, so that c = 2
#     and lambda_w* = sqrt(n2_cont / 2).
# (b) Hung and Wang (2009, Section 3) give the ratio n / m1 of the
#     two-endpoint sample size under equal effect sizes and
#     independence to the single-endpoint sample size.  For K = 2 it
#     equals 1 / {1 - r_max(1, beta)}.
# The value-by-value reproduction of the published values is
# inst/reproduce/reproduce_Sozu2015.R and reproduce_HungWang2009.R.
# =============================================================================
if (RUN[["pub"]]) {
  published_obj <- run_step("Published tables", {
    tab_file <- system.file("extdata", "sozu2015_C2_table.csv",
                            package = "twocpcont")
    if (!nzchar(tab_file) || !file.exists(tab_file)) {
      stop("sozu2015_C2_table.csv not found.", call. = FALSE)
    }
    sozu_tab <- utils::read.csv(tab_file)

    C2_exact <- mapply(function(g, rh, pw) {
      uniroot(function(l) {
        pbivnorm(l - z_alpha, g * l - z_alpha, rho = rh) - pw
      }, lower = z_alpha, upper = z_alpha + 5, tol = 1e-12)$root - z_alpha
    }, sozu_tab$gamma1, sozu_tab$rho, sozu_tab$power)
    C2_cf <- mapply(function(g, rh, pw) {
      ss <- twocpcont_ss(delta1 = g, delta2 = 1, sd1 = 1, sd2 = 1,
                         rho = rh, r = 1, alpha = alpha, beta = 1 - pw,
                         method = "univariate", gl_nodes = gl_nodes)
      sqrt(ss$n2_cont / 2) - z_alpha
    }, sozu_tab$gamma1, sozu_tab$rho, sozu_tab$power)

    df_pub <- data.frame(
      table     = sozu_tab$table,
      power     = sozu_tab$power,
      gamma1    = sozu_tab$gamma1,
      rho       = sozu_tab$rho,
      C2_pub    = sozu_tab$C2,
      C2_exact  = C2_exact,
      C2_cf     = C2_cf,
      dev_exact = round(round(C2_exact, 3) - sozu_tab$C2, 3),
      dev_cf    = round(round(C2_cf, 3) - sozu_tab$C2, 3),
      err_cf    = C2_cf - C2_exact
    )

    df_hw <- data.frame(power = power_set)
    df_hw$ratio_hw <- ((z_alpha + qnorm(sqrt(df_hw$power))) /
                         (z_alpha + qnorm(df_hw$power))) ^ 2
    df_hw$ratio_rmax <- sapply(df_hw$power, function(pw) {
      1 / (1 - r_max(kappa = 1, alpha = alpha, beta = 1 - pw))
    })

    list(sozu = df_pub, hung_wang = df_hw, alpha = alpha,
         gl_nodes = gl_nodes)
  })
  save_data(published_obj, "published_tables")
}

# =============================================================================
# Figure 3: reduction of the sample size on the continuous scale,
# R(rho; kappa, beta) = 1 - {lambda_w(rho) / lambda_w(0)}^2 of
# eq. (eq:R_def), as a function of rho, with its limit R_max(kappa, beta)
# as rho -> 1 (eq. (eq:Rmax)). lambda_w(rho) is the exact root of the power
# equation (eq:power_biv).
# =============================================================================
if (RUN[["fig3"]]) {
  reduction_rho_obj <- run_step("Figure 3 (reduction over rho)", {
    grid3 <- expand.grid(rho = SETTINGS$rho_fig3, kappa = SETTINGS$kappa_fig3,
                         power = SETTINGS$power_fig3)
    grid3$lambda <- mapply(function(k, rr, pw) {
      lambda_exact(k, rr, 1 - pw, alpha)
    }, grid3$kappa, grid3$rho, grid3$power)
    key  <- paste(grid3$kappa, grid3$power)
    lam0 <- grid3$lambda[grid3$rho == 0]
    names(lam0) <- key[grid3$rho == 0]
    grid3$reduction <- 1 - grid3$lambda ^ 2 / lam0[key] ^ 2
    limit <- expand.grid(kappa = SETTINGS$kappa_fig3,
                         power = SETTINGS$power_fig3)
    limit$R_max <- mapply(function(k, pw) {
      r_max(kappa = k, alpha = alpha, beta = 1 - pw)
    }, limit$kappa, limit$power)
    list(data = grid3, limit = limit, alpha = alpha)
  })
  save_data(reduction_rho_obj, "fig3_reduction_rho")
}

# =============================================================================
# Record the execution environment for reproducibility
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_main.txt")
)

log_msg("All requested data written to '%s/'.", data_dir)

# A RUN assigned for this run must not carry over to the next one
rm(RUN)
