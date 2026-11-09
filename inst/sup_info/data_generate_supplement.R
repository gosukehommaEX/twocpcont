# =============================================================================
# data_generate_supplement.R
#
# Computation stage for the Supporting Information of the twocpcont article.
# It writes the results to inst/sup_info/data/, which
# table_and_figure_supplement.R reads. The structure and the helpers are
# those of inst/scripts/data_generate_main.R.
#
# Parts (logical vector RUN_SI; all TRUE by default)
#   s5      Section S9: code and output of the software example
#   decomp  Error of the noncentrality parameter by approximation step:
#           anchor only, first-order correction, proposed quadratic
#           correction, exact root of the Gauss-Legendre power equation
#           (quadrature error alone) and quadratic correction with the
#           exact bivariate normal CDF (Taylor error alone)
#   grid    Accuracy of the closed form over a wide grid of
#           (alpha, 1 - beta, kappa, rho)
#   cost    Comparison of the closed form with the conventional search,
#           the Newton-Raphson algorithm of Sugimoto et al. (2012) and
#           Brent's method (uniroot) on the 96 combinations of Table 2:
#           evaluations of the bivariate normal CDF, iterations, time and
#           accuracy; convergence of the Newton-Raphson algorithm at
#           extreme correlations
#   alloc   Noninteger and unbalanced allocation ratios
#   sens    EXPEDITION 1: sensitivity to the effect size of ADCS-ADL and to
#           rho, and power when the correlation assumed at the planning
#           stage differs from the true correlation
#
# Usage
#   Open twocpcont.Rproj, load the package, and source this file:
#     devtools::load_all()
#     source("inst/sup_info/data_generate_supplement.R")
#   Parts can be selected by assigning RUN_SI before sourcing, e.g.
#     RUN_SI <- c(s5 = FALSE, decomp = TRUE, grid = FALSE, cost = FALSE,
#                 alloc = FALSE, sens = FALSE)
#
# Running time (rough guide; written to data/data_generate_supplement.log)
#   s5, decomp, alloc, sens : seconds; grid, cost : about 1 minute each
#
# The computations are deterministic: no random numbers are used. The
# computing times of the part cost depend on the computer.
#
# Outputs (relative to the package root)
#   inst/sup_info/data/software_example.rds      Section S9 (code and output)
#   inst/sup_info/data/error_decomposition.rds   part decomp
#   inst/sup_info/data/accuracy_grid.rds         part grid
#   inst/sup_info/data/cost_comparison.rds       part cost
#   inst/sup_info/data/allocation.rds            part alloc
#   inst/sup_info/data/expedition_sensitivity.rds  part sens
#   inst/sup_info/data/csv/*.csv                 Data frames as CSV
#   inst/sup_info/data/data_generate_supplement.log
#   inst/sup_info/data/sessionInfo_data_generate_supplement.txt
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

RUN_SI_PARTS <- c("s5", "decomp", "grid", "cost", "alloc", "sens")
if (!exists("RUN_SI", inherits = FALSE)) {
  RUN_SI <- c(s5 = TRUE, decomp = TRUE, grid = TRUE, cost = TRUE,
              alloc = TRUE, sens = TRUE)
}
if (!is.logical(RUN_SI) || !setequal(names(RUN_SI), RUN_SI_PARTS)) {
  stop("RUN_SI must be a logical vector with the elements ",
       paste(RUN_SI_PARTS, collapse = ", "), ".", call. = FALSE)
}

data_dir <- file.path("inst", "sup_info", "data")
csv_dir  <- file.path(data_dir, "csv")
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_supplement.log")
cat(sprintf("data_generate_supplement.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)
log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}
run_step <- function(label, expr) {
  t0  <- proc.time()[["elapsed"]]
  val <- withCallingHandlers(expr, warning = function(w) {
    log_msg("  [warning] %s: %s", label, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  log_msg("%-44s %8.1f sec", label, proc.time()[["elapsed"]] - t0)
  val
}

META <- list(date              = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
             twocpcont_version = as.character(getNamespaceVersion("twocpcont")),
             R_version         = R.version.string,
             platform          = R.version$platform)

# Save a list as .rds and each of its data frame elements as .csv
save_data <- function(obj, name) {
  obj$meta <- META
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

# =============================================================================
# Helpers
# =============================================================================
alpha <- 0.025

# Exact noncentrality parameter of the weaker endpoint for any rho in
# (-1, 1): the root in lambda of
# Phi2(lambda - z_alpha, kappa lambda - z_alpha; rho) = 1 - beta. Since
# 2 Phi(u_w) - 1 <= Phi2 <= Phi(u_w), the root lies in
# [z_alpha + z_beta, z_alpha + qnorm(1 - beta / 2)].
lambda_exact_any <- function(kappa, rho, beta, alpha, tol = 1e-14) {
  za <- qnorm(1 - alpha)
  uniroot(function(l) {
    pbivnorm::pbivnorm(l - za, kappa * l - za, rho = rho) - (1 - beta)
  }, lower = za + qnorm(1 - beta) - 1e-8,
  upper = za + qnorm(1 - beta / 2) + 1e-8, tol = tol)$root
}

# Power function P of (eq:power_approx) with the 5-point Gauss-Legendre
# rule, and its first two derivatives in lambda (Section S3)
p_gl <- function(lam, kappa, rho, za) {
  u_w <- lam - za
  u_s <- kappa * lam - za
  pk  <- plackett_gl_full(u_w, u_s, rho, kappa, gl_nodes = 5L)
  list(P  = pnorm(u_w) * pnorm(u_s) + pk$I_GL,
       P1 = dnorm(u_w) * pnorm(u_s) + kappa * dnorm(u_s) * pnorm(u_w) +
         pk$I_GL_1,
       P2 = -u_w * dnorm(u_w) * pnorm(u_s) +
         2 * kappa * dnorm(u_w) * dnorm(u_s) -
         kappa ^ 2 * u_s * dnorm(u_s) * pnorm(u_w) + pk$I_GL_2)
}

# Exact power Phi2(u_w, u_s; rho) and its first two derivatives in lambda:
# with s = sqrt(1 - rho^2), a = (u_s - rho u_w) / s, b = (u_w - rho u_s) / s
# and phi2 the bivariate normal density at (u_w, u_s),
#   P1 = phi(u_w) Phi(a) + kappa phi(u_s) Phi(b),
#   P2 = -u_w phi(u_w) Phi(a) - kappa^2 u_s phi(u_s) Phi(b)
#        + (2 kappa - rho - rho kappa^2) phi2
p_exact <- function(lam, kappa, rho, za) {
  u_w  <- lam - za
  u_s  <- kappa * lam - za
  s    <- sqrt(1 - rho ^ 2)
  a    <- (u_s - rho * u_w) / s
  b    <- (u_w - rho * u_s) / s
  phi2 <- dnorm(u_w) * dnorm(a) / s
  list(P  = pbivnorm::pbivnorm(u_w, u_s, rho = rho),
       P1 = dnorm(u_w) * pnorm(a) + kappa * dnorm(u_s) * pnorm(b),
       P2 = -u_w * dnorm(u_w) * pnorm(a) -
         kappa ^ 2 * u_s * dnorm(u_s) * pnorm(b) +
         (2 * kappa - rho - rho * kappa ^ 2) * phi2)
}

# Anchor lambda_w^(0) of (eq:lambda_anchor) with w* of (eq:w_star)
anchor_lambda <- function(kappa, beta, alpha) {
  za  <- qnorm(1 - alpha)
  tp  <- 1 - beta
  p   <- sqrt(tp)
  u_p <- qnorm(p)
  w   <- 0.5 + (kappa - 1) * (u_p + za) * dnorm(u_p) /
    ((1 + kappa) * abs(2 * log(p)) * p)
  za + qnorm(tp ^ min(1, max(0, w)))
}

# Root with the smaller absolute value of C0 + P1 d + P2 d^2 / 2 = 0
quad_step <- function(C0, P1, P2) -2 * C0 / (P1 + sqrt(P1 ^ 2 - 2 * P2 * C0))

# lambda_w* of the package: with delta_w = sd_w = 1 and r = 1, c = 2
lambda_package <- function(kappa, rho, beta, alpha) {
  sqrt(twocpcont_ss(delta1 = kappa, delta2 = 1, sd1 = 1, sd2 = 1, rho = rho,
                    r = 1, alpha = alpha, beta = beta,
                    method = "univariate")$n2_cont / 2)
}

# =============================================================================
# Section S9: software example for the EXPEDITION 1 design with rho = 0.5
# =============================================================================
if (RUN_SI[["s5"]]) {
  # Code shown in Section S9, one element per console input. The first
  # element is displayed but not evaluated, because the package has already
  # been loaded above.
  chunks <- c(
    "library(twocpcont)",
    paste0("ss <- twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,\n",
           "                   rho = 0.5, r = 1, alpha = 0.025, beta = 0.13,\n",
           "                   method = \"univariate\")"),
    "ss",
    "ss$n2_cont",
    paste0("twocpcont_power(n1 = ss$n1, n2 = ss$n2, delta1 = 1.8, delta2 = 3.1,\n",
           "                sd1 = 9, sd2 = 12, rho = 0.5, alpha = 0.025)"),
    paste0("sapply(c(0, 0.3, 0.5, 0.8), function(rho) {\n",
           "  twocpcont_ss(delta1 = 1.8, delta2 = 3.1, sd1 = 9, sd2 = 12,\n",
           "               rho = rho, alpha = 0.025, beta = 0.13,\n",
           "               method = \"univariate\")$N\n",
           "})"),
    "kappa <- (3.1 / 12) / (1.8 / 9)",
    "r_max(kappa, alpha = 0.025, beta = 0.13)",
    "kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.13)"
  )

  # Evaluate each chunk and build the console transcript
  old_opts <- options(width = 72, digits = 7)
  ex_env <- new.env(parent = globalenv())
  transcript <- character(0)
  for (i in seq_along(chunks)) {
    code_lines <- strsplit(chunks[i], "\n", fixed = TRUE)[[1]]
    prompts <- c("> ", rep("+ ", length(code_lines) - 1L))
    transcript <- c(transcript, paste0(prompts, code_lines))
    if (i == 1L) next
    out <- utils::capture.output({
      res <- withVisible(eval(parse(text = chunks[i]), envir = ex_env))
      if (res$visible) print(res$value)
    })
    # Drop leading and trailing blank lines produced by the print methods
    while (length(out) > 0L && !nzchar(trimws(out[1L]))) out <- out[-1L]
    while (length(out) > 0L && !nzchar(trimws(out[length(out)]))) {
      out <- out[-length(out)]
    }
    transcript <- c(transcript, out)
  }
  options(old_opts)

  software_example <- list(chunks = chunks, transcript = transcript,
                           R_version = R.version.string,
                           twocpcont_version =
                             as.character(getNamespaceVersion("twocpcont")))
  saveRDS(software_example, file.path(data_dir, "software_example.rds"))
  message("Saved ", file.path(data_dir, "software_example.rds"), " (",
          length(transcript), " transcript lines)")
}

# =============================================================================
# Error of the noncentrality parameter by approximation step
# =============================================================================
if (RUN_SI[["decomp"]]) {
  decomp <- run_step("Error decomposition", {
    dg <- expand.grid(kappa = seq(1, 3, by = 0.01),
                      rho = c(0, 0.3, 0.5, 0.8, 0.95),
                      power = c(0.80, 0.90))
    za <- qnorm(1 - alpha)
    res <- t(mapply(function(k, rr, pw) {
      tp   <- pw
      lam0 <- anchor_lambda(k, 1 - pw, alpha)
      gl   <- p_gl(lam0, k, rr, za)
      ex   <- p_exact(lam0, k, rr, za)
      lam_glroot <- uniroot(function(l) p_gl(l, k, rr, za)$P - tp,
                            lower = za + qnorm(tp) - 0.5,
                            upper = za + qnorm(1 - (1 - tp) / 2) + 0.5,
                            tol = 1e-14)$root
      c(lam_anchor    = lam0,
        lam_first     = lam0 - (gl$P - tp) / gl$P1,
        lam_proposed  = lam0 + quad_step(gl$P - tp, gl$P1, gl$P2),
        lam_package   = lambda_package(k, rr, 1 - pw, alpha),
        lam_gl_root   = lam_glroot,
        lam_exact_quad = lam0 + quad_step(ex$P - tp, ex$P1, ex$P2),
        lam_exact     = lambda_exact_any(k, rr, 1 - pw, alpha),
        disc_gl       = gl$P1 ^ 2 - 2 * gl$P2 * (gl$P - tp),
        P_gl = gl$P, P1_gl = gl$P1, P2_gl = gl$P2,
        P_ex = ex$P, P1_ex = ex$P1, P2_ex = ex$P2)
    }, dg$kappa, dg$rho, dg$power))
    dg <- cbind(dg, res)
    for (v in c("anchor", "first", "proposed", "gl_root", "exact_quad")) {
      dg[[paste0("err_", v)]] <- dg[[paste0("lam_", v)]] - dg$lam_exact
      dg[[paste0("rel_n2_", v)]] <-
        (dg[[paste0("lam_", v)]] / dg$lam_exact) ^ 2 - 1
    }
    list(data = dg, alpha = alpha)
  })
  save_data(decomp, "error_decomposition")
}

# =============================================================================
# Accuracy of the closed form over a wide grid of (alpha, 1 - beta, kappa,
# rho). The relative error of the continuous n2, (lambda_w* / lambda_w)^2 - 1,
# does not depend on delta, sigma or r; the absolute error is this relative
# error times n2.
# =============================================================================
if (RUN_SI[["grid"]]) {
  acc <- run_step("Accuracy over a wide grid", {
    ag <- expand.grid(kappa = seq(1, 3, by = 0.02),
                      rho = c(-0.5, -0.3, 0, 0.3, 0.5, 0.8, 0.9, 0.95, 0.99),
                      power = round(seq(0.70, 0.95, by = 0.05), 2),
                      alpha = c(0.005, 0.01, 0.025, 0.05))
    res <- t(mapply(function(k, rr, pw, a) {
      za     <- qnorm(1 - a)
      lam_cf <- lambda_package(k, rr, 1 - pw, a)
      lam_ex <- lambda_exact_any(k, rr, 1 - pw, a)
      lam0   <- anchor_lambda(k, 1 - pw, a)
      gl     <- p_gl(lam0, k, rr, za)
      c(lam_cf = lam_cf, lam_exact = lam_ex,
        power_cf = pbivnorm::pbivnorm(lam_cf - za, k * lam_cf - za, rho = rr),
        disc = gl$P1 ^ 2 - 2 * gl$P2 * (gl$P - pw))
    }, ag$kappa, ag$rho, ag$power, ag$alpha))
    ag <- cbind(ag, res)
    ag$err_lambda <- ag$lam_cf - ag$lam_exact
    ag$rel_n2     <- (ag$lam_cf / ag$lam_exact) ^ 2 - 1
    ag$power_diff <- ag$power_cf - ag$power
    list(data = ag)
  })
  save_data(acc, "accuracy_grid")
}

# =============================================================================
# Comparison of the methods on the 96 combinations of Table 2
# =============================================================================
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
cells <- do.call(rbind, lapply(c(0.80, 0.85, 0.90), function(pw) {
  do.call(rbind, lapply(scenarios, function(sc) {
    data.frame(delta1 = sc$delta1, delta2 = sc$delta2, sd1 = sc$sd1,
               sd2 = sc$sd2, r = sc$r, rho = c(0, 0.3, 0.5, 0.8), power = pw)
  }))
}))

# Brent's method (uniroot) on the exact power equation (eq:power_biv), with
# the output of twocpcont_ss() so that the computing times are comparable
UNIROOT_TOL <- 1e-8
ss_uniroot <- function(delta1, delta2, sd1, sd2, rho, r = 1, alpha = 0.025,
                       beta = 0.2) {
  g1 <- delta1 / sd1
  g2 <- delta2 / sd2
  gamma_w <- min(g1, g2)
  kappa   <- max(g1, g2) / gamma_w
  c_scale <- (1 + 1 / r) / gamma_w ^ 2
  za  <- qnorm(1 - alpha)
  lam <- uniroot(function(l) {
    pbivnorm::pbivnorm(l - za, kappa * l - za, rho = rho) - (1 - beta)
  }, lower = za + qnorm(1 - beta), upper = za + qnorm(1 - beta / 2),
  tol = UNIROOT_TOL)$root
  n2_cont <- c_scale * lam ^ 2
  n2 <- ceiling(n2_cont)
  n1 <- ceiling(r * n2)
  data.frame(delta1 = delta1, delta2 = delta2, sd1 = sd1, sd2 = sd2,
             rho = rho, r = r, alpha = alpha, beta = beta,
             method = "uniroot", n1 = n1, n2 = n2, N = n1 + n2,
             n2_cont = n2_cont)
}

METHODS <- c("univariate", "bivariate", "newton", "uniroot")
ss_method <- function(cell, m) {
  if (m == "uniroot") {
    ss_uniroot(cell$delta1, cell$delta2, cell$sd1, cell$sd2, cell$rho,
               cell$r, alpha, 1 - cell$power)
  } else {
    twocpcont_ss(delta1 = cell$delta1, delta2 = cell$delta2,
                 sd1 = cell$sd1, sd2 = cell$sd2, rho = cell$rho, r = cell$r,
                 alpha = alpha, beta = 1 - cell$power, method = m)
  }
}

if (RUN_SI[["cost"]]) {
  cost <- run_step("Comparison of the methods", {
    # Evaluations of the bivariate normal CDF, counted with trace() on
    # pbivnorm::pbivnorm(), which every method uses
    phi2_counter <- new.env()
    phi2_counter$n <- 0L
    suppressMessages(trace("pbivnorm", where = asNamespace("pbivnorm"),
                           tracer = quote(phi2_counter$n <-
                                            phi2_counter$n + 1L),
                           print = FALSE))
    per_cell <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
      cell <- cells[i, ]
      out <- lapply(METHODS, function(m) {
        phi2_counter$n <- 0L
        ss <- ss_method(cell, m)
        data.frame(method = m, N = ss$N, n2 = ss$n2, n2_cont = ss$n2_cont,
                   n_phi2 = phi2_counter$n)
      })
      out <- do.call(rbind, out)
      g1 <- cell$delta1 / cell$sd1
      g2 <- cell$delta2 / cell$sd2
      kappa <- max(g1, g2) / min(g1, g2)
      ck <- newton_ck(kappa = kappa, rho = cell$rho, alpha = alpha,
                      beta = 1 - cell$power)
      c_scale <- (1 + 1 / cell$r) / min(g1, g2) ^ 2
      n2_exact <- c_scale *
        lambda_exact_any(kappa, cell$rho, 1 - cell$power, alpha) ^ 2
      cbind(cell_id = i, cell[rep(1, nrow(out)), ], kappa = kappa, out,
            newton_iterations = ck$iterations, newton_n_eval = ck$n_eval,
            newton_converged = ck$converged, n2_exact = n2_exact,
            row.names = NULL)
    }))
    suppressMessages(untrace("pbivnorm", where = asNamespace("pbivnorm")))
    per_cell$rel_n2 <- per_cell$n2_cont / per_cell$n2_exact - 1

    # The counts are meaningful only if trace() caught the calls made inside
    # the package: the closed form makes none, the other methods make some,
    # and the count for the Newton-Raphson algorithm equals its own count
    is_m <- function(m) per_cell$method == m
    if (any(per_cell$n_phi2[is_m("univariate")] != 0L) ||
        any(per_cell$n_phi2[!is_m("univariate")] == 0L) ||
        any(per_cell$n_phi2[is_m("newton")] !=
              per_cell$newton_n_eval[is_m("newton")])) {
      stop("The evaluations of pbivnorm::pbivnorm() were not counted ",
           "correctly by trace().", call. = FALSE)
    }

    # Computing time per design: each method is run over the 96
    # combinations 'reps' times, and the elapsed time is divided by the
    # number of runs. The rows of 'cells' are extracted beforehand, so that
    # the time of the indexing is not counted.
    cell_list <- lapply(seq_len(nrow(cells)), function(i) as.list(cells[i, ]))
    reps <- c(univariate = 200L, bivariate = 5L, newton = 50L,
              uniroot = 50L)
    timing <- do.call(rbind, lapply(METHODS, function(m) {
      t0 <- proc.time()[["elapsed"]]
      for (k in seq_len(reps[[m]])) {
        for (cell in cell_list) ss_method(cell, m)
      }
      el <- proc.time()[["elapsed"]] - t0
      data.frame(method = m, reps = reps[[m]], runs = reps[[m]] * nrow(cells),
                 elapsed_sec = el,
                 ms_per_design = 1000 * el / (reps[[m]] * nrow(cells)))
    }))

    # Convergence of the Newton-Raphson algorithm at extreme correlations
    ext <- expand.grid(kappa = c(1, 1.5, 2, 3),
                       rho = c(-0.99, -0.9, 0.95, 0.99, 0.999),
                       power = c(0.80, 0.90))
    ext_res <- t(mapply(function(k, rr, pw) {
      ck <- newton_ck(kappa = k, rho = rr, alpha = alpha, beta = 1 - pw)
      lam_ex <- lambda_exact_any(k, rr, 1 - pw, alpha)
      lam_cf <- lambda_package(k, rr, 1 - pw, alpha)
      c(newton_converged = ck$converged, newton_iterations = ck$iterations,
        newton_n_eval = ck$n_eval,
        newton_rel_n2 = (ck$lambda / lam_ex) ^ 2 - 1,
        cf_rel_n2 = (lam_cf / lam_ex) ^ 2 - 1)
    }, ext$kappa, ext$rho, ext$power))
    ext <- cbind(ext, ext_res)

    list(data = per_cell, timing = timing, extreme = ext,
         uniroot_tol = UNIROOT_TOL, alpha = alpha)
  })
  save_data(cost, "cost_comparison")
}

# =============================================================================
# Noninteger and unbalanced allocation ratios
# =============================================================================
if (RUN_SI[["alloc"]]) {
  alloc <- run_step("Allocation ratios", {
    patterns <- data.frame(delta1 = c(0.3, 0.3, 0.7, 0.7),
                           delta2 = c(0.3, 0.3, 0.5, 0.5),
                           sd1 = 1, sd2 = c(1, 1.5, 1, 1.5))
    al <- expand.grid(r = c(0.25, 0.5, 1.5, 2.5, 4), pattern = 1:4,
                      rho = c(0, 0.5, 0.8), power = c(0.80, 0.90))
    res <- do.call(rbind, lapply(seq_len(nrow(al)), function(i) {
      a <- al[i, ]
      p <- patterns[a$pattern, ]
      ss <- lapply(c("bivariate", "univariate", "newton"), function(m) {
        twocpcont_ss(delta1 = p$delta1, delta2 = p$delta2, sd1 = p$sd1,
                     sd2 = p$sd2, rho = a$rho, r = a$r, alpha = alpha,
                     beta = 1 - a$power, method = m)
      })
      pw_prop <- twocpcont_power(n1 = ss[[2]]$n1, n2 = ss[[2]]$n2,
                                 delta1 = p$delta1, delta2 = p$delta2,
                                 sd1 = p$sd1, sd2 = p$sd2, rho = a$rho,
                                 alpha = alpha)$powerCoprimary
      # The exact continuous n2, c lambda_w^2, rounded as in (eq:ss_final):
      # it separates the effect of the ceiling n1 = ceiling(r n2) from that
      # of the approximation of lambda_w
      g1 <- p$delta1 / p$sd1
      g2 <- p$delta2 / p$sd2
      n2_cont_exact <- (1 + 1 / a$r) / min(g1, g2) ^ 2 *
        lambda_exact_any(max(g1, g2) / min(g1, g2), a$rho, 1 - a$power,
                         alpha) ^ 2
      n2_round_exact <- ceiling(n2_cont_exact)
      data.frame(a, p, n1_conv = ss[[1]]$n1, n2_conv = ss[[1]]$n2,
                 N_conv = ss[[1]]$N, n1_prop = ss[[2]]$n1,
                 n2_prop = ss[[2]]$n2, N_prop = ss[[2]]$N,
                 N_newton = ss[[3]]$N, power_prop = pw_prop,
                 n2_cont_prop = ss[[2]]$n2_cont,
                 n2_cont_exact = n2_cont_exact,
                 N_round_exact = ceiling(a$r * n2_round_exact) +
                   n2_round_exact,
                 row.names = NULL)
    }))
    list(data = res, alpha = alpha)
  })
  save_data(alloc, "allocation")
}

# =============================================================================
# EXPEDITION 1: sensitivity to the effect size of ADCS-ADL and to rho, and
# power when the correlation assumed at the planning stage is not the true
# correlation
# =============================================================================
if (RUN_SI[["sens"]]) {
  sens <- run_step("EXPEDITION 1 sensitivity", {
    beta_t <- 0.13
    es1    <- 1.8 / 9
    # Effect size of ADCS-ADL stated in the protocol and that of Table 4
    es2_protocol <- 0.25
    es2_design   <- 3.1 / 12
    sg <- expand.grid(es2 = round(seq(0.20, 0.30, by = 0.005), 3),
                      rho = c(0, 0.3, 0.5, 0.8))
    sg$kappa <- sg$es2 / es1
    sg$N <- mapply(function(e2, rr) {
      twocpcont_ss(delta1 = 1.8, delta2 = e2 * 12, sd1 = 9, sd2 = 12,
                   rho = rr, r = 1, alpha = alpha, beta = beta_t,
                   method = "univariate")$N
    }, sg$es2, sg$rho)
    sg$lambda <- mapply(function(k, rr) {
      lambda_exact_any(k, rr, beta_t, alpha)
    }, sg$kappa, sg$rho)
    N0   <- sg$N[sg$rho == 0]
    lam0 <- sg$lambda[sg$rho == 0]
    names(N0) <- names(lam0) <- sg$es2[sg$rho == 0]
    sg$reduction_int  <- 1 - sg$N / N0[as.character(sg$es2)]
    sg$reduction_cont <- 1 - sg$lambda ^ 2 / lam0[as.character(sg$es2)] ^ 2
    rmax <- data.frame(es2 = unique(sg$es2), kappa = unique(sg$kappa))
    rmax$R_max <- sapply(rmax$kappa, function(k) {
      r_max(kappa = k, alpha = alpha, beta = beta_t)
    })
    thresholds <- data.frame(
      epsilon = c(0.05, 0.01),
      kappa_star = sapply(c(0.05, 0.01), function(e) {
        kappa_star(epsilon = e, alpha = alpha, beta = beta_t)
      }))

    # Power at the true rho of the sample size planned at rho_plan
    designs <- data.frame(
      design = c("kappa1", "expedition1"),
      delta1 = c(0.3, 1.8), delta2 = c(0.3, 3.1), sd1 = c(1, 9),
      sd2 = c(1, 12), power = c(0.80, 0.87))
    pl <- expand.grid(design = designs$design, rho_plan = c(0, 0.3, 0.5, 0.8),
                      rho_true = c(-0.3, 0, 0.3, 0.5, 0.8),
                      stringsAsFactors = FALSE)
    pl <- merge(pl, designs, by = "design", sort = FALSE)
    plan <- t(mapply(function(d1, d2, s1, s2, pw, rp, rt) {
      ss <- twocpcont_ss(delta1 = d1, delta2 = d2, sd1 = s1, sd2 = s2,
                         rho = rp, r = 1, alpha = alpha, beta = 1 - pw,
                         method = "univariate")
      c(N_plan = ss$N,
        power_true = twocpcont_power(n1 = ss$n1, n2 = ss$n2, delta1 = d1,
                                     delta2 = d2, sd1 = s1, sd2 = s2,
                                     rho = rt, alpha = alpha)$powerCoprimary)
    }, pl$delta1, pl$delta2, pl$sd1, pl$sd2, pl$power, pl$rho_plan,
    pl$rho_true))
    pl <- cbind(pl, plan)
    pl <- pl[order(pl$design, pl$rho_plan, pl$rho_true), ]
    rownames(pl) <- NULL

    list(data = sg, r_max = rmax, thresholds = thresholds, planning = pl,
         es1 = es1, es2_protocol = es2_protocol, es2_design = es2_design,
         beta = beta_t, alpha = alpha)
  })
  save_data(sens, "expedition_sensitivity")
}

# =============================================================================
# Record the execution environment for reproducibility
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_supplement.txt")
)
log_msg("All requested data written to '%s/'.", data_dir)

# A RUN_SI assigned for this run must not carry over to the next one
rm(RUN_SI)
