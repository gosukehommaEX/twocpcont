# =============================================================================
# validate_claims_numerical.R
#
# Checks of the statements of Section 3 (numerical investigation) of the
# twocpcont article that are not single numbers: which design combinations
# meet the target power, monotonicity in the correlation, the direction of
# the conventional search, the sign of the discriminant of (eq:Delta_star),
# bounds on the deviations from the published tables, and the structure of
# those tables. The numbers cited in the text are checked against the ledger
# by inst/scripts/number_check.R.
#
# Claims about the results of Section 3 are checked on the data written by
# inst/scripts/data_generate_main.R; the starting value of the search and the
# discriminant are computed here.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file after data_generate_main.R:
#     devtools::load_all()
#     source("inst/validation/validate_claims_numerical.R")
#   For a test of the script itself, SCRIPTS_DATA_DIR (a copy of
#   inst/scripts/data) and VALIDATION_OUT_DIR can be assigned before
#   sourcing; both are removed when the script starts.
#
# Output  inst/validation/results/validate_claims_numerical.csv
#         (item, claim, verdict PASS, FAIL or INFO, detail; INFO when a
#         check cannot be made, e.g. without the manuscript folder)
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
# ok = NA gives the verdict INFO (the check could not be made)
check <- function(item, claim, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    item = item, claim = claim,
    verdict = if (identical(ok, NA)) "INFO" else if (isTRUE(ok)) "PASS" else
      "FAIL",
    detail = detail, stringsAsFactors = FALSE)
}

tab2 <- read_scripts_data("tab2_sample_size")
grid <- tab2$data
alpha     <- tab2$alpha
power_set <- tab2$power_set
z_alpha   <- qnorm(1 - alpha)

# =============================================================================
# Section 3.1: setting
# =============================================================================
# The R version stated in Section 3.1 of the manuscript ("R version x.y.z")
# must be the version with which every data file of inst/scripts/data was
# computed. Without the manuscript folder (e.g., in a clone of the
# repository) the verdict is INFO.
manuscript_main <- file.path("..", "..", "Manuscript", "Main_body", "main.tex")
r_data <- unique(vapply(
  list.files(scripts_data_dir, pattern = "\\.rds$", full.names = TRUE),
  function(f) readRDS(f)$meta$R_version, character(1L)))
r_data_num <- sub("^R version ([0-9.]+) .*$", "\\1", r_data)
if (file.exists(manuscript_main)) {
  ms <- readLines(manuscript_main, encoding = "UTF-8", warn = FALSE)
  ms <- gsub("(?<!\\\\)%.*$", "", ms, perl = TRUE)
  r_text <- unique(unlist(regmatches(
    ms, gregexpr("(?<=R\\} version )[0-9]+\\.[0-9]+\\.[0-9]+", ms,
                 perl = TRUE))))
  check("r_version",
        paste0("The R version stated in the manuscript is the version with ",
               "which all data of inst/scripts/data were computed"),
        length(r_data) == 1L && length(r_text) == 1L &&
          identical(r_text, r_data_num),
        sprintf("manuscript: %s; data: %s",
                if (length(r_text) == 0L) "(not found)" else
                  paste(r_text, collapse = ", "),
                paste(r_data, collapse = "; ")))
} else {
  check("r_version",
        paste0("The R version stated in the manuscript is the version with ",
               "which all data of inst/scripts/data were computed"),
        NA, sprintf("%s not found; data: %s", manuscript_main,
                    paste(r_data, collapse = "; ")))
}

N_fix  <- 600
n1_seq <- 1:(N_fix - 1)
check("equal_allocation_minimizes_N",
      "1/n1 + 1/n2 at fixed N is minimized at n1 = n2",
      n1_seq[which.min(1 / n1_seq + 1 / (N_fix - n1_seq))] == N_fix / 2,
      sprintf("N = %d", N_fix))

# =============================================================================
# Section 3.2: sample size comparison and achieved power
# =============================================================================
same <- grid$N_prop == grid$N_conv
mis  <- grid[!same, , drop = FALSE]
check("single_mismatch_cell",
      paste0("N_prop differs from N_conv only at (0.3, 0.3, 1, 1, r = 1, ",
             "rho = 0.5), 1 - beta = 0.90, where N_prop < N_conv"),
      nrow(mis) == 1L && mis$delta1 == 0.3 && mis$delta2 == 0.3 &&
        mis$sd1 == 1 && mis$sd2 == 1 && mis$r == 1 && mis$rho == 0.5 &&
        abs(mis$power - 0.90) < 1e-9 && mis$N_prop < mis$N_conv,
      sprintf("%d mismatch(es); N_prop = %s, N_conv = %s", nrow(mis),
              paste(mis$N_prop, collapse = ", "),
              paste(mis$N_conv, collapse = ", ")))

shortfall <- mis$power - mis$pow
check("mismatch_shortfall_about_1e-6",
      "The achieved power in that cell is short of the target by about 1e-6",
      length(shortfall) == 1L && shortfall > 0 &&
        round(log10(shortfall)) == -6,
      sprintf("shortfall = %s", paste(sprintf("%.2e", shortfall),
                                      collapse = ", ")))

check("other_cells_meet_target",
      "In the other combinations the achieved power is at least the target",
      all(grid$pow[same] >= grid$power[same]),
      sprintf("cells = %d, min(power - target) = %.2e", sum(same),
              min(grid$pow[same] - grid$power[same])))

series <- unique(grid[, c("delta1", "delta2", "sd1", "sd2", "r", "power")])
N_series <- lapply(seq_len(nrow(series)), function(i) {
  s <- series[i, ]
  d <- grid[grid$delta1 == s$delta1 & grid$delta2 == s$delta2 &
              grid$sd1 == s$sd1 & grid$sd2 == s$sd2 & grid$r == s$r &
              grid$power == s$power, ]
  d$N_prop[order(d$rho)]
})
check("N_prop_nonincreasing_in_rho",
      "N_prop is non-increasing in rho in every pattern and target power",
      length(N_series) == 24L && all(sapply(N_series, function(v) {
        all(diff(v) <= 0)
      })),
      sprintf("series = %d", length(N_series)))

k21 <- series$delta1 == 0.7 & series$delta2 == 0.5 & series$sd2 == 1.5 &
  series$r == 1
check("kappa21_constant_in_rho",
      "For (0.7, 0.5, 1, 1.5, r = 1) with kappa = 2.1, N does not change with rho",
      sum(k21) == 3L && all(sapply(N_series[k21], function(v) {
        length(unique(v)) == 1L
      })),
      paste(sapply(N_series[k21], paste, collapse = "/"), collapse = "; "))

red_frac <- function(v) (v[1] - v[length(v)]) / v[1]
mid <- (series$delta1 == 0.7 & series$delta2 == 0.5 & series$sd2 == 1) |
  (series$delta1 == 0.3 & series$delta2 == 0.3 & series$sd2 == 1.5)
red_mid <- sapply(N_series[mid], red_frac)
check("reduction_kappa14_15_at_most_4pct",
      "For kappa = 1.4 and 1.5 the reduction from rho = 0 to 0.8 is at most 4%",
      sum(mid) == 12L && max(red_mid) <= 0.04 + 1e-12,
      sprintf("max reduction = %.4f%%", 100 * max(red_mid)))

r2 <- grid$r == 2
check("r2_all_identical",
      "All combinations with r = 2 give identical sample sizes",
      sum(r2) == 48L && all(same[r2]),
      sprintf("cells = %d, identical = %d", sum(r2), sum(same[r2])))

# =============================================================================
# Section 3.3: published tables and computational cost
# =============================================================================
pub <- read_scripts_data("published_tables")$sozu
lo  <- pub$rho <= 0.8
check("exact_root_reproduces_all_published",
      "The exact root reproduces all tabulated values to three decimals",
      all(abs(pub$dev_exact) < 1e-9),
      sprintf("values = %d, reproduced = %d", nrow(pub),
              sum(abs(pub$dev_exact) < 1e-9)))

off_lo <- lo & abs(pub$dev_cf) > 1e-9
check("closed_form_other_values_differ_by_0.001",
      paste0("For rho <= 0.8, the closed form differs from the values it ",
             "does not reproduce by one unit in the third decimal"),
      sum(off_lo) > 0 && all(abs(abs(pub$dev_cf[off_lo]) - 0.001) < 1e-9),
      sprintf("values not reproduced = %d, deviations = %s", sum(off_lo),
              paste(sprintf("%+.3f", pub$dev_cf[off_lo]), collapse = ", ")))

max_err <- max(abs(pub$err_cf[lo]))
check("closed_form_lambda_error_below_4e-4",
      "For rho <= 0.8, |lambda_w* - exact root| < 4e-4",
      max_err < 4e-4, sprintf("max error = %.3e", max_err))

rel_n2 <- max(abs(((pub$C2_cf[lo] + z_alpha) /
                     (pub$C2_exact[lo] + z_alpha)) ^ 2 - 1))
check("closed_form_n2_change_below_0.03pct",
      "For rho <= 0.8, c (lambda_w*)^2 changes by less than 0.03%",
      rel_n2 < 3e-4, sprintf("max relative change = %.4f%%", 100 * rel_n2))

sozu_tab <- utils::read.csv(file.path("inst", "extdata",
                                      "sozu2015_C2_table.csv"))
check("published_table_structure",
      paste0("The tables cover 1 - beta in {0.80, 0.90}, 15 effect size ",
             "ratios between 1.00 and 2.00 and 7 correlations"),
      nrow(sozu_tab) == 210L &&
        setequal(unique(sozu_tab$power), c(0.80, 0.90)) &&
        length(unique(sozu_tab$gamma1)) == 15L &&
        abs(min(sozu_tab$gamma1) - 1) < 1e-9 &&
        abs(max(sozu_tab$gamma1) - 2) < 1e-9 &&
        setequal(unique(sozu_tab$rho), c(0, 0.2, 0.3, 0.5, 0.7, 0.8, 0.95)),
      sprintf("powers = %s; ratios = %d; correlations = %s",
              paste(unique(sozu_tab$power), collapse = ", "),
              length(unique(sozu_tab$gamma1)),
              paste(unique(sozu_tab$rho), collapse = ", ")))

tab4 <- read_scripts_data("tab4_real_example")
check("expedition1_not_in_tables",
      paste0("The design of Section 4 (1 - beta = 0.87, kappa = 1.292) is ",
             "not covered by the tables"),
      all(abs(unique(sozu_tab$power) - (1 - tab4$beta)) > 1e-9) &&
        all(abs(unique(sozu_tab$gamma1) - tab4$kappa) > 1e-3),
      sprintf("kappa = %.4f between tabulated %.2f and %.2f", tab4$kappa,
              max(sozu_tab$gamma1[sozu_tab$gamma1 < tab4$kappa]),
              min(sozu_tab$gamma1[sozu_tab$gamma1 > tab4$kappa])))

ev <- read_scripts_data("phi2_evaluations")$data
top5 <- order(-ev$n_phi2_conv)[1:5]
check("phi2_largest_counts_at_kappa_1.5",
      "The largest numbers of Phi2 evaluations occur at kappa = 1.5",
      all(abs(ev$kappa[top5] - 1.5) < 1e-9),
      sprintf("kappa of the five largest counts: %s",
              paste(sprintf("%.2f", ev$kappa[top5]), collapse = ", ")))
check("closed_form_no_phi2_evaluation",
      "The closed form evaluates Phi2 at no point",
      nrow(ev) == 96L && all(ev$n_phi2_prop == 0L),
      sprintf("cells = %d, max count = %d", nrow(ev), max(ev$n_phi2_prop)))

# For rho >= 0 the starting value (eq:n2_init) attains the target power, so
# the search of Algorithm 1 proceeds downward
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
check("search_proceeds_downward",
      "The starting value attains the target power in all combinations",
      length(start_ok) == 96L && all(start_ok),
      sprintf("cells = %d, downward = %d", length(start_ok), sum(start_ok)))

# The discriminant of (eq:Delta_star) is nonnegative in every design of
# Sections 3 and 4: the 96 cells of Table 2, the 210 tabulated designs of
# Sozu et al. (2015), the four EXPEDITION 1 designs, and the rho = 0 designs
# behind the closed-form R_max (kappa in [1, 3] by 0.001, 1 - beta in
# {0.70, 0.80, 0.85, 0.87, 0.90})
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
disc_grid  <- mapply(disc_fun, grid$kappa, grid$rho, grid$power)
disc_pub   <- mapply(disc_fun, pub$gamma1, pub$rho, pub$power)
disc_trial <- sapply(tab4$rho_set, function(rr) {
  disc_fun(tab4$kappa, rr, 1 - tab4$beta)
})
disc_rho0 <- unlist(lapply(c(0.70, 0.80, 0.85, 0.87, 0.90), function(tp) {
  sapply(seq(1, 3, by = 0.001), disc_fun, rr = 0, tp = tp)
}))
disc_all <- c(disc_grid, disc_pub, disc_trial, disc_rho0)
check("discriminant_nonnegative",
      paste0("The discriminant of (eq:Delta_star) is nonnegative in every ",
             "design of Sections 3 and 4"),
      length(disc_grid) == 96L && length(disc_pub) == 210L &&
        all(disc_all >= 0),
      sprintf("designs = %d, min discriminant = %.4f", length(disc_all),
              min(disc_all)))

# =============================================================================
# Section 3.4: threshold effect size ratio
# =============================================================================
kst <- read_scripts_data("tab3_kappa_star")$data
ks05 <- kst$kappa_star[abs(kst$epsilon - 0.05) < 1e-9]
check("kstar_05_at_most_1.4",
      paste0("kappa_star(0.05, beta) is at most 1.4 at the three target ",
             "powers"),
      length(ks05) == 3L && max(ks05) <= 1.4,
      sprintf("kappa_star(0.05) = %s",
              paste(sprintf("%.4f", ks05), collapse = ", ")))

# =============================================================================
# Results
# =============================================================================
res <- do.call(rbind, checks)
out_file <- file.path(out_dir, "validate_claims_numerical.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf(
  "validate_claims_numerical.R: PASS %d, FAIL %d, INFO %d; written to %s",
  sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"),
  sum(res$verdict == "INFO"), out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("validate_claims_numerical: FAIL in %s.",
               paste(res$item[res$verdict == "FAIL"], collapse = ", ")),
       call. = FALSE)
}
