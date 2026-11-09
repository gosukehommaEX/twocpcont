# =============================================================================
# validate_newton_mvtnorm.R
#
# Validation of newton_ck() against the R function CKsolution of Sugimoto,
# Sozu and Hamasaki (2012), Pharmaceutical Statistics 11(2): 118-128,
# Supporting Information, Appendix A.2.1, transcribed below as printed
# (only the indentation is changed). CKsolution computes the bivariate
# normal probabilities with mvtnorm::pmvnorm and the initial value with
# mvtnorm::qmvnorm, whereas newton_ck() uses pbivnorm::pbivnorm and
# uniroot(); the two are independent routes for the same algorithm.
#
# Designs: the 210 designs of the tables of Sozu et al. (2015, Tables 4.3
# and 4.4; inst/extdata/sozu2015_C2_table.csv) and the 48 distinct
# combinations of (kappa, rho, 1 - beta) of Table 2 of the article.
#
# Verdicts (one row per design)
#   PASS  newton_ck() met its stopping rule, its power is at least the
#         target, and its constant differs from that of CKsolution by
#         less than 1e-4 (each stops within 1e-5 of the target power, i.e.
#         within about 4e-5 of the root on the scale of C)
#   FAIL  otherwise
# The exact root (uniroot() on pbivnorm::pbivnorm with tolerance 1e-13) is
# reported for information.
#
# Usage
#   Requires the mvtnorm package (install.packages("mvtnorm")).
#   Open twocpcont.Rproj (the working directory must be the package root),
#   load the package, and source this file:
#     devtools::load_all()
#     source("inst/validation/validate_newton_mvtnorm.R")
#   For a test of the script itself, VALIDATION_OUT_DIR can be assigned
#   before sourcing; it is removed when the script starts.
#
# Output  inst/validation/results/validate_newton_mvtnorm.csv
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
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("The mvtnorm package is required: install.packages(\"mvtnorm\").",
       call. = FALSE)
}

out_dir <- file.path("inst", "validation", "results")
if (exists("VALIDATION_OUT_DIR", envir = globalenv())) {
  out_dir <- get("VALIDATION_OUT_DIR", envir = globalenv())
  rm("VALIDATION_OUT_DIR", envir = globalenv())
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# CKsolution: Sugimoto et al. (2012), Supporting Information, Appendix A.2.1
# (transcribed as printed; R version 2.11.0 at the time of publication)
# -----------------------------------------------------------------------------
CKsolution <- function(alpha,power,rho,gamma,K){
  if(det(rho) <= 0){print("no positive definite");break}
  library(mvtnorm)
  z_a = qnorm(1-alpha); ndel = 0.001; rn = round(runif(1)*1000);
  # # Initial estimation of CK
  CK = qmvnorm(power, corr=rho, tail="lower.tail")$quantile
  # # Begin: Newton-Raphson algorithm to find CK
  for(j in 1:1000){
    set.seed(rn)
    C1k = CK*gamma + z_a*(gamma - 1)
    pow1 = pmvnorm(lower=rep(-Inf,K), upper=c(C1k,CK), corr=rho)[1]
    G = power - pow1
    if(abs(G) < 0.00001 & G <= 0){break}
    F1k = rep(0,K-1)
    for(l in 1:(K-1)){
      vndel = rep(0,K-1); vndel[l] = ndel
      F1k[l] = pmvnorm(lower=c(rep(-Inf,l-1),C1k[l],rep(-Inf,K-l)),
                       upper=c(C1k+vndel,CK), corr=rho)[1]/ndel }
    FK = pmvnorm(lower=c(rep(-Inf,K-1),CK),upper=c(C1k,CK+ndel),corr=rho)[1]/ndel
    dG = -t(F1k)%*%gamma - FK
    CK = CK - G/dG }
  # # End: Newton-Raphson algorithm to find CK
  return(c(CK))}

# -----------------------------------------------------------------------------
# Designs
# -----------------------------------------------------------------------------
alpha   <- 0.025
z_alpha <- qnorm(1 - alpha)
sozu_tab <- utils::read.csv(file.path("inst", "extdata",
                                      "sozu2015_C2_table.csv"))
designs <- rbind(
  data.frame(source = sprintf("Sozu et al. (2015), Table %s", sozu_tab$table),
             kappa = sozu_tab$gamma1, rho = sozu_tab$rho,
             power = sozu_tab$power, stringsAsFactors = FALSE),
  data.frame(source = "Table 2 of the article",
             expand.grid(kappa = c(1, 0.7 / 0.5, 0.3 / 0.2, 0.7 / (0.5 / 1.5)),
                         rho = c(0, 0.3, 0.5, 0.8),
                         power = c(0.80, 0.85, 0.90)),
             stringsAsFactors = FALSE)
)

# -----------------------------------------------------------------------------
# Comparison
# -----------------------------------------------------------------------------
res <- do.call(rbind, lapply(seq_len(nrow(designs)), function(i) {
  d <- designs[i, ]
  rho_mat <- matrix(c(1, d$rho, d$rho, 1), ncol = 2)
  C_orig <- CKsolution(alpha, d$power, rho_mat, d$kappa, 2)
  ck <- newton_ck(kappa = d$kappa, rho = d$rho, alpha = alpha,
                  beta = 1 - d$power)
  C_exact <- uniroot(function(cc) {
    pbivnorm::pbivnorm(d$kappa * cc + z_alpha * (d$kappa - 1), cc,
                       rho = d$rho) - d$power
  }, lower = -2, upper = 5, tol = 1e-13)$root
  pw_pkg <- pbivnorm::pbivnorm(d$kappa * ck$C + z_alpha * (d$kappa - 1),
                               ck$C, rho = d$rho)
  ok <- ck$converged && pw_pkg >= d$power && abs(ck$C - C_orig) < 1e-4
  data.frame(item = sprintf("%s: kappa = %.4f, rho = %.2f, 1 - beta = %.2f",
                            d$source, d$kappa, d$rho, d$power),
             C_CKsolution = C_orig, C_newton_ck = ck$C, C_exact = C_exact,
             difference = ck$C - C_orig, iterations = ck$iterations,
             verdict = if (ok) "PASS" else "FAIL",
             stringsAsFactors = FALSE)
}))

out_file <- file.path(out_dir, "validate_newton_mvtnorm.csv")
utils::write.csv(res, out_file, row.names = FALSE)
message(sprintf(paste0(
  "validate_newton_mvtnorm.R: PASS %d, FAIL %d; max |newton_ck - ",
  "CKsolution| = %.2e, max (newton_ck - exact) = %.2e; written to %s"),
  sum(res$verdict == "PASS"), sum(res$verdict == "FAIL"),
  max(abs(res$difference)), max(res$C_newton_ck - res$C_exact), out_file))
if (any(res$verdict == "FAIL")) {
  stop(sprintf("validate_newton_mvtnorm: %d FAIL verdict(s); see %s.",
               sum(res$verdict == "FAIL"), out_file), call. = FALSE)
}
