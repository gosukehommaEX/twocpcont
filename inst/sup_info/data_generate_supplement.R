# =============================================================================
# data_generate_supplement.R
#
# Computation stage for the Supporting Information of the twocpcont article.
# It writes the results to inst/sup_info/data/, which
# table_and_figure_supplement.R reads. The structure and the helpers are
# those of inst/scripts/data_generate_main.R.
#
# Section S5 (software implementation) shows R code together with its
# output. The code is evaluated here, one console input at a time, and the
# console transcript is saved, so that the code and output printed in the
# Supporting Information always come from the current implementation and
# are never copied by hand.
#
# Usage
#   Open twocpcont.Rproj, load the package, and source this file:
#     devtools::load_all()
#     source("inst/sup_info/data_generate_supplement.R")
#
# Outputs (relative to the package root)
#   inst/sup_info/data/software_example.rds   Section S5 (code and output)
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

data_dir <- file.path("inst", "sup_info", "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# =============================================================================
# Section S5: software example for the EXPEDITION 1 design with rho = 0.5
# =============================================================================
# Code shown in Section S5, one element per console input. The first
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

writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_supplement.txt")
)
