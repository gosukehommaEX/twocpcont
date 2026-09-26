# Build the flat code supplement for the CSSC submission from the
# package sources, so that the supplement never drifts from the
# package.
#
# Run from the package root (the folder containing DESCRIPTION):
#   source("tools/make_supplement.R")
#
# Output:
#   ../../Manuscript/CSSC/Rcode/twocpcont_supplement/
#     README.txt   generated from tools/supplement_README.txt, with the
#                  number of R files and the current sessionInfo()
#     Rcode/       copies of R/*.R, inst/reproduce/*.R, and
#                  inst/extdata/*.csv (published values read by the
#                  reproduction scripts)

if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[1, 1] != "twocpcont") {
  stop("Run this script from the root of the twocpcont package.")
}

supp_dir  <- file.path("..", "..", "Manuscript", "CSSC", "Rcode",
                       "twocpcont_supplement")
rcode_dir <- file.path(supp_dir, "Rcode")
dir.create(rcode_dir, recursive = TRUE, showWarnings = FALSE)

# Remove old R and csv files so that deleted files do not linger
old_files <- list.files(rcode_dir, pattern = "\\.(R|csv)$", full.names = TRUE)
if (length(old_files) > 0L) file.remove(old_files)

src_files <- c(list.files("R", pattern = "\\.R$", full.names = TRUE),
               list.files(file.path("inst", "reproduce"),
                          pattern = "\\.R$", full.names = TRUE))
data_files <- list.files(file.path("inst", "extdata"),
                         pattern = "\\.csv$", full.names = TRUE)
all_files <- c(src_files, data_files)
ok <- file.copy(all_files, rcode_dir, overwrite = TRUE)
if (!all(ok)) stop("Failed to copy: ", paste(all_files[!ok], collapse = ", "))

# Confirm that every copy is identical to its source
dst_files <- file.path(rcode_dir, basename(all_files))
same <- unname(tools::md5sum(all_files)) == unname(tools::md5sum(dst_files))
if (!all(same)) {
  stop("Copies differ from sources: ",
       paste(basename(all_files)[!same], collapse = ", "))
}

# README with the number of files and the session information.
# sessionInfo() is taken in a fresh R process, so that packages loaded
# in the current session (e.g., by devtools::load_all()) are not listed.
si_script <- tempfile(fileext = ".R")
writeLines(c("suppressPackageStartupMessages({",
             "  library(pbivnorm)",
             "  library(ggplot2)",
             "})",
             "print(sessionInfo())"), si_script)
rscript <- file.path(R.home("bin"), "Rscript")
si <- system2(rscript, shQuote(si_script), stdout = TRUE)
unlink(si_script)
if (!any(grepl("pbivnorm", si))) stop("Failed to obtain sessionInfo().")
si <- paste0("  ", si)
tmpl <- readLines(file.path("tools", "supplement_README.txt"),
                  encoding = "UTF-8")
tmpl <- sub("{{N_FILES}}", as.character(length(src_files)), tmpl,
            fixed = TRUE)
pos  <- which(tmpl == "{{SESSION_INFO}}")
if (length(pos) != 1L) stop("Placeholder {{SESSION_INFO}} not found.")
readme <- c(tmpl[seq_len(pos - 1L)], si, tmpl[-seq_len(pos)])
writeLines(readme, file.path(supp_dir, "README.txt"), useBytes = TRUE)

message(sprintf(paste0("Supplement written to %s (%d R files and %d csv ",
                       "files, all identical)."),
                normalizePath(supp_dir), length(src_files),
                length(data_files)))
