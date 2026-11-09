# =============================================================================
# copy_results.R
#
# Copies the figures (800 dpi EPS) and LaTeX tables of the manuscript and of
# the Supporting Information to the manuscript folders and checks that each
# copy is identical to its source (MD5 checksum). Use this script instead of
# copying the files by hand.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root) and
#   source this file after table_and_figure_manuscript.R and
#   table_and_figure_supplement.R:
#     source("inst/scripts/copy_results.R")
#   To copy to another manuscript folder (e.g., a new submission), assign
#   DEST_DIR before sourcing this file. A folder of a submitted version should
#   not be overwritten; create a new folder instead.
#
# Inputs (relative to the package root)
#   inst/scripts/results/   Figure*.eps, table*.tex   -> DEST_DIR/Main_body/
#   inst/sup_info/results/  FigureS*.eps, tableS*.tex -> DEST_DIR/Supp_info/
# Only EPS figures are copied: the manuscript uses
# \includegraphics{Figure1.eps}, and a PDF of a figure must not be placed in
# the manuscript folders (tools/build.py of the manuscript checks this).
#
# Output
#   The files in DEST_DIR/Main_body and DEST_DIR/Supp_info and the log
#   DEST_DIR/tools/logs/copy_results.log with the MD5 checksums
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "twocpcont")) {
  stop("Run this script from the root of the twocpcont package ",
       "(open twocpcont.Rproj first).", call. = FALSE)
}

if (!exists("DEST_DIR")) {
  DEST_DIR <- file.path("..", "..", "Manuscript")
}
TARGETS <- list(
  list(src = file.path("inst", "scripts", "results"),
       dest = file.path(DEST_DIR, "Main_body")),
  list(src = file.path("inst", "sup_info", "results"),
       dest = file.path(DEST_DIR, "Supp_info"))
)

src <- character(0)
dest <- character(0)
for (tg in TARGETS) {
  files <- list.files(tg$src, pattern = "\\.(eps|tex)$", full.names = TRUE)
  if (length(files) == 0L) next
  if (!dir.exists(tg$dest)) {
    stop(sprintf("The destination folder '%s' does not exist.", tg$dest),
         call. = FALSE)
  }
  src <- c(src, files)
  dest <- c(dest, file.path(tg$dest, basename(files)))
}
if (length(src) == 0L) {
  stop("No figures or tables found; run the rendering scripts first.",
       call. = FALSE)
}
if (anyDuplicated(basename(src))) {
  stop("Files with the same name exist in the result folders: ",
       paste(unique(basename(src)[duplicated(basename(src))]),
             collapse = ", "), call. = FALSE)
}

ok <- file.copy(src, dest, overwrite = TRUE, copy.date = TRUE)
if (!all(ok)) {
  stop("Copy failed: ", paste(basename(src)[!ok], collapse = ", "),
       call. = FALSE)
}

md5_src  <- unname(tools::md5sum(src))
md5_dest <- unname(tools::md5sum(dest))
res <- data.frame(file   = basename(src),
                  source = dirname(src),
                  dest   = basename(dirname(dest)),
                  md5    = md5_src,
                  status = ifelse(md5_src == md5_dest, "OK", "MISMATCH"),
                  stringsAsFactors = FALSE)

log_dir <- file.path(DEST_DIR, "tools", "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_path <- file.path(log_dir, "copy_results.log")
out <- c(sprintf("copy_results.R -- %s", format(Sys.time())),
         sprintf("Destination: %s", normalizePath(DEST_DIR)),
         "",
         utils::capture.output(print(res, row.names = FALSE, right = FALSE,
                                     max = nrow(res) * ncol(res) + 1L)))
writeLines(out, log_path)
writeLines(out)

if (any(res$status != "OK")) {
  stop("MD5 mismatch for: ",
       paste(res$file[res$status != "OK"], collapse = ", "), call. = FALSE)
}
message(sprintf("%d files copied to %s and verified by MD5.", nrow(res),
                DEST_DIR))
