================================================================
Reproduction scripts for the manuscript

  "A Noniterative Closed-Form Sample Size Formula for Clinical
   Trials with Two Co-Primary Continuous Endpoints"

  Author:  Gosuke Homma
  Contact: my.name.is.gosuke@gmail.com
================================================================


1. CONTENTS
----------------------------------------------------------------

  run_table_and_figure_manuscript.R
      Performs all numerical computations for the tables and
      figures and saves the results (rds) into the subfolder
      table_and_figure_manuscript/.
  create_table_and_figure_manuscript.R
      Reads the rds files and writes the tables (tex, csv) and
      figures (pdf, eps) into the same subfolder.  It performs
      no numerical computation.
  verify_intext_numbers.R
      Checks the numbers cited in the text of the manuscript.
      Each check prints a PASS / FAIL line, and the log is
      written to the subfolder verify_intext_numbers_manuscript/.
  appendix_software_example.R
      Runs the software example of the appendix and writes the
      code and its output (appendix_software_example.tex) into
      table_and_figure_manuscript/.


2. HOW TO REPRODUCE
----------------------------------------------------------------

(1) Install twocpcont and the packages used by the scripts:

      install.packages(c("pbivnorm", "ggplot2", "remotes"))
      remotes::install_github("gosukehommaEX/twocpcont")

(2) Copy the four scripts to a writable folder, because the
    scripts write their outputs into subfolders of the working
    directory.  The scripts are located at

      system.file("reproduce", package = "twocpcont")

(3) Set the working directory to that folder and run

      source("run_table_and_figure_manuscript.R")
      source("create_table_and_figure_manuscript.R")
      source("verify_intext_numbers.R")
      source("appendix_software_example.R")

All PASS / FAIL lines printed by verify_intext_numbers.R should read
PASS.  The computations are deterministic, so no random seed is
needed.

If the R source files of twocpcont (GL_nodes_and_weights.R,
plackett_gl_full.R, twocpcont_power.R, twocpcont_ss.R,
print_twocpcont_ss.R, print_twocpcont_power.R, r_max.R and
kappa_star.R) are placed in the same folder as the scripts, as in
the code supplement of the manuscript, the scripts source them
instead of loading the installed package.  In the same way,
run_table_and_figure_manuscript.R reads sozu2015_C2_table.csv (the
published tables of Sozu et al., 2015) from the working directory
when present and from inst/extdata of the installed package
otherwise.

================================================================
End of README
================================================================
