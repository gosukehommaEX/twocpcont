================================================================
README for the code supplement of the manuscript

  "A Noniterative Closed-Form Sample Size Formula for Clinical
   Trials with Two Co-Primary Continuous Endpoints"

  Author:  Gosuke Homma
  Contact: my.name.is.gosuke@gmail.com
================================================================


1. OVERVIEW
----------------------------------------------------------------

This supplement contains all R code required to reproduce every
numerical result, table and figure in the manuscript. The same
code is available as the R package "twocpcont" on GitHub
(https://github.com/gosukehommaEX/twocpcont), but the supplement
is self-contained and does not require installing the package.


2. CONTENTS
----------------------------------------------------------------

README.txt   This file.

Rcode/       Folder containing {{N_FILES}} R source files and one
             csv file.

The files below define the functions used in the analysis. Each
file contains a single function, documented in roxygen2 format.

  GL_nodes_and_weights.R
      Gauss-Legendre nodes and weights on [-1, 1].
  plackett_gl_full.R
      Plackett correlation integral and its first two
      derivatives, evaluated by Gauss-Legendre quadrature.
  twocpcont_power.R
      Marginal and co-primary power for given sample sizes.
  twocpcont_ss.R
      Required sample size by the proposed closed-form formula
      or the conventional sequential search.
  r_max.R
      Maximum reduction rate Rmax(kappa, beta).
  kappa_star.R
      Threshold effect-size ratio kappa*(epsilon, beta).
  print_twocpcont_power.R
      print method for class "twocpcont_power".
  print_twocpcont_ss.R
      print method for class "twocpcont_ss".

The four files below are the reproduction scripts.

  run_table_and_figure_manuscript.R
      Performs all numerical computations and saves the results
      (rds) into the subfolder table_and_figure_manuscript/.
  create_table_and_figure_manuscript.R
      Reads the rds files and writes the tables (tex, csv) and
      figures (pdf, eps) into the same subfolder. It performs no
      numerical computation.
  verify_intext_numbers.R
      Checks the numbers cited in the text of the manuscript.
      Each check prints a PASS / FAIL line, and the log is written
      to the subfolder verify_intext_numbers_manuscript/.
  appendix_software_example.R
      Runs the software example of the appendix and writes the
      code and its output (appendix_software_example.tex) into
      table_and_figure_manuscript/.

The csv file below is read by run_table_and_figure_manuscript.R.

  sozu2015_C2_table.csv
      Constant C_2 tabulated in Tables 4.3 and 4.4 of Sozu et al.
      (2015), used to compare the exact and closed-form values
      with the published tables.


3. HOW TO REPRODUCE
----------------------------------------------------------------

(1) Copy the Rcode/ folder to a writable local directory.

(2) In R or RStudio, set the working directory to that Rcode/
    folder:

      setwd("path/to/Rcode")

(3) Install the two required CRAN packages if not yet installed:

      install.packages("pbivnorm")
      install.packages("ggplot2")

(4) Run the four scripts in the following order:

      source("run_table_and_figure_manuscript.R")
      source("create_table_and_figure_manuscript.R")
      source("verify_intext_numbers.R")
      source("appendix_software_example.R")

The scripts source the function files in Rcode/. All PASS / FAIL
lines printed by verify_intext_numbers.R should read PASS. The
computations are deterministic, so no random seed is needed, and
no manual edits to the scripts are required. All paths used in
the scripts are relative.


4. SESSION INFORMATION
----------------------------------------------------------------

Output of sessionInfo() after loading the two required packages:

{{SESSION_INFO}}

================================================================
End of README
================================================================
