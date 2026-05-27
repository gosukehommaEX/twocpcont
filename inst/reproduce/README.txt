================================================================
README for the code supplement of the manuscript

  "A Noniterative Closed-Form Sample Size Formula for Clinical
   Trials with Two Co-Primary Continuous Endpoints"

  Author:  Gosuke Homma
  Contact: my.name.is.gosuke@gmail.com
================================================================


1. OVERVIEW
----------------------------------------------------------------

This zip file contains all R code required to reproduce every
numerical result, table and figure in the manuscript. The same
code is also available as an R package "twocpcont" on GitHub
(https://github.com/gosukehommaEX/twocpcont), but the present zip
file is fully self-contained and does not require any
installation from GitHub.


2. CONTENTS
----------------------------------------------------------------

README.txt   This file.

Rcode/       Folder containing eleven R source files.

The eight files below define the functions used in the analysis.
Each file contains a single function, documented in roxygen2
format.

  GL_nodes_and_weights.R
      Gauss-Legendre nodes and weights on [0, rho].
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

The three files below are the reproduction scripts.

  run_table_and_figure_manuscript.R
      Master script. Performs all numerical computations and
      saves intermediate results (csv, rds) inside Rcode/.
  create_table_and_figure_manuscript.R
      Reads the intermediate results and renders the tables
      (tex) and figures (pdf, eps) inside Rcode/.
  verify_intext_numbers.R
      Verifies the in-text numerical claims of Sections 3 and 4
      of the manuscript. Each check prints a PASS / FAIL line.


3. HOW TO REPRODUCE
----------------------------------------------------------------

(1) Unzip and copy the Rcode/ folder to a local directory of
    your choice.

(2) In R or RStudio, set the working directory to that Rcode/
    folder:

      setwd("path/to/Rcode")

(3) Install the two required CRAN packages if not yet installed:

      install.packages("pbivnorm")
      install.packages("ggplot2")

(4) Run the three scripts in the following order:

      source("run_table_and_figure_manuscript.R")
      source("create_table_and_figure_manuscript.R")
      source("verify_intext_numbers.R")

All outputs are written into two subfolders that are created
automatically inside Rcode/:

  Rcode/table_and_figure_manuscript/
      Tables (csv, tex), intermediate results (rds) and figures
      (pdf, eps) produced by the first two scripts.
  Rcode/verify_intext_numbers_manuscript/
      The PASS / FAIL log (txt) produced by the third script.

The verification script also prints PASS / FAIL lines to the
console; all should read PASS.

The proposed method is deterministic, no random seed is needed,
and no manual edits to the scripts are required. All paths used
in the scripts are relative.


4. SESSION INFORMATION
----------------------------------------------------------------

Output of sessionInfo() after loading the two required packages:

  > library(pbivnorm)
  > library(ggplot2)
  > sessionInfo()
  R version 4.6.0 (2026-04-24 ucrt)
  Platform: x86_64-w64-mingw32/x64
  Running under: Windows 11 x64 (build 26200)

  Matrix products: default
    LAPACK version 3.12.1

  locale:
  [1] LC_COLLATE=Japanese_Japan.utf8
  [2] LC_CTYPE=Japanese_Japan.utf8
  [3] LC_MONETARY=Japanese_Japan.utf8
  [4] LC_NUMERIC=C
  [5] LC_TIME=Japanese_Japan.utf8

  time zone: Asia/Tokyo
  tzcode source: internal

  attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base

  other attached packages:
  [1] ggplot2_4.0.3  pbivnorm_0.6.0

================================================================
End of README
================================================================
