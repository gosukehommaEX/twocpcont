# Scripts for the figures and tables of the article

This folder and `inst/sup_info/` contain the code that produces every figure,
table and number of the article "A noniterative closed-form sample size
formula for clinical trials with two co-primary continuous endpoints".
Computation and rendering are separated, so that figures and tables can be
changed without repeating the computation.

All scripts are run from the package root (open `twocpcont.Rproj`).

## Folder structure

```
inst/scripts/
  data_generate_main.R           computation  -> data/*.rds, data/csv/*.csv,
                                                 data/data_generate_main.log
  table_and_figure_manuscript.R  rendering    -> results/figure*.eps,
                                                 results/table*.tex,
                                                 results/manuscript_numbers.csv
  number_check.R                 check of the numbers in the manuscript
                                              -> number_check_output.txt
  number_locations.csv           positions of the numbers in the manuscript
                                 (maintained by hand)
  copy_results.R                 copy of figures and tables to the manuscript
                                 folders with MD5 verification
inst/sup_info/
  data_generate_supplement.R     computation  -> data/*.rds, data/csv/*.csv,
                                                 data/data_generate_supplement.log
  table_and_figure_supplement.R  rendering    -> results/FigureS*.eps,
                                                 results/tableS*.tex,
                                                 results/appendix_software_example.tex,
                                                 results/supplement_numbers.csv
inst/reproduce/                  reproduction of published values
                                 (Sozu et al., 2015; Hung and Wang, 2009;
                                 Sugimoto et al., 2012)
inst/validation/                 checks of the statements of the article that
                                 are not single numbers
```

## Order of execution

```r
devtools::load_all()
source("inst/scripts/data_generate_main.R")              # about 1 minute
source("inst/scripts/table_and_figure_manuscript.R")     # a few seconds
source("inst/sup_info/data_generate_supplement.R")       # a few minutes
source("inst/sup_info/table_and_figure_supplement.R")    # a few seconds
source("inst/reproduce/reproduce_Sozu2015.R")            # a few seconds
source("inst/reproduce/reproduce_HungWang2009.R")        # a few seconds
source("inst/reproduce/reproduce_Sugimoto2012.R")        # a few seconds
source("inst/validation/validate_newton_mvtnorm.R")      # less than 1 minute
source("inst/validation/validate_claims_methods.R")      # about 1 minute
source("inst/validation/validate_claims_numerical.R")    # about 1 minute
source("inst/validation/validate_claims_example.R")      # a few seconds
source("inst/validation/validate_error_decomposition.R") # less than 1 minute
source("inst/scripts/copy_results.R")
source("inst/scripts/number_check.R")
```

The times are rough guides; the time of each step of
`data_generate_main.R` is written to `data/data_generate_main.log`, and that
of each part of `data_generate_supplement.R` to
`inst/sup_info/data/data_generate_supplement.log`. The computing times in
Table S1 depend on the computer.

## Random numbers

No random numbers are used: every result is a deterministic function of the
design parameters, and no seed is set.

## Results and their sources

Each table, figure and number of the article, the script that produces it and
the data file it is computed from.

| Result in the article | Produced by | From |
|---|---|---|
| Table 1 (Gauss-Legendre nodes and weights) | `table_and_figure_manuscript.R` | `data/tab1_gl_nodes.rds` |
| Table 2 (sample sizes of the two methods, 96 combinations) | `table_and_figure_manuscript.R` | `data/tab2_sample_size.rds` |
| Table 3 (threshold effect size ratio) | `table_and_figure_manuscript.R` | `data/tab3_kappa_star.rds` |
| Table 4 (EXPEDITION 1) | `table_and_figure_manuscript.R` | `data/tab4_real_example.rds` |
| Figure 1 (anchor weight) | `table_and_figure_manuscript.R` | `data/fig1_anchor_weight.rds` |
| Figure 2 (maximum reduction rate) | `table_and_figure_manuscript.R` | `data/fig2_reduction_curve.rds` |
| Figure 3 (reduction of the sample size over the correlation) | `table_and_figure_manuscript.R` | `data/fig3_reduction_rho.rds` |
| Comparison with the published tables (Section 3.3) | `table_and_figure_manuscript.R` (ledger) | `data/published_tables.rds` |
| Evaluations of the bivariate normal CDF (Section 3.3) | `table_and_figure_manuscript.R` (ledger) | `data/phi2_evaluations.rds` |
| Numbers in the text | `table_and_figure_manuscript.R` (`add_num()`, `add_diff()`) | `results/manuscript_numbers.csv` |
| Code and output of Section S5 | `table_and_figure_supplement.R` | `inst/sup_info/data/software_example.rds` |
| Table S1 (comparison with the Newton-Raphson algorithm and Brent's method) | `table_and_figure_supplement.R` | `inst/sup_info/data/cost_comparison.rds` |
| Figure S1 (error of the continuous n2 by approximation step) | `table_and_figure_supplement.R` | `inst/sup_info/data/error_decomposition.rds` |
| Table S2 and Figure S2 (accuracy over a wide grid) | `table_and_figure_supplement.R` | `inst/sup_info/data/accuracy_grid.rds` |
| Table S3 (allocation ratios) | `table_and_figure_supplement.R` | `inst/sup_info/data/allocation.rds` |
| Table S4 and Figure S3 (EXPEDITION 1: planning correlation and effect size of ADCS-ADL) | `table_and_figure_supplement.R` | `inst/sup_info/data/expedition_sensitivity.rds` |
| Numbers in the Supporting Information | `table_and_figure_supplement.R` (`add_num()`) | `inst/sup_info/results/supplement_numbers.csv` |
| Statements in the text that are not single numbers (bounds, monotonicity, which combinations meet the target power) | `inst/validation/validate_claims_*.R` | `inst/validation/results/*.csv` |
| Reproduction of the published values of Sozu et al. (2015), Hung and Wang (2009) and Sugimoto et al. (2012) | `inst/reproduce/reproduce_*.R` | `inst/reproduce/results/*.csv` |
