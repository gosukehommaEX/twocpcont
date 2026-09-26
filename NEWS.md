# twocpcont (development version)

## Changes

- `r_max()` gains an argument `exact`. With `exact = TRUE` (the new
  default), the noncentrality parameter at rho = 0 is the exact root of
  the factorized power equation, found by `uniroot()`, so that
  `kappa_star()` is the exact inverse of `r_max()`. With
  `exact = FALSE`, the closed-form approximation of version 0.1.0 is
  used. The two agree to within 5e-5 on the reduction-rate scale and
  coincide at kappa = 1.
- The documentation of `r_max()` now derives the monotonicity of
  R_max in kappa, which makes R_max(1, beta) an upper bound over all
  kappa >= 1.
- `twocpcont_ss()` returns a new column `n2_cont`, the closed-form
  value of n2 before the ceiling operation, c (lambda_w*)^2, when
  `method = "univariate"` (`NA` when `method = "bivariate"`). The
  existing columns and the computed sample sizes are unchanged.
- The print methods of `twocpcont_ss()` and `twocpcont_power()` label
  `delta1` and `delta2` as "Mean difference" instead of "Effect size",
  because "effect size" denotes the standardized effect size
  delta / sd elsewhere in the package and the manuscript.
- The published tables of the constant C_2 in Sozu et al. (2015,
  Tables 4.3 and 4.4) are included as
  `inst/extdata/sozu2015_C2_table.csv`, and
  `run_table_and_figure_manuscript.R` compares them with the exact and
  closed-form values.

## Bug fixes

- The unit test comparing `twocpcont_power()` across the two methods
  compared `NULL` with `NULL`, because `$power` partially matches three
  columns. It now compares the `powerCoprimary` column.

## Tests

- The single test file `test-twocpcont.R` was replaced by one test file
  per function (8 files, 73 tests). New tests cover the analytic
  derivatives and the boundary limits of `plackett_gl_full()`, the
  exactness of the Gauss-Legendre rule, agreement of the two methods
  on the 32-cell grid of the numerical study at 80%, 85% and 90%
  power, invariance of the sample size to the order and the scale of
  the endpoints, the monotonicity of `r_max()` and `kappa_star()`, and
  the consistency of `kappa_star()` with `r_max()`.
- At 90% power, the closed-form formula gives N = 556 instead of 558
  for the cell (0.3, 0.3, 1, 1, r = 1, rho = 0.5), where the exact
  continuous n2 is 278.0009. This boundary case is recorded in the
  tests.
- Tests reproduce values published by other authors: the 210
  tabulated values of C_2 in Sozu et al. (2015, Tables 4.3 and 4.4),
  by the exact root (all 210 to three decimals) and by the closed form
  (within 0.001 for rho <= 0.8 and 0.002 for rho = 0.95), and the
  sample size ratio n / m1 of Hung and Wang (2009) for two endpoints,
  which equals 1 / {1 - r_max(1, beta)}.

## Documentation

- The Description field and `README.md` no longer claim that the
  maximum reduction rate and the threshold effect-size ratio are
  unavailable from iterative methods.
- `plackett_gl_full()` is no longer marked as internal, since it is
  exported.

## Reproducibility

- The scripts under `inst/reproduce/` now load the installed package
  when the R source files are not found in the working directory, so
  they run both from the package and from the flat code supplement.
- The anchor weight curve for Figure 1 is now computed in
  `run_table_and_figure_manuscript.R`, so that
  `create_table_and_figure_manuscript.R` performs no numerical
  computation.
- `inst/reproduce/README.txt` describes the package layout.
- The numerical study now uses target powers 0.80, 0.85 and 0.90 (96
  cells in Table 2, which reports N_conv, N_prop and the achieved power
  for each target power and marks achieved powers below the target).
  Figures 1 and 2 use the same three powers, and Table 3 uses the
  tolerances 0.05, 0.10 and 0.15.
- Each table file written by `create_table_and_figure_manuscript.R`
  is now a complete table environment with caption, label and table
  notes (threeparttable), with the numbers in the captions taken from
  the cached results. The manuscript only needs `\input{}`.
- `tools/make_supplement.R` builds the flat code supplement of the
  manuscript from the package sources and checks that the copies are
  identical.
- `verify_intext_numbers.R` now checks the numbers cited in Sections
  2.4, 3 and 4 of the revised manuscript, including the 96-cell grid of
  Table 2 and the exact definition of the maximum reduction rate.
- `tools/sync_manuscript_assets.R` copies the generated tables and
  figures into the folder of the manuscript and checks that the copies
  are identical.

## Infrastructure

- Added `.gitattributes` to normalize line endings, and excluded
  `.Rhistory`, `.RData` and `.gitattributes` from the package build.
- The output folders of the reproduction scripts and `tools/` are
  excluded from the package build; the output folders are also
  excluded from git.

# twocpcont 0.1.0

Initial release.

## New features

- `twocpcont_ss()` computes the required sample size for a two-group
  superiority trial with two co-primary continuous endpoints under the
  known-variance assumption. Two methods are available: a sequential
  search using the bivariate normal distribution (`method = "bivariate"`,
  the conventional benchmark) and the proposed noniterative closed-form
  formula (`method = "univariate"`).
- `twocpcont_power()` computes the marginal and co-primary power for
  given sample sizes, with the same two methods.
- `r_max()` returns the maximum reduction rate achievable through
  correlation adjustment, in closed form.
- `kappa_star()` returns the threshold effect-size ratio at which the
  required sample size becomes insensitive to the between-endpoint
  correlation up to a user-specified tolerance.
- `plackett_gl_full()` evaluates the Plackett correlation integral and
  its first two derivatives via Gauss-Legendre quadrature.
- `GL_nodes_and_weights()` returns the Gauss-Legendre nodes and weights
  on the interval `[-1, 1]` for any positive integer number of nodes
  (default: 5).
- `print` methods are provided for objects of class `"twocpcont_ss"`
  and `"twocpcont_power"`.

## Reproducibility

- Scripts that reproduce the numerical results, tables, and figures of
  the accompanying manuscript are bundled under `inst/reproduce/`.
