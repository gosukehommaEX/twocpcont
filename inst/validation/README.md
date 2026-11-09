# Validation by independent routes

One script per comparison of independent routes (closed form against
simulation, the package against another package, the package against an
independent implementation): `validate_<content>.R`. Each script writes
`results/validate_<content>.csv` with the verdict of each comparison (PASS,
EXPLAINED, FAIL or INFO). The reasons of the EXPLAINED verdicts are
registered in `explained.csv`; the summary of the verdicts is written by the
scripts, never by hand.

The statements of the article that are not single numbers (bounds,
monotonicity, which design combinations meet the target power) are checked by
the scripts below; the numbers themselves are checked against the ledger by
`inst/scripts/number_check.R`.

| Script | Statements checked |
|---|---|
| `validate_newton_mvtnorm.R` | `newton_ck()` against the R function `CKsolution` of Sugimoto et al. (2012), transcribed from their Supporting Information and run with the `mvtnorm` package, on the 210 designs of Sozu et al. (2015) and the 48 designs of Table 2 |
| `validate_claims_methods.R` | Section 2 and Section S2: accuracy of the 5-point Gauss-Legendre rule against adaptive quadrature, the anchor weight against the exact root of (S2), and the properties of `r_max()` and `kappa_star()` (closed form against `uniroot()`) |
| `validate_claims_numerical.R` | Section 3: the combinations of Table 2 that meet the target power, monotonicity in the correlation, the direction of the conventional search, the sign of the discriminant of the correction term, and the deviations from the published tables |
| `validate_claims_example.R` | Section 4, Section 5 and Section S5: Table 4, the thresholds of the example, the protocol's table of marginal powers, and the power checks |
| `validate_error_decomposition.R` | Figures S1 and S2 and Table S2: the analytic P, P_1 and P_2 at the anchor against numerical derivatives, the steps from the anchor, the agreement of the quadratic correction with `twocpcont_ss()`, the roots of the power equations, and the sign of the discriminant of the correction term |

Run from the package root after `devtools::load_all()` and
`inst/scripts/data_generate_main.R` (and, for
`validate_error_decomposition.R`, `inst/sup_info/data_generate_supplement.R`),
e.g. `source("inst/validation/validate_claims_methods.R")`.
