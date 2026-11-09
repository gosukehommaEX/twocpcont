# Reproduction of published results

One script per published article whose tables or figures are reproduced with
the functions of the package: `reproduce_<Author><Year>.R`. Each script
writes `results/reproduce_<Author><Year>.csv` with one row per published
value (published value, value of the package, difference, verdict PASS,
EXPLAINED, FAIL or INFO). A value that cannot be reproduced is EXPLAINED only
with a reason registered in `../validation/explained.csv`.

| Script | Published values |
|---|---|
| `reproduce_Sozu2015.R` | The 210 values of the constant C_2 in Sozu et al. (2015, Tables 4.3 and 4.4), stored in `inst/extdata/sozu2015_C2_table.csv`, by the exact root of the power equation and by the closed form |
| `reproduce_Sugimoto2012.R` | The constant C = 0.9988124 and the sample size n = 142.93 of the example of Sugimoto et al. (2012, Supporting Information, Appendix A.2.1) by their Newton-Raphson algorithm (`newton_ck()`), and the powers of the sequential search in their Table I |
| `reproduce_HungWang2009.R` | The sample size ratio n / m1 of Hung and Wang (2009) for two endpoints, which equals 1 / {1 - R_max(1, beta)} |

Run from the package root after `devtools::load_all()`, e.g.
`source("inst/reproduce/reproduce_Sozu2015.R")`.
