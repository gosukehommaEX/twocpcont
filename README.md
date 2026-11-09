# twocpcont

<!-- badges: start -->
[![License: GPL (>= 2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
<!-- badges: end -->

`twocpcont` provides a noniterative, closed-form sample size formula
for two-group superiority clinical trials with two co-primary
continuous endpoints under the known-variance assumption.

The required total sample size is expressed as an explicit function of
the significance level, the target power, the standardized effect
sizes, the between-endpoint correlation, and the allocation ratio.
The bivariate normal cumulative distribution function that appears in
the co-primary power is rendered tractable by combining Plackett's
identity with a 5-point Gauss-Legendre quadrature for the correlation
integral, and the resulting univariate power equation is inverted in
closed form by two analytic Taylor linearizations. The package also
implements the conventional sequential-search method based on the
bivariate normal distribution for benchmarking purposes.

The package also provides two design summaries:

- the maximum reduction rate, an upper bound on the fractional sample
  size reduction attainable through correlation adjustment, and
- a threshold effect-size ratio above which the between-endpoint
  correlation changes the required sample size by at most a
  user-specified fraction.

## Installation

You can install the development version of `twocpcont` from
[GitHub](https://github.com/gosukehommaEX/twocpcont) with:

``` r
# install.packages("remotes")
remotes::install_github("gosukehommaEX/twocpcont", build_vignettes = TRUE)
```

## Main functions

| Function | Description |
|---|---|
| `twocpcont_ss()` | Required sample size by the proposed closed-form formula (`method = "univariate"`) or by the conventional sequential search (`method = "bivariate"`) |
| `twocpcont_power()` | Marginal and co-primary power for given sample sizes |
| `r_max()` | Maximum reduction rate of the sample size attainable through the correlation |
| `kappa_star()` | Threshold effect-size ratio above which the reduction is at most a given tolerance |
| `plackett_gl_full()` | Plackett correlation integral and its first two derivatives by Gauss-Legendre quadrature |
| `GL_nodes_and_weights()` | Gauss-Legendre nodes and weights on `[-1, 1]` |

## Example

### Sample size calculation

The proposed noniterative closed-form formula:

``` r
library(twocpcont)

twocpcont_ss(
  delta1 = 1, delta2 = 1,
  sd1 = 2, sd2 = 2,
  rho = 0.5,
  r = 1,
  alpha = 0.025, beta = 0.2,
  method = "univariate"
)
```

The conventional sequential search using the bivariate normal
distribution (for benchmarking):

``` r
twocpcont_ss(
  delta1 = 1, delta2 = 1,
  sd1 = 2, sd2 = 2,
  rho = 0.5,
  r = 1,
  alpha = 0.025, beta = 0.2,
  method = "bivariate"
)
```

### Power calculation

``` r
twocpcont_power(
  n1 = 100, n2 = 100,
  delta1 = 1, delta2 = 1,
  sd1 = 2, sd2 = 2,
  rho = 0.5,
  alpha = 0.025,
  method = "univariate"
)
```

### Maximum reduction rate

``` r
r_max(kappa = 1, beta = 0.2)
```

### Threshold effect-size ratio

``` r
kappa_star(epsilon = 0.05, beta = 0.2)
```

## Reproducing the manuscript results

The scripts that reproduce all numerical results, tables, and figures
of the accompanying manuscript are in `inst/scripts/` (manuscript) and
`inst/sup_info/` (Supporting Information) of the source package, with
the reproduction of published values in `inst/reproduce/` and the
checks of the statements of the manuscript in `inst/validation/`. Clone
the repository, open `twocpcont.Rproj`, and run the scripts from the
package root in the order given in `inst/scripts/README.md`, starting
with:

``` r
devtools::load_all()
source("inst/scripts/data_generate_main.R")
source("inst/scripts/table_and_figure_manuscript.R")
```

## Author

Gosuke Homma
(<my.name.is.gosuke@gmail.com>)

## License

GPL (>= 2)
