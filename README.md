# twocpcont

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
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

Beyond the sample size itself, the closed-form structure yields two
analytic design quantities that are not directly available from the
conventional iterative method:

- a maximum reduction rate giving a closed-form upper bound on the
  fractional sample size reduction attainable through correlation
  adjustment, and
- a threshold effect-size ratio identifying design configurations in
  which the precise value of the between-endpoint correlation has
  negligible impact on the required sample size.

## Installation

You can install the development version of `twocpcont` from
[GitHub](https://github.com/gosukehommaEX/twocpcont) with:

``` r
# install.packages("remotes")
remotes::install_github("gosukehommaEX/twocpcont")
```

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
in the accompanying manuscript are bundled with the package under
`inst/reproduce/`. After installation, the path can be obtained by:

``` r
system.file("reproduce", package = "twocpcont")
```

A `README.txt` in that directory explains the execution order.

## Author

Gosuke Homma
(<my.name.is.gosuke@gmail.com>)

## License

MIT (c) Gosuke Homma. See the `LICENSE` and `LICENSE.md` files for
details.
