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
  on the interval `[0, rho]` for any positive integer number of nodes
  (default: 5).
- `print` methods are provided for objects of class `"twocpcont_ss"`
  and `"twocpcont_power"`.

## Reproducibility

- Scripts that reproduce the numerical results, tables, and figures of
  the accompanying manuscript are bundled under `inst/reproduce/`.
