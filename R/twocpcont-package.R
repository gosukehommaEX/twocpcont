#' twocpcont: Noniterative Closed-Form Sample Size for Trials with Two
#' Co-Primary Continuous Endpoints
#'
#' @description
#' The package computes the required sample size of a two-group superiority
#' trial with two co-primary continuous endpoints under the known-variance
#' assumption, either by the noniterative closed-form formula
#' (\code{\link{twocpcont_ss}} with \code{method = "univariate"}) or by the
#' conventional sequential search based on the bivariate normal
#' distribution (\code{method = "bivariate"}), and the constant of the
#' formula of Sugimoto et al. (2012) is computed by their Newton-Raphson
#' algorithm with \code{\link{newton_ck}} (\code{method = "newton"}).
#' The marginal and co-primary
#' powers for given sample sizes are computed by
#' \code{\link{twocpcont_power}}. The two design summaries, the maximum
#' reduction rate of the sample size attainable through the between-endpoint
#' correlation and the threshold effect-size ratio above which this reduction
#' is at most a given tolerance, are computed by \code{\link{r_max}} and
#' \code{\link{kappa_star}}. The Gauss-Legendre approximation of the
#' Plackett correlation integral that the closed form uses is available as
#' \code{\link{plackett_gl_full}} and \code{\link{GL_nodes_and_weights}}.
#'
#' @keywords internal
"_PACKAGE"
