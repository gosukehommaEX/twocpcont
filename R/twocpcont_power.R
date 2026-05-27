#' Power Calculation for Two Co-Primary Continuous Endpoints
#'
#' Calculates the power for a two-group superiority trial with two
#' co-primary continuous endpoints under the known-variance assumption.
#' Two methods are available: the conventional method based on the
#' bivariate normal distribution, and a univariate approximation based
#' on the Plackett identity with Gauss-Legendre quadrature.
#'
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param delta1 Mean difference for the first endpoint.
#' @param delta2 Mean difference for the second endpoint.
#' @param sd1 Common standard deviation for the first endpoint.
#' @param sd2 Common standard deviation for the second endpoint.
#' @param rho Common correlation between the two endpoints.
#' @param alpha One-sided significance level (e.g., 0.025).
#' @param method Character string specifying the calculation method:
#'   \code{"bivariate"} (default) uses the bivariate normal distribution
#'   (conventional benchmark); \code{"univariate"} uses the
#'   Plackett-identity approximation with Gauss-Legendre quadrature.
#' @param gl_nodes Integer; number of Gauss-Legendre quadrature points
#'   used when \code{method = "univariate"} (default: 5).  Ignored when
#'   \code{method = "bivariate"}.
#'
#' @return An object of class \code{"twocpcont_power"} (which inherits
#'   from \code{"data.frame"}) with the following columns:
#'   \item{n1}{Sample size for group 1.}
#'   \item{n2}{Sample size for group 2.}
#'   \item{delta1}{Mean difference for endpoint 1.}
#'   \item{delta2}{Mean difference for endpoint 2.}
#'   \item{sd1}{Standard deviation for endpoint 1.}
#'   \item{sd2}{Standard deviation for endpoint 2.}
#'   \item{rho}{Correlation between endpoints.}
#'   \item{alpha}{One-sided significance level.}
#'   \item{method}{Calculation method used.}
#'   \item{power1}{Marginal power for endpoint 1.}
#'   \item{power2}{Marginal power for endpoint 2.}
#'   \item{powerCoprimary}{Overall co-primary power.}
#'
#' @details
#' \strong{Notation.} Let \eqn{w} index the weaker endpoint
#' (\eqn{\gamma_w = \delta_w / \sigma_w \leq \delta_s / \sigma_s
#' = \gamma_s}) and \eqn{s} the stronger endpoint.  Let
#' \eqn{\lambda_k = \delta_k / (\sigma_k \sqrt{1 / n_1 + 1 / n_2})}
#' denote the non-centrality parameter for endpoint \eqn{k}, and
#' \eqn{u_k = \lambda_k - z_\alpha} its shifted counterpart.
#'
#' \strong{Bivariate method (\code{method = "bivariate"}):}
#'
#' The co-primary power is the bivariate normal CDF
#' \deqn{1 - \beta = \Phi_2(u_w, u_s; \rho),}
#' computed via \code{pbivnorm::pbivnorm}.
#'
#' \strong{Univariate approximation (\code{method = "univariate"}):}
#'
#' The co-primary power is computed via the Plackett identity combined
#' with Gauss-Legendre quadrature:
#' \deqn{1 - \beta = \Phi(u_w) \Phi(u_s)
#'   + I_{\mathrm{GL}}(u_w, u_s, \rho).}
#'
#' @examples
#' # Bivariate normal CDF (conventional benchmark)
#' twocpcont_power(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, alpha = 0.025,
#'   method = "bivariate"
#' )
#'
#' # Plackett-Gauss-Legendre approximation with the default 5-point rule
#' twocpcont_power(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, alpha = 0.025,
#'   method = "univariate"
#' )
#'
#' # Plackett-Gauss-Legendre approximation with the 7-point rule
#' twocpcont_power(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, alpha = 0.025,
#'   method = "univariate", gl_nodes = 7
#' )
#'
#' @seealso \code{\link{twocpcont_ss}}
#'
#' @importFrom stats qnorm pnorm
#' @importFrom pbivnorm pbivnorm
#'
#' @export
twocpcont_power <- function(n1, n2, delta1, delta2, sd1, sd2, rho,
                            alpha = 0.025,
                            method = c("bivariate", "univariate"),
                            gl_nodes = 5L) {

  method <- match.arg(method)

  # --- Identify weaker endpoint and define kappa ---
  gamma1 <- delta1 / sd1
  gamma2 <- delta2 / sd2
  if (gamma1 <= gamma2) {
    delta_w <- delta1; delta_s <- delta2
    sigma_w <- sd1;    sigma_s <- sd2
  } else {
    delta_w <- delta2; delta_s <- delta1
    sigma_w <- sd2;    sigma_s <- sd1
  }
  gamma_w <- delta_w / sigma_w
  gamma_s <- delta_s / sigma_s
  kappa   <- gamma_s / gamma_w

  z_alpha <- qnorm(1 - alpha)
  se_inv  <- 1 / sqrt(1 / n1 + 1 / n2)
  lam_w   <- gamma_w * se_inv
  lam_s   <- gamma_s * se_inv

  u_w <- lam_w - z_alpha
  u_s <- lam_s - z_alpha

  if (method == "bivariate") {
    power_co <- pbivnorm::pbivnorm(x = u_w, y = u_s, rho = rho)
    power1   <- pnorm(u_w)
    power2   <- pnorm(u_s)
  } else {
    pk     <- plackett_gl_full(u_w, u_s, rho, kappa, gl_nodes = gl_nodes)
    Phi_w  <- pnorm(u_w)
    Phi_s  <- pnorm(u_s)
    power_co <- Phi_w * Phi_s + pk$I_GL
    power1   <- Phi_w
    power2   <- Phi_s
  }

  # --- Reorder marginal powers back to the input endpoint order ---
  if (gamma1 <= gamma2) {
    power_e1 <- power1
    power_e2 <- power2
  } else {
    power_e1 <- power2
    power_e2 <- power1
  }

  result <- data.frame(
    n1             = n1,
    n2             = n2,
    delta1         = delta1,
    delta2         = delta2,
    sd1            = sd1,
    sd2            = sd2,
    rho            = rho,
    alpha          = alpha,
    method         = method,
    power1         = power_e1,
    power2         = power_e2,
    powerCoprimary = power_co
  )
  class(result) <- c("twocpcont_power", "data.frame")

  return(result)
}
