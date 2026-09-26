#' Sample Size Calculation for Two Co-Primary Continuous Endpoints
#'
#' Calculates the required sample size for a two-group superiority trial
#' with two co-primary continuous endpoints under the known-variance
#' assumption.  Two methods are available: the conventional sequential
#' search using the bivariate normal distribution, and the proposed
#' noniterative closed-form formula based on an anchor weight derived
#' from a linear Taylor expansion of the marginal-power equation and a
#' correction term derived from a quadratic Taylor expansion of the
#' Plackett-approximated power function.
#'
#' @param delta1 Mean difference for the first endpoint.
#' @param delta2 Mean difference for the second endpoint.
#' @param sd1 Common standard deviation for the first endpoint.
#' @param sd2 Common standard deviation for the second endpoint.
#' @param rho Common correlation between the two endpoints.
#' @param r Allocation ratio \eqn{n_1 / n_2} (default: 1 for equal
#'   allocation).
#' @param alpha One-sided significance level (e.g., 0.025).
#' @param beta Target type II error rate (e.g., 0.2 for 80\% power).
#' @param method Character string specifying the calculation method:
#'   \code{"bivariate"} (default) uses a sequential search with the
#'   bivariate normal distribution (conventional benchmark);
#'   \code{"univariate"} uses the proposed noniterative closed-form
#'   formula.
#' @param gl_nodes Integer; number of Gauss-Legendre quadrature points
#'   used in the Plackett correlation integral when
#'   \code{method = "univariate"} (default: 5).  Ignored when
#'   \code{method = "bivariate"}.
#'
#' @return An object of class \code{"twocpcont_ss"} (which inherits from
#'   \code{"data.frame"}) with the following columns:
#'   \item{delta1}{Mean difference for endpoint 1.}
#'   \item{delta2}{Mean difference for endpoint 2.}
#'   \item{sd1}{Standard deviation for endpoint 1.}
#'   \item{sd2}{Standard deviation for endpoint 2.}
#'   \item{rho}{Correlation between endpoints.}
#'   \item{r}{Allocation ratio \eqn{n_1 / n_2}.}
#'   \item{alpha}{One-sided significance level.}
#'   \item{beta}{Type II error rate.}
#'   \item{method}{Calculation method used.}
#'   \item{n1}{Required sample size for group 1.}
#'   \item{n2}{Required sample size for group 2.}
#'   \item{N}{Total required sample size (\eqn{n_1 + n_2}).}
#'   \item{n2_cont}{Closed-form value of \eqn{n_2} before the ceiling
#'     operation, \eqn{c (\lambda_w^{*})^2}, when
#'     \code{method = "univariate"}; \code{NA} when
#'     \code{method = "bivariate"}.}
#'
#' @details
#' \strong{Notation.} Let \eqn{w} index the weaker endpoint
#' (\eqn{\gamma_w = \delta_w / \sigma_w \leq \delta_s / \sigma_s
#' = \gamma_s}) and \eqn{s} the stronger endpoint, with effect size
#' ratio \eqn{\kappa = \gamma_s / \gamma_w \geq 1}.  Let
#' \eqn{\lambda_w = \delta_w / (\sigma_w \sqrt{1 / n_1 + 1 / n_2})} and
#' \eqn{u_k = \lambda_k - z_\alpha} for \eqn{k} in \eqn{\{w, s\}}, so
#' that \eqn{u_s = \kappa \lambda_w - z_\alpha
#' = \kappa u_w + (\kappa - 1) z_\alpha}.
#'
#' \strong{Bivariate method (\code{method = "bivariate"}):}
#'
#' Uses a sequential search algorithm that increments or decrements
#' \eqn{n_2} until the minimum value satisfying the co-primary power
#' target is found.  The co-primary power is evaluated by
#' \code{\link{twocpcont_power}}.  This method serves as the benchmark
#' for the proposed closed-form formula.
#'
#' \strong{Noniterative closed-form formula (\code{method = "univariate"}):}
#'
#' Writing \eqn{p = \sqrt{1 - \beta}}, \eqn{u_p = \Phi^{-1}(p)}, and the
#' co-primary power as \eqn{P(\lambda_w) = \Phi(u_w) \Phi(u_s) +
#' I_{\mathrm{GL}}(u_w, u_s, \rho)}, the closed-form anchor weight is
#' \deqn{w^{*} = \frac{1}{2} +
#'   \frac{(\kappa - 1)(u_p + z_\alpha) \phi(u_p)}
#'        {(1 + \kappa) |2 \log p| p},}
#' clipped to \eqn{[0, 1]}.  The anchor non-centrality is
#' \eqn{\lambda_w^{(0)} = z_\alpha +
#' \Phi^{-1}((1 - \beta)^{w^{*}})}, and the correction term is
#' \deqn{\Delta^{*} = \frac{-P_1(\lambda_w^{(0)}) +
#'   \sqrt{P_1(\lambda_w^{(0)})^2 -
#'   2 P_2(\lambda_w^{(0)})
#'   \{P(\lambda_w^{(0)}) - (1 - \beta)\}}}{P_2(\lambda_w^{(0)})},}
#' where \eqn{P_1 = dP / d\lambda_w} and
#' \eqn{P_2 = d^2 P / d\lambda_w^2}.  Setting
#' \eqn{\lambda_w^{*} = \lambda_w^{(0)} + \Delta^{*}}, the final sample
#' size is \eqn{n_2 = \lceil c (\lambda_w^{*})^2 \rceil},
#' \eqn{n_1 = \lceil r n_2 \rceil}, \eqn{N = n_1 + n_2}, with
#' \eqn{c = \sigma_w^2 (1 + 1 / r) / \delta_w^2}.
#'
#' @examples
#' # Bivariate normal sequential search (conventional benchmark)
#' twocpcont_ss(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   method = "bivariate"
#' )
#'
#' # Noniterative closed-form formula with the default 5-point rule
#' twocpcont_ss(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   method = "univariate"
#' )
#'
#' # Noniterative closed-form formula with the 7-point rule
#' twocpcont_ss(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   method = "univariate", gl_nodes = 7
#' )
#'
#' @seealso \code{\link{twocpcont_power}}
#'
#' @importFrom stats qnorm pnorm dnorm
#'
#' @export
twocpcont_ss <- function(delta1, delta2, sd1, sd2, rho,
                         r = 1, alpha = 0.025, beta = 0.2,
                         method = c("bivariate", "univariate"),
                         gl_nodes = 5L) {

  method       <- match.arg(method)
  target_power <- 1 - beta

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

  # Scale factor c in eq. (eq:scale)
  c_scale <- sigma_w ^ 2 * (1 + 1 / r) / delta_w ^ 2

  if (method == "univariate") {

    z_alpha <- qnorm(1 - alpha)

    # ---- Anchor weight w_star (eq. (eq:w_star)) ----
    # w_star = 1/2 + (kappa - 1) (u_p + z_alpha) phi(u_p) /
    #               [(1 + kappa) |2 log p| p]
    p       <- sqrt(target_power)
    u_p     <- qnorm(p)
    phi_u_p <- dnorm(u_p)
    L_abs   <- abs(2 * log(p))   # |2 log p|, equal to -log(1 - beta)

    w_star <- 0.5 +
      (kappa - 1) * (u_p + z_alpha) * phi_u_p /
      ((1 + kappa) * L_abs * p)
    w_star <- min(1, max(0, w_star))

    # ---- Anchor noncentrality lambda_w^(0) (eq. (eq:lambda_anchor)) ----
    lam_w_0 <- z_alpha + qnorm(target_power ^ w_star)

    # ---- Quadratic-Taylor correction Delta_star (eq. (eq:Delta_star)) ----
    u_w_0 <- lam_w_0 - z_alpha
    u_s_0 <- kappa * lam_w_0 - z_alpha
    Phi_w <- pnorm(u_w_0)
    Phi_s <- pnorm(u_s_0)
    phi_w <- dnorm(u_w_0)
    phi_s <- dnorm(u_s_0)

    pk <- plackett_gl_full(u_w_0, u_s_0, rho, kappa, gl_nodes = gl_nodes)

    # P(lambda_w^(0)) - (1 - beta), P_1, P_2 at the anchor
    P_at <- Phi_w * Phi_s + pk$I_GL
    P_1  <- phi_w * Phi_s + kappa * phi_s * Phi_w + pk$I_GL_1
    P_2  <- -u_w_0 * phi_w * Phi_s +
            kappa * phi_w * phi_s -
            kappa ^ 2 * u_s_0 * phi_s * Phi_w +
            kappa * phi_s * phi_w + pk$I_GL_2

    C0 <- P_at - target_power

    # Positive square-root branch in eq. (eq:Delta_star).  At the
    # anchor, the Vieta product of the two roots is small while their
    # difference is of order |P_1| / |P_2|, so the positive branch
    # selects the root with smaller absolute value.
    Delta_star <- (-P_1 + sqrt(P_1 ^ 2 - 2 * P_2 * C0)) / P_2

    lam_w_star <- lam_w_0 + Delta_star

    n2_cont <- c_scale * lam_w_star ^ 2
    n2      <- ceiling(n2_cont)
    n1      <- ceiling(r * n2)
    N       <- n1 + n2

  } else {
    # ============================================================
    # method == "bivariate": sequential search benchmark
    # ============================================================

    z_alpha   <- qnorm(1 - alpha)
    z_beta_co <- qnorm(sqrt(target_power))

    n2_init <- ceiling(max(
      sd1 ^ 2 * (1 + 1 / r) / delta1 ^ 2 * (z_alpha + z_beta_co) ^ 2,
      sd2 ^ 2 * (1 + 1 / r) / delta2 ^ 2 * (z_alpha + z_beta_co) ^ 2
    ))

    power_fun <- function(n2_try) {
      n1_try <- ceiling(r * n2_try)
      twocpcont_power(
        n1 = n1_try, n2 = n2_try,
        delta1 = delta1, delta2 = delta2,
        sd1 = sd1, sd2 = sd2,
        rho = rho, alpha = alpha,
        method = "bivariate"
      )$powerCoprimary
    }

    n2 <- n2_init
    if (power_fun(n2) >= target_power) {
      while (n2 > 1L && power_fun(n2 - 1L) >= target_power) {
        n2 <- n2 - 1L
      }
    } else {
      while (power_fun(n2) < target_power) {
        n2 <- n2 + 1L
      }
    }

    n1      <- ceiling(r * n2)
    N       <- n1 + n2
    n2_cont <- NA_real_

  }

  # --- Assemble result ---
  result <- data.frame(
    delta1 = delta1,
    delta2 = delta2,
    sd1    = sd1,
    sd2    = sd2,
    rho    = rho,
    r      = r,
    alpha  = alpha,
    beta   = beta,
    method = method,
    n1     = n1,
    n2     = n2,
    N       = N,
    n2_cont = n2_cont
  )
  class(result) <- c("twocpcont_ss", "data.frame")

  return(result)
}
