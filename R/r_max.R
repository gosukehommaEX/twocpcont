#' Maximum Reduction Rate in Sample Size Due to Between-Endpoint Correlation
#'
#' Computes the maximum reduction rate \eqn{R_{\max}(\kappa, \beta)} that
#' the between-endpoint correlation can induce in the required sample
#' size for a two-group superiority trial with two co-primary continuous
#' endpoints under the known-variance assumption.  The reduction rate is
#' defined as
#' \deqn{R(\rho; \kappa, \beta) = 1 - N(\rho) / N(0)
#'   = 1 - (\lambda_w^{*}(\rho))^2 / (\lambda_w^{*}(0))^2,}
#' and its maximum is attained in the limit \eqn{\rho \to 1}.
#'
#' @param kappa Effect size ratio \eqn{\gamma_s / \gamma_w \geq 1}.
#' @param alpha One-sided significance level (default: 0.025).
#' @param beta Target type II error rate (default: 0.2 for 80\% power).
#'
#' @return Numeric scalar; the maximum reduction rate
#'   \eqn{R_{\max}(\kappa, \beta)}, lying in \eqn{[0, R_{\max}(1, \beta)]}.
#'
#' @details
#' \strong{Derivation.}  As \eqn{\rho \to 1}, the bivariate normal
#' cumulative distribution function satisfies
#' \eqn{\Phi_2(u_w, u_s; \rho) \to \Phi(\min(u_w, u_s)) = \Phi(u_w)}
#' since \eqn{\kappa \geq 1} implies \eqn{u_w \leq u_s}.  The power
#' equation therefore reduces to \eqn{\Phi(\lambda_w - z_\alpha)
#' = 1 - \beta}, giving the single-endpoint noncentrality parameter
#' \eqn{\lambda_w^{*}(\kappa, \beta, 1) = z_\alpha + z_\beta},
#' independent of \eqn{\kappa}.
#'
#' At \eqn{\rho = 0}, the Plackett correlation integral
#' \eqn{I_{\mathrm{GL}}} and both of its derivatives vanish, so the
#' quadratic Taylor correction coefficients of the closed-form formula
#' reduce to their marginal-only components.  Writing
#' \eqn{u_w^{(0)} = \lambda_w^{(0)} - z_\alpha},
#' \eqn{u_s^{(0)} = \kappa \lambda_w^{(0)} - z_\alpha} and
#' \eqn{\lambda_w^{(0)} = z_\alpha + \Phi^{-1}((1 - \beta)^{w^{*}})},
#' \deqn{P(\lambda_w^{(0)}) - (1 - \beta)
#'   = \Phi(u_w^{(0)}) \Phi(u_s^{(0)}) - (1 - \beta),}
#' \deqn{P_1(\lambda_w^{(0)}) = \phi(u_w^{(0)}) \Phi(u_s^{(0)})
#'         + \kappa \phi(u_s^{(0)}) \Phi(u_w^{(0)}),}
#' \deqn{P_2(\lambda_w^{(0)})
#'   = - u_w^{(0)} \phi(u_w^{(0)}) \Phi(u_s^{(0)})
#'   - \kappa^{2} u_s^{(0)} \phi(u_s^{(0)}) \Phi(u_w^{(0)})
#'   + 2 \kappa \phi(u_w^{(0)}) \phi(u_s^{(0)}),}
#' and the closed-form noncentrality parameter at \eqn{\rho = 0} is
#' \deqn{\lambda_w^{*}(\kappa, \beta, 0) = \lambda_w^{(0)}
#'   + \frac{- P_1 + \sqrt{P_1^{2} - 2 P_2 [P(\lambda_w^{(0)}) - (1 - \beta)]}}
#'          {P_2}.}
#' Combining the two limits yields the closed-form expression
#' \deqn{R_{\max}(\kappa, \beta)
#'   = 1 - \frac{(z_\alpha + z_\beta)^{2}}
#'         {[\lambda_w^{*}(\kappa, \beta, 0)]^{2}}.}
#' No iteration is required.
#'
#' \strong{At \eqn{\kappa = 1}.}  Since \eqn{w^{*} = 1 / 2} and the
#' linear Taylor equation is exact,
#' \eqn{P(\lambda_w^{(0)}) - (1 - \beta) = 0} and the quadratic
#' correction vanishes, giving
#' \eqn{\lambda_w^{*}(1, \beta, 0) = z_\alpha + \Phi^{-1}(\sqrt{1 - \beta})}
#' and
#' \deqn{R_{\max}(1, \beta) = 1 - \frac{(z_\alpha + z_\beta)^{2}}
#'   {[z_\alpha + \Phi^{-1}(\sqrt{1 - \beta})]^{2}}.}
#' This is also a universal upper bound on \eqn{R_{\max}(\kappa, \beta)}
#' over all \eqn{\kappa \geq 1}.
#'
#' @examples
#' # At kappa = 1 and 1 - beta = 0.80, R_max is about 23.8%.
#' r_max(kappa = 1.0, alpha = 0.025, beta = 0.2)
#'
#' # For kappa = 1.5 with the same power, R_max is about 2.3%.
#' r_max(kappa = 1.5, alpha = 0.025, beta = 0.2)
#'
#' @seealso \code{\link{twocpcont_ss}}, \code{\link{kappa_star}}
#'
#' @importFrom stats qnorm pnorm dnorm
#'
#' @export
r_max <- function(kappa, alpha = 0.025, beta = 0.2) {

  if (kappa < 1) {
    stop("kappa must be >= 1 (the weaker endpoint is indexed by w).")
  }

  target_power <- 1 - beta
  z_alpha <- qnorm(1 - alpha)
  z_beta  <- qnorm(target_power)

  # ---- Anchor weight w_star (eq. (eq:w_star)), same as in twocpcont_ss ----
  p       <- sqrt(target_power)
  u_p     <- qnorm(p)
  phi_u_p <- dnorm(u_p)
  L_abs   <- abs(2 * log(p))

  w_star <- 0.5 +
    (kappa - 1) * (u_p + z_alpha) * phi_u_p /
    ((1 + kappa) * L_abs * p)
  w_star <- min(1, max(0, w_star))

  # ---- Anchor noncentrality at rho = 0 ----
  lam_w_0 <- z_alpha + qnorm(target_power ^ w_star)

  # ---- Quadratic Taylor correction at rho = 0 ----
  # (I_GL and both derivatives vanish at rho = 0)
  u_w_0 <- lam_w_0 - z_alpha
  u_s_0 <- kappa * lam_w_0 - z_alpha
  Phi_w <- pnorm(u_w_0)
  Phi_s <- pnorm(u_s_0)
  phi_w <- dnorm(u_w_0)
  phi_s <- dnorm(u_s_0)

  C0  <- Phi_w * Phi_s - target_power
  P_1 <- phi_w * Phi_s + kappa * phi_s * Phi_w
  P_2 <- -u_w_0 * phi_w * Phi_s +
         kappa * phi_w * phi_s -
         kappa ^ 2 * u_s_0 * phi_s * Phi_w +
         kappa * phi_s * phi_w

  # Closed-form lambda_w*(rho = 0); at kappa = 1, C0 = 0 and Delta = 0
  if (abs(C0) < .Machine$double.eps * 16) {
    lam_w_at_rho0 <- lam_w_0
  } else {
    Delta_star    <- (-P_1 + sqrt(P_1 ^ 2 - 2 * P_2 * C0)) / P_2
    lam_w_at_rho0 <- lam_w_0 + Delta_star
  }

  # Maximum reduction rate
  1 - (z_alpha + z_beta) ^ 2 / lam_w_at_rho0 ^ 2
}
