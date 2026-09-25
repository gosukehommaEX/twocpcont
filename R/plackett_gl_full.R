#' Plackett Integral and Analytic Derivatives via Gauss-Legendre Quadrature
#'
#' Computes the Plackett correlation integral
#' \eqn{I_{\mathrm{GL}}(u_w, u_s, \rho)
#'   = (\rho / 2) \sum_{i = 1}^{\mathrm{gl\_nodes}}
#'     w_i \phi_2(u_w, u_s; t_i)},
#' which approximates
#' \eqn{\int_0^\rho \phi_2(u_w, u_s; t) dt}, together with its first
#' and second derivatives with respect to \eqn{\lambda_w}, where
#' \eqn{u_w = \lambda_w - z_\alpha} and
#' \eqn{u_s = \kappa \lambda_w - z_\alpha}.  For interior \eqn{\rho}
#' (\eqn{|\rho| < 1}), all three quantities are evaluated using the
#' same \code{gl_nodes}-point Gauss-Legendre quadrature on
#' \eqn{[0, \rho]}.  At the boundary values \eqn{\rho = 0} and
#' \eqn{|\rho| = 1}, the analytic limits are returned directly without
#' quadrature.
#'
#' @param u_w Standardized argument for the weaker endpoint.
#' @param u_s Standardized argument for the stronger endpoint.
#' @param rho Correlation between the two endpoints, in \eqn{[-1, 1]}.
#' @param kappa Effect size ratio \eqn{\gamma_s / \gamma_w \geq 1}.
#' @param gl_nodes Integer; number of Gauss-Legendre quadrature points
#'   used to approximate the Plackett integral (default: 5).
#'
#' @return A list with three components:
#'   \item{I_GL}{Approximation of the Plackett correlation integral
#'     \eqn{I_{\mathrm{GL}}(u_w, u_s, \rho)}.}
#'   \item{I_GL_1}{First derivative of the Plackett integral with
#'     respect to \eqn{\lambda_w}.}
#'   \item{I_GL_2}{Second derivative of the Plackett integral with
#'     respect to \eqn{\lambda_w}.}
#'
#' @details
#' Let \eqn{\phi_2(u_w, u_s; t)} denote the standard bivariate normal
#' density with correlation \eqn{t}.  For interior \eqn{\rho}, the
#' standard \code{gl_nodes}-point Gauss-Legendre nodes \eqn{x_i} and
#' weights \eqn{w_i} on \eqn{[-1, 1]} from
#' \code{\link{GL_nodes_and_weights}} are mapped to \eqn{[0, \rho]} by
#' \eqn{t_i = \rho (x_i + 1) / 2} with Jacobian \eqn{\rho / 2}.  The
#' integral and derivatives are
#' \deqn{I_{\mathrm{GL}} = \frac{\rho}{2}
#'   \sum_{i = 1}^{\mathrm{gl\_nodes}} w_i \phi_2(u_w, u_s; t_i),}
#' \deqn{I_{\mathrm{GL}, 1} = -\frac{\rho}{2}
#'   \sum_{i = 1}^{\mathrm{gl\_nodes}} w_i \phi_2(u_w, u_s; t_i)
#'   \frac{(1 - \kappa t_i) u_w + (\kappa - t_i) u_s}{1 - t_i^2},}
#' \deqn{I_{\mathrm{GL}, 2} = \frac{\rho}{2}
#'   \sum_{i = 1}^{\mathrm{gl\_nodes}}
#'   \frac{w_i \phi_2(u_w, u_s; t_i)}{1 - t_i^2}
#'   \left\{\frac{[(1 - \kappa t_i) u_w + (\kappa - t_i) u_s]^2}
#'              {1 - t_i^2}
#'        - (1 + \kappa^2 - 2 \kappa t_i)\right\}.}
#'
#' At \eqn{\rho = 0}, the integral and both derivatives vanish.  At
#' \eqn{\rho = 1} with \eqn{\kappa \geq 1} (so \eqn{u_w \leq u_s}), the
#' bivariate normal CDF satisfies \eqn{\Phi_2(u_w, u_s; 1) =
#' \Phi(u_w)}, hence \eqn{I_{\mathrm{GL}}(u_w, u_s, 1) =
#' \Phi(u_w)(1 - \Phi(u_s))} and differentiation with respect to
#' \eqn{\lambda_w} gives
#' \deqn{I_{\mathrm{GL}, 1}\bigg|_{\rho = 1}
#'   = \phi(u_w)(1 - \Phi(u_s)) - \kappa \Phi(u_w) \phi(u_s),}
#' \deqn{I_{\mathrm{GL}, 2}\bigg|_{\rho = 1}
#'   = -u_w \phi(u_w)(1 - \Phi(u_s)) - 2 \kappa \phi(u_w) \phi(u_s)
#'   + \kappa^2 u_s \Phi(u_w) \phi(u_s).}
#' At \eqn{\rho = -1} with \eqn{\kappa \geq 1}, \eqn{\Phi_2(u_w, u_s;
#' -1) = \max(0, \Phi(u_w) + \Phi(u_s) - 1)} and the corresponding
#' analytic limits are returned.  The Gauss-Legendre rule loses
#' accuracy as \eqn{|\rho|} approaches unity because the integrand
#' \eqn{1 / \sqrt{1 - t^2}} becomes unbounded at the endpoint
#' \eqn{t = \pm 1}.
#'
#' @examples
#' # Default 5-point rule at a sample anchor
#' plackett_gl_full(u_w = 0.84, u_s = 0.84, rho = 0.3, kappa = 1)
#'
#' # 7-point rule for higher accuracy
#' plackett_gl_full(u_w = 0.84, u_s = 0.84, rho = 0.3, kappa = 1,
#'                  gl_nodes = 7)
#'
#' # Boundary handling at rho = 1
#' plackett_gl_full(u_w = 1.25, u_s = 1.25, rho = 1, kappa = 1)
#'
#' @export
plackett_gl_full <- function(u_w, u_s, rho, kappa, gl_nodes = 5L) {

  if (rho == 0) {
    return(list(I_GL = 0, I_GL_1 = 0, I_GL_2 = 0))
  }

  # Analytic limit at rho = 1 (with kappa >= 1 implying u_w <= u_s)
  if (rho == 1) {
    Phi_w <- pnorm(u_w)
    Phi_s <- pnorm(u_s)
    phi_w <- dnorm(u_w)
    phi_s <- dnorm(u_s)
    I_GL_lim   <- Phi_w * (1 - Phi_s)
    I_GL_1_lim <- phi_w * (1 - Phi_s) - kappa * Phi_w * phi_s
    I_GL_2_lim <- -u_w * phi_w * (1 - Phi_s) -
      2 * kappa * phi_w * phi_s +
      kappa ^ 2 * u_s * Phi_w * phi_s
    return(list(I_GL = I_GL_lim, I_GL_1 = I_GL_1_lim, I_GL_2 = I_GL_2_lim))
  }

  # Analytic limit at rho = -1 (with kappa >= 1 implying u_w <= u_s)
  if (rho == -1) {
    Phi_w <- pnorm(u_w)
    Phi_s <- pnorm(u_s)
    phi_w <- dnorm(u_w)
    phi_s <- dnorm(u_s)
    if (Phi_w + Phi_s <= 1) {
      # Phi_2 = 0; I_GL = -Phi(u_w) Phi(u_s)
      I_GL_lim   <- -Phi_w * Phi_s
      I_GL_1_lim <- -phi_w * Phi_s - kappa * Phi_w * phi_s
      I_GL_2_lim <- u_w * phi_w * Phi_s - 2 * kappa * phi_w * phi_s +
        kappa ^ 2 * u_s * Phi_w * phi_s
    } else {
      # Phi_2 = Phi(u_w) + Phi(u_s) - 1
      I_GL_lim   <- Phi_w + Phi_s - 1 - Phi_w * Phi_s
      I_GL_1_lim <- phi_w + kappa * phi_s -
        phi_w * Phi_s - kappa * Phi_w * phi_s
      I_GL_2_lim <- -u_w * phi_w - kappa ^ 2 * u_s * phi_s +
        u_w * phi_w * Phi_s - 2 * kappa * phi_w * phi_s +
        kappa ^ 2 * u_s * Phi_w * phi_s
    }
    return(list(I_GL = I_GL_lim, I_GL_1 = I_GL_1_lim, I_GL_2 = I_GL_2_lim))
  }

  gl       <- GL_nodes_and_weights(gl_nodes = gl_nodes)
  half_rho <- rho / 2

  # Map [-1, 1] Gauss-Legendre nodes x_i to integration variable t_i
  # on [0, rho]: t_i = (rho / 2) (x_i + 1)
  t_pts <- half_rho * (gl$x + 1)
  w_pts <- gl$w

  # Standard bivariate normal density phi_2(u_w, u_s; t_i) at each node
  one_m_t2  <- 1 - t_pts ^ 2
  quad_form <- (u_w ^ 2 - 2 * t_pts * u_w * u_s + u_s ^ 2) / one_m_t2
  phi_2     <- exp(-0.5 * quad_form) / (2 * pi * sqrt(one_m_t2))

  num     <- (1 - kappa * t_pts) * u_w + (kappa - t_pts) * u_s
  bracket <- num ^ 2 / one_m_t2 - (1 + kappa ^ 2 - 2 * kappa * t_pts)

  I_GL   <- half_rho * sum(w_pts * phi_2)
  I_GL_1 <- -half_rho * sum(w_pts * phi_2 * num / one_m_t2)
  I_GL_2 <- half_rho * sum(w_pts * phi_2 * bracket / one_m_t2)

  list(I_GL = I_GL, I_GL_1 = I_GL_1, I_GL_2 = I_GL_2)
}
