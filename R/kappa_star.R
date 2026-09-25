#' Threshold Effect Size Ratio Beyond Which Correlation Has Negligible Effect
#'
#' Computes the threshold effect size ratio
#' \eqn{\kappa^{*}(\epsilon, \beta)} above which the between-endpoint
#' correlation reduces the required sample size by at most a fraction
#' \eqn{\epsilon} regardless of the actual value of \eqn{\rho}.  In
#' practical terms, when \eqn{\kappa \geq \kappa^{*}(\epsilon, \beta)},
#' the iterative correlation adjustment can be skipped without losing
#' more than a fraction \eqn{\epsilon} of sample size efficiency.
#'
#' @param epsilon Tolerance for the maximum reduction rate, in
#'   \eqn{(0, 1)} (e.g., \code{0.01} for a 1\% tolerance).
#' @param alpha One-sided significance level (default: 0.025).
#' @param beta Target type II error rate (default: 0.2 for 80\% power).
#'
#' @return Numeric scalar; the threshold
#'   \eqn{\kappa^{*}(\epsilon, \beta) \geq 1} such that
#'   \eqn{R_{\max}(\kappa, \beta) \leq \epsilon} for all
#'   \eqn{\kappa \geq \kappa^{*}}.
#'
#' @details
#' \strong{Derivation.}  Using the exact \eqn{\rho = 0} factorization
#' \eqn{\Phi(v) \Phi(\kappa v + (\kappa - 1) z_\alpha) = 1 - \beta} and
#' the limit \eqn{\lambda_w^{*}(\kappa, \beta, 1) = z_\alpha + z_\beta},
#' the equation \eqn{R_{\max}(\kappa^{*}, \beta) = \epsilon} can be
#' inverted in closed form.  Setting
#' \eqn{v^{*} = z_\beta + \nu} with
#' \eqn{\nu = (z_\alpha + z_\beta) (1 / \sqrt{1 - \epsilon} - 1)} and
#' \eqn{G = \kappa^{*} (z_\alpha + z_\beta + \nu) - z_\alpha}, the
#' equation becomes
#' \deqn{\Phi(z_\beta + \nu) \Phi(G) = 1 - \beta,}
#' which gives \eqn{G = \Phi^{-1}((1 - \beta) / \Phi(z_\beta + \nu))} and
#' hence
#' \deqn{\kappa^{*}(\epsilon, \beta)
#'   = \frac{z_\alpha + \Phi^{-1}((1 - \beta) / \Phi(z_\beta + \nu))}
#'          {z_\alpha + z_\beta + \nu}.}
#' This is a fully closed-form expression involving only the standard
#' normal cumulative distribution function and its inverse; no iteration
#' is required.  Because the same exact \eqn{\rho = 0} relation defines
#' \code{\link{r_max}} with \code{exact = TRUE} (the default),
#' \code{r_max(kappa_star(epsilon, alpha, beta), alpha, beta)} returns
#' \eqn{\epsilon} up to rounding error.
#'
#' \strong{Interpretation.}  For typical design choices, the threshold
#' is moderate.  For example, at \eqn{1 - \beta = 0.80} and
#' \eqn{\epsilon = 0.01}, the threshold is approximately
#' \eqn{\kappa^{*} \approx 1.61}.  In trials where the effect size
#' ratio exceeds this threshold, the sample size is essentially
#' determined by the weaker endpoint and depends only weakly on
#' \eqn{\rho}.
#'
#' @examples
#' # Threshold at 1% tolerance, 80% power: kappa_star is about 1.61.
#' kappa_star(epsilon = 0.01, alpha = 0.025, beta = 0.2)
#'
#' # Threshold at 5% tolerance, 90% power: kappa_star is about 1.24.
#' kappa_star(epsilon = 0.05, alpha = 0.025, beta = 0.1)
#'
#' @seealso \code{\link{r_max}}, \code{\link{twocpcont_ss}}
#'
#' @importFrom stats qnorm pnorm
#'
#' @export
kappa_star <- function(epsilon, alpha = 0.025, beta = 0.2) {

  if (epsilon <= 0 || epsilon >= 1) {
    stop("epsilon must lie in (0, 1).")
  }

  z_alpha <- qnorm(1 - alpha)
  z_beta  <- qnorm(1 - beta)

  nu <- (z_alpha + z_beta) * (1 / sqrt(1 - epsilon) - 1)
  G  <- qnorm((1 - beta) / pnorm(z_beta + nu))

  (z_alpha + G) / (z_alpha + z_beta + nu)
}
