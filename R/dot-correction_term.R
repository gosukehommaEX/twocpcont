#' Correction term of the closed-form formula
#'
#' Root with the smaller absolute value of the quadratic Taylor equation
#' \eqn{C_0 + P_1 \Delta + P_2 \Delta^2 / 2 = 0} at the anchor, computed in
#' the form \eqn{\Delta^{*} = -2 C_0 / (P_1 + \sqrt{P_1^2 - 2 P_2 C_0})},
#' which equals \eqn{(-P_1 + \sqrt{P_1^2 - 2 P_2 C_0}) / P_2} and remains
#' defined when \eqn{P_2} is close to zero. If the discriminant is
#' negative, the quadratic equation has no real root, and the first-order
#' correction \eqn{-C_0 / P_1} (one Newton step from the anchor) is
#' returned with a warning.
#'
#' @param C0 \eqn{P(\lambda_w^{(0)}) - (1 - \beta)} at the anchor.
#' @param P_1 First derivative of \eqn{P} with respect to \eqn{\lambda_w}
#'   at the anchor.
#' @param P_2 Second derivative of \eqn{P} with respect to
#'   \eqn{\lambda_w} at the anchor.
#'
#' @return Numeric scalar; the correction \eqn{\Delta^{*}}.
#'
#' @noRd
.correction_term <- function(C0, P_1, P_2) {
  disc <- P_1 ^ 2 - 2 * P_2 * C0
  if (!is.na(disc) && disc < 0) {
    warning("The discriminant of the quadratic correction is negative; ",
            "the first-order correction -(P - (1 - beta)) / P_1 is used.",
            call. = FALSE)
    return(-C0 / P_1)
  }
  -2 * C0 / (P_1 + sqrt(disc))
}
