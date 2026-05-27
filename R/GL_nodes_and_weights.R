#' Gauss-Legendre Quadrature Nodes and Weights
#'
#' Computes the nodes and weights for \code{gl_nodes}-point Gauss-Legendre
#' quadrature on the interval \eqn{[-1, 1]} using the Golub-Welsch
#' algorithm, which reduces the problem to a symmetric tridiagonal
#' eigenvalue decomposition.
#'
#' @param gl_nodes Integer; the number of quadrature points (default: 5).
#'   Must be a positive integer.  Choosing a larger value increases the
#'   accuracy of the Plackett correlation integral approximation used in
#'   the sample size formula at a modest computational cost.
#'
#' @return A list with two components:
#'   \item{x}{Numeric vector of length \code{gl_nodes} containing the
#'     quadrature nodes (zeros of the \code{gl_nodes}-th Legendre
#'     polynomial) in increasing order.}
#'   \item{w}{Numeric vector of length \code{gl_nodes} containing the
#'     corresponding quadrature weights.}
#'
#' @details
#' The \code{gl_nodes}-point Gauss-Legendre rule integrates polynomials
#' of degree up to \eqn{2 \mathrm{gl\_nodes} - 1} exactly:
#' \deqn{\int_{-1}^{1} f(x) dx \approx
#'   \sum_{i = 1}^{\mathrm{gl\_nodes}} w_i f(x_i).}
#' The nodes \eqn{x_1, \ldots, x_{\mathrm{gl\_nodes}}} are the
#' eigenvalues of the
#' \eqn{\mathrm{gl\_nodes} \times \mathrm{gl\_nodes}} symmetric
#' tridiagonal matrix with zero diagonal and subdiagonal entries
#' \eqn{b_j = j / \sqrt{4 j^2 - 1}},
#' \eqn{j = 1, \ldots, \mathrm{gl\_nodes} - 1}.  The weights are
#' \eqn{w_i = 2 v_{i, 1}^2}, where \eqn{v_{i, 1}} is the first component
#' of the \eqn{i}-th eigenvector.  The notation \eqn{(x_i, w_i)} matches
#' that used in the proposed sample size formula and in the Plackett
#' integral approximation \eqn{I_{\mathrm{GL}}}.
#'
#' @examples
#' # Default 5-point rule
#' gl5 <- GL_nodes_and_weights()
#' gl5$x
#' gl5$w
#'
#' # 7-point rule
#' gl7 <- GL_nodes_and_weights(7)
#' gl7$x
#' gl7$w
#'
#' # Verify: integrate x^8 over [-1, 1] (exact = 2 / 9)
#' gl <- GL_nodes_and_weights(5)
#' sum(gl$w * gl$x ^ 8)
#'
#' @export
GL_nodes_and_weights <- function(gl_nodes = 5L) {

  gl_nodes <- as.integer(gl_nodes)
  if (gl_nodes < 1L) stop("gl_nodes must be a positive integer.")

  if (gl_nodes == 1L) {
    return(list(x = 0, w = 2))
  }

  # Golub-Welsch algorithm: construct the symmetric tridiagonal Jacobi
  # matrix for the Legendre polynomials.  The diagonal is zero; the
  # subdiagonal entries are b_j = j / sqrt(4 j^2 - 1),
  # j = 1, ..., gl_nodes - 1.
  j    <- seq_len(gl_nodes - 1L)
  beta <- j / sqrt(4 * j ^ 2 - 1)

  # Build the tridiagonal matrix
  J <- matrix(0, nrow = gl_nodes, ncol = gl_nodes)
  for (i in seq_len(gl_nodes - 1L)) {
    J[i, i + 1L] <- beta[i]
    J[i + 1L, i] <- beta[i]
  }

  # Eigenvalue decomposition
  eig <- eigen(J, symmetric = TRUE)

  # Nodes: eigenvalues (sorted in increasing order)
  ord <- order(eig$values)
  x   <- eig$values[ord]
  # Weights: 2 * (first component of each eigenvector)^2
  w   <- 2 * eig$vectors[1, ord] ^ 2

  list(x = x, w = w)
}
