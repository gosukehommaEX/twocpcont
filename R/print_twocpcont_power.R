#' Print Method for twocpcont_power Objects
#'
#' Provides a clean, formatted display of power calculation results
#' from \code{\link{twocpcont_power}}.
#'
#' @param x An object of class \code{"twocpcont_power"}.
#' @param digits Integer; number of significant digits for numeric
#'   values (default: 4).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the original object \code{x}.
#'
#' @details
#' The output is divided into three sections: \code{[Design]},
#' \code{[Computation]}, and \code{[Results]}.  The method label
#' distinguishes the conventional bivariate approach from the proposed
#' univariate Plackett-GL approximation.
#'
#' @seealso \code{\link{twocpcont_power}}
#'
#' @examples
#' result <- twocpcont_power(
#'   n1 = 100, n2 = 100,
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, alpha = 0.025,
#'   method = "univariate"
#' )
#' print(result)
#'
#' @export
print.twocpcont_power <- function(x, digits = 4, ...) {

  # --- Helper: format a single numeric value ---
  fmt <- function(val, d = digits) {
    if (is.na(val)) return("NA")
    if (abs(val - round(val)) < 1e-10) return(sprintf("%.0f", val))
    formatC(signif(val, digits = d), format = "fg", flag = "#")
  }

  # --- Build label and value strings ---
  method_label <- if (x$method == "bivariate") {
    "Bivariate normal (conventional)"
  } else {
    "Plackett-GL approx. (proposed)"
  }

  # All labels that will appear in output
  all_labels <- c(
    "Sample size (n1, n2)",
    "Effect size (d1, d2)",
    "Std. deviation (s1, s2)",
    "Correlation (rho)",
    "Significance level",
    "Method",
    "Marginal power (ep. 1)",
    "Marginal power (ep. 2)",
    "Co-primary power"
  )

  # All value strings that will appear in output
  all_values <- c(
    paste(fmt(x$n1),     fmt(x$n2),     sep = ", "),
    paste(fmt(x$delta1), fmt(x$delta2), sep = ", "),
    paste(fmt(x$sd1),    fmt(x$sd2),    sep = ", "),
    fmt(x$rho),
    fmt(x$alpha),
    method_label,
    fmt(x$power1),
    fmt(x$power2),
    fmt(x$powerCoprimary)
  )

  # Column widths derived from actual content
  w <- max(nchar(all_labels))
  v <- max(nchar(all_values))

  # Separator: indent(2) + label(w) + " : "(3) + value(v)
  sep <- paste0("  ", strrep("-", w + 3L + v))

  # --- Helper: one aligned row ---
  row <- function(label, ...) {
    values <- paste(c(...), collapse = ", ")
    cat(sprintf("  %-*s : %s\n", w, label, values))
  }

  # --- Output ---
  cat("\n")
  cat("  Power calculation: two co-primary continuous endpoints\n")
  cat(sep, "\n")

  cat("  [Design]\n")
  row("Sample size (n1, n2)",    fmt(x$n1),     fmt(x$n2))
  row("Effect size (d1, d2)",    fmt(x$delta1), fmt(x$delta2))
  row("Std. deviation (s1, s2)", fmt(x$sd1),    fmt(x$sd2))
  row("Correlation (rho)",       fmt(x$rho))
  row("Significance level",      fmt(x$alpha))
  cat(sep, "\n")

  cat("  [Computation]\n")
  row("Method", method_label)
  cat(sep, "\n")

  cat("  [Results]\n")
  row("Marginal power (ep. 1)", fmt(x$power1))
  row("Marginal power (ep. 2)", fmt(x$power2))
  row("Co-primary power",       fmt(x$powerCoprimary))
  cat(sep, "\n")
  cat("\n")

  invisible(x)
}
