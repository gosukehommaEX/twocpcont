#' Print Method for twocpcont_ss Objects
#'
#' Provides a clean, formatted display of sample size calculation
#' results from \code{\link{twocpcont_ss}}.
#'
#' @param x An object of class \code{"twocpcont_ss"}.
#' @param digits Integer; number of significant digits for numeric
#'   values (default: 4).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the original object \code{x}.
#'
#' @details
#' The output is divided into three sections: \code{[Design]},
#' \code{[Computation]}, and \code{[Results]}.
#'
#' @seealso \code{\link{twocpcont_ss}}
#'
#' @examples
#' result <- twocpcont_ss(
#'   delta1 = 0.5, delta2 = 0.5,
#'   sd1 = 1, sd2 = 1,
#'   rho = 0.3, r = 1,
#'   alpha = 0.025, beta = 0.2,
#'   method = "univariate"
#' )
#' print(result)
#'
#' @export
print.twocpcont_ss <- function(x, digits = 4, ...) {

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
    "Noniterative closed-form (proposed)"
  }

  # All labels that will appear in output
  all_labels <- c(
    "Mean difference (d1, d2)",
    "Std. deviation (s1, s2)",
    "Correlation (rho)",
    "Allocation ratio (r)",
    "Significance level",
    "Target power",
    "Method",
    "Sample size (n1, n2)",
    "Total sample size (N)"
  )

  # All value strings that will appear in output
  all_values <- c(
    paste(fmt(x$delta1), fmt(x$delta2), sep = ", "),
    paste(fmt(x$sd1),    fmt(x$sd2),    sep = ", "),
    fmt(x$rho),
    fmt(x$r),
    fmt(x$alpha),
    fmt(1 - x$beta),
    method_label,
    paste(fmt(x$n1), fmt(x$n2), sep = ", "),
    fmt(x$N)
  )

  # Column widths derived from actual content
  w <- max(nchar(all_labels))
  v <- max(nchar(all_values))

  # Header title
  title <- "Sample size calculation: two co-primary continuous endpoints"

  # Separator width: max of row width and title length
  line_width <- max(w + 3L + v, nchar(title))
  sep <- paste0("  ", strrep("-", line_width))

  # --- Helper: one aligned row ---
  row <- function(label, ...) {
    values <- paste(c(...), collapse = ", ")
    cat(sprintf("  %-*s : %s\n", w, label, values))
  }

  # --- Output ---
  cat("\n")
  cat(sprintf("  %s\n", title))
  cat(sep, "\n")

  cat("  [Design]\n")
  row("Mean difference (d1, d2)", fmt(x$delta1), fmt(x$delta2))
  row("Std. deviation (s1, s2)", fmt(x$sd1),    fmt(x$sd2))
  row("Correlation (rho)",       fmt(x$rho))
  row("Allocation ratio (r)",    fmt(x$r))
  row("Significance level",      fmt(x$alpha))
  row("Target power",            fmt(1 - x$beta))
  cat(sep, "\n")

  cat("  [Computation]\n")
  row("Method", method_label)
  cat(sep, "\n")

  cat("  [Results]\n")
  row("Sample size (n1, n2)",  fmt(x$n1), fmt(x$n2))
  row("Total sample size (N)", fmt(x$N))
  cat(sep, "\n")
  cat("\n")

  invisible(x)
}
