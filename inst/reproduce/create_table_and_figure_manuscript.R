# Generate manuscript LaTeX tables and PDF/EPS figures from cached
# numerical results.
#
# Reads the .rds files produced by run_table_and_figure_manuscript.R
# (located under "table_and_figure_manuscript/") and writes the
# manuscript tables and figures to the same subdirectory:
#
#   table1_gl_nodes.tex               : Gauss-Legendre rule
#   table2_sample_size_comparison.tex : Sample size comparison
#   table3_kappa_star.tex             : Threshold kappa_star table
#   table4_real_example.tex           : Real trial example (Doody 2014)
#   figure1_anchor_weight.{pdf,eps}   : Anchor weight w* vs kappa
#   figure2_reduction_curve.{pdf,eps} : R_max(kappa) at three power levels
#   published_table_comparison_sozu2015.csv : C_2 of Sozu et al. (2015)
#                                       versus exact and closed-form values
#   published_table_comparison_hungwang2009.csv : n / m1 of Hung and
#                                       Wang (2009) versus r_max(1, beta)
#
# Each table file is a complete table environment (caption, label,
# tabular and table notes), so that the manuscript only needs
# \input{table1_gl_nodes} etc.  The manuscript preamble must load the
# threeparttable package.  Captions refer to equation labels of the
# manuscript (eq:kappa_star) and to the bibliography key Doody2014.
#
# This script does not perform any numerical computation; it only
# formats results and draws figures.
#
# Before running this script, set the working directory to the folder
# that contains the "table_and_figure_manuscript/" subdirectory, that
# is, the folder in which run_table_and_figure_manuscript.R was run.

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------------------------------------
# Output directory for all tables and figures (also where the
# cached .rds files were saved by run_table_and_figure_manuscript.R)
# ------------------------------------------------------------
output_dir <- "table_and_figure_manuscript"
if (!dir.exists(output_dir)) {
  stop("Directory not found: ", output_dir,
       ".  Run run_table_and_figure_manuscript.R first.")
}
out_path <- function(name) file.path(output_dir, name)

# ------------------------------------------------------------
# Load cached results
# ------------------------------------------------------------
required_rds <- c(
  "gl_nodes.rds",
  "anchor_weight.rds",
  "table2_sample_size_comparison.rds",
  "reduction_curve.rds",
  "kappa_star_table.rds",
  "table4_real_example.rds",
  "published_table_comparison.rds"
)
missing_rds <- required_rds[!file.exists(out_path(required_rds))]
if (length(missing_rds) > 0) {
  stop("Missing cached results under ", output_dir, "/: ",
       paste(missing_rds, collapse = ", "),
       ".  Run run_table_and_figure_manuscript.R first.")
}

gl_obj         <- readRDS(out_path("gl_nodes.rds"))
anchor_obj     <- readRDS(out_path("anchor_weight.rds"))
table2_obj     <- readRDS(out_path("table2_sample_size_comparison.rds"))
reduction_obj  <- readRDS(out_path("reduction_curve.rds"))
kappa_star_obj <- readRDS(out_path("kappa_star_table.rds"))
table4_obj     <- readRDS(out_path("table4_real_example.rds"))
published_obj  <- readRDS(out_path("published_table_comparison.rds"))

# ------------------------------------------------------------
# Common figure helper: save .pdf and .eps with cairo devices
# ------------------------------------------------------------
save_fig <- function(plot, name, width, height) {
  ggsave(out_path(paste0(name, ".pdf")),
         plot = plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(out_path(paste0(name, ".eps")),
         plot = plot, width = width, height = height,
         units = "mm", device = cairo_ps,
         fallback_resolution = 600)
}

theme_article <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(size = 9))

# ------------------------------------------------------------
# Common table helper: wrap a tabular in a complete table environment
# with caption, label and optional table notes, and write it to file.
# Each element of notes is a complete \item line.
# ------------------------------------------------------------
write_table_env <- function(tabular, file_out, caption, label,
                            notes = NULL, size = "\\small",
                            tabcolsep = NULL) {
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    size,
    if (!is.null(tabcolsep)) {
      sprintf("\\setlength{\\tabcolsep}{%s}", tabcolsep)
    },
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{threeparttable}",
    tabular,
    if (!is.null(notes)) {
      c("\\begin{tablenotes}", "\\footnotesize", notes,
        "\\end{tablenotes}")
    },
    "\\end{threeparttable}",
    "\\end{table}"
  )
  writeLines(lines, file_out)
  cat(sprintf("Wrote: %s\n", basename(file_out)))
}

# Format a vector of numbers as "$a$, $b$, and $c$" for captions.  Each
# number is set in math mode separately, so that "and" stays in text.
fmt_list <- function(x, fmt) {
  s <- paste0("$", sprintf(fmt, x), "$")
  if (length(s) <= 2L) return(paste(s, collapse = " and "))
  paste0(paste(s[-length(s)], collapse = ", "), ", and ", s[length(s)])
}

# Common color/linetype scheme used in Figure 1 and Figure 2
color_set    <- c("#1B9E77", "#444444", "#E7298A")
linetype_set <- c("dashed",  "solid",   "dotdash")

# ============================================================
# Table 1: Gauss-Legendre nodes and weights
# ============================================================
gl_nodes_n <- length(gl_obj$x)
tab1_lines <- c(
  "\\begin{tabular}{rrr}",
  "\\hline",
  "$i$ & $x_i$ & $\\omega_i$ \\\\",
  "\\hline"
)
for (i in seq_len(gl_nodes_n)) {
  tab1_lines <- c(
    tab1_lines,
    # Math mode so that negative nodes are typeset with a minus sign
    sprintf("%d & $%+.10f$ & $%.10f$ \\\\", i,
            gl_obj$x[i], gl_obj$w[i])
  )
}
tab1_lines <- c(tab1_lines, "\\hline", "\\end{tabular}")
write_table_env(
  tab1_lines, out_path("table1_gl_nodes.tex"),
  caption = sprintf(paste0("Nodes ($x_i$) and weights ($\\omega_i$) of the ",
                           "%d-point Gauss--Legendre quadrature on ",
                           "$[-1, 1]$, computed by the Golub--Welsch ",
                           "algorithm."), gl_nodes_n),
  label = "tab:gl_nodes", size = "\\normalsize"
)

# ============================================================
# Table 2: Sample size comparison
# ============================================================
build_ss_tex <- function(df, power_set, file_out,
                         alpha = table2_obj$alpha,
                         gl_nodes = table2_obj$gl_nodes) {
  # One row per (pattern, rho); for each target power, the columns are
  # N_conv, N_prop and the achieved power of N_prop.  A dagger marks
  # achieved powers below the target.
  n_pow <- length(power_set)
  lines <- c(
    sprintf("\\begin{tabular}{cccc%s}",
            paste(rep("rrr", n_pow), collapse = "")),
    "\\hline",
    paste0(" & & & & ",
           paste(sprintf("\\multicolumn{3}{c}{$1 - \\beta = %.2f$}",
                         power_set), collapse = " & "),
           " \\\\"),
    paste(sprintf("\\cline{%d-%d}", 5 + 3 * (seq_len(n_pow) - 1),
                  7 + 3 * (seq_len(n_pow) - 1)), collapse = " "),
    paste0("$(\\delta_1, \\delta_2)$ & $(\\sigma_1, \\sigma_2)$ & ",
           "$r$ & $\\rho$ & ",
           paste(rep("$N_{\\rm conv}$ & $N_{\\rm prop}$ & power", n_pow),
                 collapse = " & "),
           " \\\\"),
    "\\hline"
  )
  keys <- unique(df[, c("delta1", "delta2", "sd1", "sd2", "r", "rho")])
  prev_key <- ""
  for (i in seq_len(nrow(keys))) {
    k <- keys[i, ]
    cur_key <- sprintf("%s|%s|%s|%s|%s",
                       k$delta1, k$delta2, k$sd1, k$sd2, k$r)
    if (cur_key != prev_key) {
      delta_str <- sprintf("(%.1f, %.1f)", k$delta1, k$delta2)
      sd_str    <- sprintf("(%.1f, %.1f)", k$sd1, k$sd2)
      r_str     <- sprintf("%d", as.integer(k$r))
    } else {
      delta_str <- ""
      sd_str    <- ""
      r_str     <- ""
    }
    prev_key <- cur_key
    cells <- sapply(power_set, function(pw) {
      d <- df[df$delta1 == k$delta1 & df$delta2 == k$delta2 &
                df$sd1 == k$sd1 & df$sd2 == k$sd2 & df$r == k$r &
                df$rho == k$rho & abs(df$power - pw) < 1e-9, ]
      if (nrow(d) != 1L) stop("Table 2: cell not found or not unique.")
      mark <- if (d$pow < pw) "$^{\\dagger}$" else ""
      sprintf("%d & %d & %.4f%s", d$N_conv, d$N_prop, d$pow, mark)
    })
    lines <- c(lines, sprintf("%s & %s & %s & %.1f & %s \\\\",
                              delta_str, sd_str, r_str, k$rho,
                              paste(cells, collapse = " & ")))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")

  caption <- sprintf(paste0(
    "Total sample size ($N$) of the conventional iterative method ",
    "($N_{\\rm conv}$) and the proposed closed-form formula ",
    "($N_{\\rm prop}$), and the achieved co-primary power at ",
    "$N_{\\rm prop}$, for target powers $1 - \\beta$ of %s. The ",
    "significance level is $\\alpha = %s$ (one-sided)."),
    fmt_list(power_set, "%.2f"), format(alpha))
  notes <- sprintf(paste0(
    "\\item[] Achieved power is the bivariate normal power at the ",
    "group sizes $(n_1, n_2)$ of the proposed formula, computed with ",
    "the \\texttt{pbivnorm} package. The proposed formula uses the ",
    "%d-point Gauss--Legendre rule."), gl_nodes)
  short <- df[df$pow < df$power, ]
  if (nrow(short) > 0L) {
    detail <- sprintf(paste0(
      "at $(\\delta_1, \\delta_2, \\sigma_1, \\sigma_2, r, \\rho) = ",
      "(%.1f, %.1f, %.1f, %.1f, %d, %.1f)$ and $1 - \\beta = %.2f$, the ",
      "achieved power at $N_{\\rm prop} = %d$ is %.7f, whereas ",
      "$N_{\\rm conv} = %d$"),
      short$delta1, short$delta2, short$sd1, short$sd2,
      as.integer(short$r), short$rho, short$power,
      as.integer(short$N_prop), short$pow, as.integer(short$N_conv))
    notes <- c(notes, paste0(
      "\\item[$\\dagger$] Achieved power below the target: ",
      paste(detail, collapse = "; "), "."))
  }
  write_table_env(lines, file_out, caption = caption,
                  label = "tab:ss_comparison", notes = notes,
                  size = "\\footnotesize", tabcolsep = "4pt")
}
build_ss_tex(table2_obj$data, table2_obj$power_set,
             out_path("table2_sample_size_comparison.tex"))
write.csv(table2_obj$data, out_path("table2_sample_size_comparison.csv"),
          row.names = FALSE)

# ============================================================
# Table 3: Threshold kappa_star for selected tolerances
# ============================================================
build_kstar_tex <- function(df, file_out,
                            alpha = kappa_star_obj$alpha) {
  power_set   <- sort(unique(df$power))
  epsilon_set <- sort(unique(df$epsilon))
  ncol_eps    <- length(epsilon_set)

  col_spec <- paste0("c", paste(rep("r", ncol_eps), collapse = ""))
  header   <- paste(
    "$1 - \\beta$",
    paste(sprintf("$\\epsilon = %.2f$", epsilon_set),
          collapse = " & "),
    sep = " & "
  )

  lines <- c(
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline"
  )
  for (p in power_set) {
    row_vals <- sapply(epsilon_set, function(e) {
      df$kappa_star[df$power == p & df$epsilon == e]
    })
    lines <- c(lines, sprintf(
      "%.2f & %s \\\\", p,
      paste(sprintf("%.4f", row_vals), collapse = " & ")
    ))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")

  caption <- sprintf(paste0(
    "Threshold effect size ratio $\\kappa^{*}(\\epsilon, \\beta)$ ",
    "in~(\\ref{eq:kappa_star}), above which the maximum reduction rate ",
    "satisfies $R_{\\max}(\\kappa, \\beta) \\leq \\epsilon$, for ",
    "tolerances $\\epsilon$ of %s and target powers $1 - \\beta$ of %s. ",
    "The significance level is $\\alpha = %s$ (one-sided)."),
    fmt_list(epsilon_set, "%.2f"), fmt_list(power_set, "%.2f"),
    format(alpha))
  write_table_env(lines, file_out, caption = caption,
                  label = "tab:kappa_star")
}
build_kstar_tex(kappa_star_obj$data, out_path("table3_kappa_star.tex"))
write.csv(kappa_star_obj$data, out_path("table3_kappa_star.csv"),
          row.names = FALSE)

# ============================================================
# Table 4: Real trial example based on Doody et al. (2014),
# EXPEDITION 1 of solanezumab
# ============================================================
build_real_example_tex <- function(obj, file_out) {
  d <- obj$data
  es1 <- obj$delta1 / obj$sd1
  es2 <- obj$delta2 / obj$sd2

  lines <- c(
    "\\begin{tabular}{lrr}",
    "\\hline",
    " & ADAS-Cog11 & ADCS-ADL \\\\",
    "\\hline",
    sprintf("Mean difference $\\delta$ & %.1f & %.1f \\\\",
            obj$delta1, obj$delta2),
    sprintf("Standard deviation $\\sigma$ & %.0f & %.0f \\\\",
            obj$sd1, obj$sd2),
    sprintf("Standardized effect size $\\delta / \\sigma$ & %.4f & %.4f \\\\",
            es1, es2),
    "\\hline",
    sprintf("Effect-size ratio $\\kappa$ & \\multicolumn{2}{c}{%.4f} \\\\",
            obj$kappa),
    sprintf("Significance level $\\alpha$ (one-sided) & \\multicolumn{2}{c}{%.3f} \\\\",
            obj$alpha),
    sprintf("Target overall power $1 - \\beta$ & \\multicolumn{2}{c}{%.2f} \\\\",
            1 - obj$beta),
    sprintf("Allocation ratio $r$ & \\multicolumn{2}{c}{%d} \\\\", obj$r),
    "\\hline",
    "\\multicolumn{3}{c}{} \\\\[-1.8ex]",
    "$\\rho$ & $N_{\\rm conv}$ & $N_{\\rm prop}$ \\\\",
    "\\hline"
  )
  for (i in seq_len(nrow(d))) {
    lines <- c(lines, sprintf(
      "%.1f & %d & %d \\\\",
      d$rho[i], d$N_conv[i], d$N_prop[i]
    ))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")

  caption <- sprintf(paste0(
    "Design assumptions and sample sizes for EXPEDITION~1, the Phase~III ",
    "trial of solanezumab for mild-to-moderate Alzheimer's disease",
    "~\\citep{Doody2014}. The upper block lists the design ",
    "parameters from the trial protocol, with overall power ",
    "$1 - \\beta = %.2f$. The lower block reports the total sample size ",
    "$N$ of the conventional iterative method ($N_{\\rm conv}$) and the ",
    "proposed closed-form formula ($N_{\\rm prop}$) for %d values of ",
    "$\\rho$. At $\\rho = 0$, the proposed formula gives $N = %d$."),
    1 - obj$beta, nrow(d), as.integer(d$N_prop[d$rho == 0]))
  write_table_env(lines, file_out, caption = caption,
                  label = "tab:real_example")
}
build_real_example_tex(table4_obj, out_path("table4_real_example.tex"))
write.csv(table4_obj$data, out_path("table4_real_example.csv"),
          row.names = FALSE)

# ============================================================
# Figure 1: Anchor weight w* as a function of kappa
# (curve computed in run_table_and_figure_manuscript.R)
# ============================================================
power_set <- anchor_obj$power_set

df_fig1 <- anchor_obj$data
df_fig1$target <- factor(sprintf("%.2f", df_fig1$power),
                         levels = sprintf("%.2f", power_set))

target_labels <- lapply(power_set, function(p) {
  bquote(1 - beta == .(sprintf("%.2f", p)))
})

p1 <- ggplot(df_fig1, aes(x = kappa, y = w_star,
                          linetype = target, colour = target)) +
  geom_hline(yintercept = c(0.5, 1.0),
             linetype = "dotted", colour = "#888888",
             linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = linetype_set, labels = target_labels) +
  scale_colour_manual(values = color_set, labels = target_labels) +
  scale_x_continuous(breaks = seq(1, 4, by = 0.5)) +
  scale_y_continuous(limits = c(0.30, 1.05),
                     breaks = seq(0.3, 1.0, by = 0.1)) +
  labs(x = expression(kappa),
       y = expression(italic(w)^"*")) +
  theme_article +
  theme(legend.text      = element_text(size = 9),
        legend.key.width = unit(1.0, "cm"))
save_fig(p1, "figure1_anchor_weight", width = 120, height = 90)
cat("Wrote: figure1_anchor_weight.{pdf,eps}\n")

# ============================================================
# Figure 2: R_max(kappa) curve at three power levels
# ============================================================
df_fig2 <- reduction_obj$data
df_fig2$target <- factor(sprintf("%.2f", df_fig2$power),
                         levels = sprintf("%.2f", power_set))
df_fig2$R_pct <- df_fig2$R_max * 100

p2 <- ggplot(df_fig2, aes(x = kappa, y = R_pct,
                          linetype = target, colour = target)) +
  geom_hline(yintercept = c(5, 10, 15),
             linetype = "dotted", colour = "#888888",
             linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = linetype_set, labels = target_labels) +
  scale_colour_manual(values = color_set, labels = target_labels) +
  scale_x_continuous(breaks = seq(1.0, 2.0, by = 0.2)) +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  labs(x = expression(kappa),
       y = expression(R[max] * " (%)")) +
  theme_article +
  theme(legend.text      = element_text(size = 9),
        legend.key.width = unit(1.0, "cm"))
save_fig(p2, "figure2_reduction_curve", width = 120, height = 90)
cat("Wrote: figure2_reduction_curve.{pdf,eps}\n")

# ============================================================
# Comparison with published values (supplementary output, not a
# manuscript table)
# ============================================================
write.csv(published_obj$sozu,
          out_path("published_table_comparison_sozu2015.csv"),
          row.names = FALSE)
cat("Wrote: published_table_comparison_sozu2015.csv\n")
write.csv(published_obj$hung_wang,
          out_path("published_table_comparison_hungwang2009.csv"),
          row.names = FALSE)
cat("Wrote: published_table_comparison_hungwang2009.csv\n")

# ============================================================
# Reporting cached metadata
# ============================================================
cat("\n----- Cached computation metadata -----\n")
cat(sprintf("Table 1 elapsed: %.3f s (%s)\n",
            gl_obj$elapsed,
            format(gl_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Anchor weight curve elapsed: %.3f s (%s)\n",
            anchor_obj$elapsed,
            format(anchor_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Table 2 elapsed: %.2f s (%s)\n",
            table2_obj$elapsed,
            format(table2_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Reduction curve elapsed: %.3f s (%s)\n",
            reduction_obj$elapsed,
            format(reduction_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("kappa_star table elapsed: %.3f s (%s)\n",
            kappa_star_obj$elapsed,
            format(kappa_star_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Table 4 elapsed: %.3f s (%s)\n",
            table4_obj$elapsed,
            format(table4_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Published table comparison elapsed: %.3f s (%s)\n",
            published_obj$elapsed,
            format(published_obj$computed, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("\nAll outputs written to: %s/\n", output_dir))
