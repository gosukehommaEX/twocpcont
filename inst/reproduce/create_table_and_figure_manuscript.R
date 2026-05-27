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
#
# This script does not perform any numerical computation; it only
# formats results and draws figures.
#
# Before running this script, set the working directory to the folder
# that contains the "table_and_figure_manuscript/" subdirectory.

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
  "table2_sample_size_comparison.rds",
  "reduction_curve.rds",
  "kappa_star_table.rds",
  "table4_real_example.rds"
)
missing_rds <- required_rds[!file.exists(out_path(required_rds))]
if (length(missing_rds) > 0) {
  stop("Missing cached results under ", output_dir, "/: ",
       paste(missing_rds, collapse = ", "),
       ".  Run run_table_and_figure_manuscript.R first.")
}

gl_obj         <- readRDS(out_path("gl_nodes.rds"))
table2_obj     <- readRDS(out_path("table2_sample_size_comparison.rds"))
reduction_obj  <- readRDS(out_path("reduction_curve.rds"))
kappa_star_obj <- readRDS(out_path("kappa_star_table.rds"))
table4_obj     <- readRDS(out_path("table4_real_example.rds"))

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
  "$i$ & $x_i$ & $w_i$ \\\\",
  "\\hline"
)
for (i in seq_len(gl_nodes_n)) {
  tab1_lines <- c(
    tab1_lines,
    sprintf("%d & %+.10f & %.10f \\\\", i,
            gl_obj$x[i], gl_obj$w[i])
  )
}
tab1_lines <- c(tab1_lines, "\\hline", "\\end{tabular}")
writeLines(tab1_lines, out_path("table1_gl_nodes.tex"))
cat("Wrote: table1_gl_nodes.tex\n")

# ============================================================
# Table 2: Sample size comparison
# ============================================================
build_ss_tex <- function(df, file_out) {
  lines <- c(
    "\\begin{tabular}{cccccrrrr}",
    "\\hline",
    paste0("$(\\delta_1, \\delta_2)$ & $(\\sigma_1, \\sigma_2)$ & ",
           "$r$ & $\\rho$ & ",
           "$N_{\\rm conv}$ & $N_{\\rm prop}$ & diff & power \\\\"),
    "\\hline"
  )
  prev_key <- ""
  for (i in seq_len(nrow(df))) {
    cur_key <- sprintf("%s|%s|%s|%s|%s",
                       df$delta1[i], df$delta2[i],
                       df$sd1[i], df$sd2[i], df$r[i])
    if (cur_key != prev_key) {
      delta_str <- sprintf("(%.1f, %.1f)", df$delta1[i], df$delta2[i])
      sd_str    <- sprintf("(%.1f, %.1f)", df$sd1[i],    df$sd2[i])
      r_str     <- sprintf("%d", as.integer(df$r[i]))
    } else {
      delta_str <- ""
      sd_str    <- ""
      r_str     <- ""
    }
    prev_key <- cur_key
    diff_str <- if (df$diff[i] == 0L) "0" else sprintf("%+d", df$diff[i])
    lines <- c(lines, sprintf(
      "%s & %s & %s & %.1f & %d & %d & %s & %.4f \\\\",
      delta_str, sd_str, r_str, df$rho[i],
      df$N_conv[i], df$N_prop[i], diff_str, df$pow[i]
    ))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")
  writeLines(lines, file_out)
  cat(sprintf("Wrote: %s\n", basename(file_out)))
}
build_ss_tex(table2_obj$data, out_path("table2_sample_size_comparison.tex"))
write.csv(table2_obj$data, out_path("table2_sample_size_comparison.csv"),
          row.names = FALSE)

# ============================================================
# Table 3: Threshold kappa_star for selected tolerances
# ============================================================
build_kstar_tex <- function(df, file_out) {
  power_set   <- sort(unique(df$power))
  epsilon_set <- sort(unique(df$epsilon))
  ncol_eps    <- length(epsilon_set)

  col_spec <- paste0("c", paste(rep("r", ncol_eps), collapse = ""))
  header   <- paste(
    "$1 - \\beta$",
    paste(sprintf("$\\epsilon = %.3f$", epsilon_set),
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
  writeLines(lines, file_out)
  cat(sprintf("Wrote: %s\n", basename(file_out)))
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
  writeLines(lines, file_out)
  cat(sprintf("Wrote: %s\n", basename(file_out)))
}
build_real_example_tex(table4_obj, out_path("table4_real_example.tex"))
write.csv(table4_obj$data, out_path("table4_real_example.csv"),
          row.names = FALSE)

# ============================================================
# Figure 1: Anchor weight w* as a function of kappa
# Implements eq. (eq:w_star) directly:
#   w_star = 1/2 + (kappa - 1) (u_p + z_alpha) phi(u_p) /
#                  [(1 + kappa) |2 log p| p]
# ============================================================
w_star_fun <- function(kappa, beta = 0.2, alpha = 0.025) {
  z_alpha <- qnorm(1 - alpha)
  target  <- 1 - beta
  p       <- sqrt(target)
  u_p     <- qnorm(p)
  phi_u_p <- dnorm(u_p)
  L_abs   <- abs(2 * log(p))
  w <- 0.5 +
    (kappa - 1) * (u_p + z_alpha) * phi_u_p /
    ((1 + kappa) * L_abs * p)
  pmin(1, pmax(0, w))
}

kappa_seq <- seq(1, 4, length.out = 301)
power_set <- c(0.70, 0.80, 0.90)
beta_set  <- 1 - power_set

df_fig1 <- do.call(rbind, lapply(seq_along(beta_set), function(k) {
  data.frame(
    kappa  = kappa_seq,
    w_star = w_star_fun(kappa_seq, beta = beta_set[k]),
    target = factor(sprintf("%.2f", power_set[k]),
                    levels = sprintf("%.2f", power_set))
  )
}))

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
  geom_hline(yintercept = c(1, 5),
             linetype = "dotted", colour = "#888888",
             linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = linetype_set, labels = target_labels) +
  scale_colour_manual(values = color_set, labels = target_labels) +
  scale_x_continuous(breaks = seq(1.0, 2.0, by = 0.2)) +
  scale_y_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 5)) +
  labs(x = expression(kappa),
       y = expression(R[max] * " (%)")) +
  theme_article +
  theme(legend.text      = element_text(size = 9),
        legend.key.width = unit(1.0, "cm"))
save_fig(p2, "figure2_reduction_curve", width = 120, height = 90)
cat("Wrote: figure2_reduction_curve.{pdf,eps}\n")

# ============================================================
# Reporting cached metadata
# ============================================================
cat("\n----- Cached computation metadata -----\n")
cat(sprintf("Table 1 elapsed: %.3f s (%s)\n",
            gl_obj$elapsed,
            format(gl_obj$computed, "%Y-%m-%d %H:%M:%S")))
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
cat(sprintf("\nAll outputs written to: %s/\n", output_dir))
