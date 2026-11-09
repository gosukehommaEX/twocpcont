# =============================================================================
# table_and_figure_manuscript.R
#
# Rendering stage for the manuscript figures and tables of the twocpcont
# article. This script reads the data written by data_generate_main.R
# (inst/scripts/data/*.rds) and produces the figures, the LaTeX tables and
# the ledger of the numbers cited in the text. It performs no numerical
# computation, so figures and tables can be adjusted without re-running
# data_generate_main.R.
#
# Usage
#   Open twocpcont.Rproj (the working directory must be the package root)
#   and source this file:
#     source("inst/scripts/table_and_figure_manuscript.R")
#
# Prerequisite
#   inst/scripts/data/*.rds produced by data_generate_main.R.
#
# Required packages: ggplot2
#
# Outputs (relative to the package root)
#   inst/scripts/results/figure1_anchor_weight.eps      Figure 1 (800 dpi EPS)
#   inst/scripts/results/figure2_reduction_curve.eps    Figure 2 (800 dpi EPS)
#   inst/scripts/results/figure3_reduction_rho.eps      Figure 3 (800 dpi EPS)
#   inst/scripts/results/table1_gl_nodes.tex            Table 1
#   inst/scripts/results/table2_sample_size_comparison.tex  Table 2
#   inst/scripts/results/table3_kappa_star.tex          Table 3
#   inst/scripts/results/table4_real_example.tex        Table 4
#   inst/scripts/results/manuscript_numbers.csv         Ledger of cited numbers
#
# In the manuscript, each table is included with \input{table1_gl_nodes} and
# each figure with \includegraphics[width=\textwidth]{figure1_anchor_weight.eps};
# the captions of the figures are written in the manuscript, those of the
# tables here. Each table file is a complete table environment (caption,
# label, tabular and table notes), so the manuscript preamble must load the
# threeparttable package. Captions refer to equation labels of the
# manuscript (eq:kappa_star) and to the bibliography key Doody2014.
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "twocpcont")) {
  stop("Run this script from the root of the twocpcont package ",
       "(open twocpcont.Rproj first).", call. = FALSE)
}

library(ggplot2)

# --------------------------------------------------------------------------- #
# Read the data produced by data_generate_main.R
# --------------------------------------------------------------------------- #
data_dir    <- file.path("inst", "scripts", "data")
results_dir <- file.path("inst", "scripts", "results")
slides_dir  <- file.path(results_dir, "slides")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
out_path <- function(name) file.path(results_dir, name)

read_data <- function(name) {
  path <- file.path(data_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run data_generate_main.R first.", path),
         call. = FALSE)
  }
  readRDS(path)
}

gl_obj         <- read_data("tab1_gl_nodes")
anchor_obj     <- read_data("fig1_anchor_weight")
table2_obj     <- read_data("tab2_sample_size")
phi2_obj       <- read_data("phi2_evaluations")
reduction_obj  <- read_data("fig2_reduction_curve")
kappa_star_obj <- read_data("tab3_kappa_star")
table4_obj     <- read_data("tab4_real_example")
published_obj  <- read_data("published_tables")
reduction_rho_obj <- read_data("fig3_reduction_rho")

# --------------------------------------------------------------------------- #
# Plotting constants
# --------------------------------------------------------------------------- #
BASE_SIZE <- 10                  # base font size of the article figures
FIG_WIDTH  <- 120 / 25.4         # inches (120 mm)
FIG_HEIGHT <- 90 / 25.4          # inches (90 mm)
DPI       <- 800                 # resolution of the EPS output

# Figures for the slides of a conference talk (4:3 beamer frame, 128 mm by
# 96 mm, text width 116 mm): saved at the size at which they are shown, so
# that BASE_SIZE_SLIDE is the size of the letters on the slide
SLIDE_WIDTH     <- 4.6           # inches (116 mm, \textwidth of the frame)
SLIDE_HEIGHT    <- 2.9           # inches (below the frame title and subtitle)
BASE_SIZE_SLIDE <- 10            # points

# --------------------------------------------------------------------------- #
# Rendering helpers
# --------------------------------------------------------------------------- #
theme_paper <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(size = 9))
}

# Save a figure as an 800 dpi EPS file, the only format used in the
# manuscript. The cairo device embeds the fonts; fallback_resolution is the
# resolution of the parts rasterised for transparency (72 dpi by default).
# A PDF for viewing is written next to it only when preview = TRUE, and is
# never copied to the manuscript.
use_cairo <- isTRUE(capabilities("cairo"))
save_figure <- function(plot, name, width, height, dir = results_dir,
                        preview = FALSE) {
  f_eps <- file.path(dir, paste0(name, ".eps"))
  if (use_cairo) {
    ggsave(f_eps, plot, device = grDevices::cairo_ps, width = width,
           height = height, units = "in", dpi = DPI,
           fallback_resolution = DPI, onefile = FALSE)
  } else {
    ggsave(f_eps, plot, device = "eps", width = width, height = height,
           units = "in", dpi = DPI)
  }
  if (preview) {
    ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = width,
           height = height, units = "in")
  }
  message("Figure saved: ", f_eps)
}

# Save a figure for the slides at the size of the frame
save_slide_figure <- function(plot, name) {
  if (!dir.exists(slides_dir)) dir.create(slides_dir, recursive = TRUE)
  save_figure(plot + theme_paper(BASE_SIZE_SLIDE), name, width = SLIDE_WIDTH,
              height = SLIDE_HEIGHT, dir = slides_dir)
}

# The .tex files are written with LF line endings on every operating
# system, so that they are identical (MD5) wherever they are produced
write_tex <- function(lines, filename) {
  con <- file(filename, open = "wb")
  on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
  message("LaTeX table saved: ", filename)
}

# Fixed-point formatting; NA becomes "---"
fmt <- function(x, digits = 3, na = "---") {
  ifelse(is.na(x), na, sprintf(paste0("%.", digits, "f"), x))
}

# Wrap a tabular in a complete table environment with caption, label and
# optional table notes, and write it to file. Each element of notes is a
# complete \item line.
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
  write_tex(lines, file_out)
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

# ============================================================
# Figure 1: Anchor weight w* as a function of kappa
# (curve computed in data_generate_main.R)
# ============================================================
power_set <- anchor_obj$power_set

df_fig1 <- anchor_obj$data
df_fig1$target <- factor(sprintf("%.2f", df_fig1$power),
                         levels = sprintf("%.2f", power_set))

target_labels <- lapply(power_set, function(p) {
  bquote(1 - beta == .(sprintf("%.2f", p)))
})

# Closed-form w* (solid) and the exact root of (S2) (dashed) for each
# target power. The axis and the legend name no symbol, so that the figure
# does not depend on the notation of the weight; the caption gives it.
df_fig1_long <- rbind(
  data.frame(kappa = df_fig1$kappa, target = df_fig1$target,
             w = df_fig1$w_star, curve = "closed"),
  data.frame(kappa = df_fig1$kappa, target = df_fig1$target,
             w = df_fig1$w_root, curve = "exact")
)
df_fig1_long$curve <- factor(df_fig1_long$curve, levels = c("closed", "exact"))

p1 <- ggplot(df_fig1_long, aes(x = kappa, y = w, colour = target,
                               linetype = curve,
                               group = interaction(target, curve))) +
  geom_hline(yintercept = c(0.5, 1.0),
             linetype = "dotted", colour = "#888888",
             linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = color_set, labels = target_labels) +
  scale_linetype_manual(values = c(closed = "solid", exact = "dashed"),
                        labels = c(closed = "Closed form",
                                   exact = "Exact root")) +
  scale_x_continuous(breaks = seq(1, 4, by = 0.5)) +
  scale_y_continuous(limits = c(0.30, 1.05),
                     breaks = seq(0.3, 1.0, by = 0.1)) +
  labs(x = expression(kappa), y = "Anchor weight") +
  theme_paper() +
  theme(legend.text      = element_text(size = 9),
        legend.key.width = unit(1.0, "cm"),
        legend.box       = "vertical",
        legend.spacing.y = unit(2, "pt"),
        legend.margin    = margin(0, 0, 0, 0))
save_figure(p1, "figure1_anchor_weight", width = FIG_WIDTH,
            height = 100 / 25.4)

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
  theme_paper() +
  theme(legend.text      = element_text(size = 9),
        legend.key.width = unit(1.0, "cm"))
save_figure(p2, "figure2_reduction_curve", width = FIG_WIDTH,
            height = FIG_HEIGHT)

# ============================================================
# Figure 3: reduction of the sample size on the continuous scale as a
# function of rho, with its limit R_max as rho -> 1 (points at rho = 1)
# (computed in data_generate_main.R)
# ============================================================
df_fig3  <- reduction_rho_obj$data
lim_fig3 <- reduction_rho_obj$limit
kappa_levels <- sprintf("%.1f", sort(unique(df_fig3$kappa)))
df_fig3$kappa_f  <- factor(sprintf("%.1f", df_fig3$kappa),
                           levels = kappa_levels)
lim_fig3$kappa_f <- factor(sprintf("%.1f", lim_fig3$kappa),
                           levels = kappa_levels)
df_fig3$panel  <- sprintf("%.2f", df_fig3$power)
lim_fig3$panel <- sprintf("%.2f", lim_fig3$power)
# Okabe-Ito colours (colour-blind friendly)
kappa_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2",
                "#D55E00")

p3 <- ggplot(df_fig3, aes(x = rho, y = 100 * reduction, colour = kappa_f)) +
  geom_line(linewidth = 0.6) +
  geom_point(data = lim_fig3, aes(x = 1, y = 100 * R_max, colour = kappa_f),
             size = 1.5) +
  facet_grid(. ~ panel, labeller = label_bquote(cols = 1 - beta == .(panel))) +
  scale_colour_manual(values = kappa_cols, name = expression(kappa)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) +
  labs(x = expression(rho), y = "Reduction in sample size (%)") +
  theme_paper() +
  theme(legend.title = element_text(size = 10),
        legend.text  = element_text(size = 9))
save_figure(p3, "figure3_reduction_rho", width = 160 / 25.4,
            height = FIG_HEIGHT)

# =============================================================================
# Ledger of the numbers cited in the text: results/manuscript_numbers.csv
#   key          identifier used in number_locations.csv
#   value        unrounded value
#   digits       number of decimals in the text
#   text         value formatted with 'digits' (to be transcribed)
#   description  what the number is
#   source       data file and element
# Percentages are entered as fractions; number_locations.csv gives scale
# 100 where the text writes them in percent.
# =============================================================================
ledger <- list()
add_num <- function(key, value, digits, description, source) {
  if (key %in% names(ledger)) {
    stop(sprintf("Duplicated ledger key '%s'.", key), call. = FALSE)
  }
  if (length(value) != 1L || !is.finite(value)) {
    stop(sprintf("Ledger value '%s' must be a single finite number.", key),
         call. = FALSE)
  }
  ledger[[key]] <<- data.frame(key = key, value = value, digits = digits,
                               text = sprintf(paste0("%.", digits, "f"),
                                              value),
                               description = description, source = source,
                               stringsAsFactors = FALSE)
}

# A difference that the text states between two values shown in a table is
# computed from the values as the table displays them (formatted with the
# decimals of the table by sprintf(), as fmt() does), so that the reader
# obtains it from the table and the text needs no remark on rounding.
#   a, b          unrounded values shown in the table
#   digits_table  number of decimals in the table
add_diff <- function(key, a, b, digits_table, digits, description, source) {
  shown <- function(x) as.numeric(sprintf(paste0("%.", digits_table, "f"), x))
  add_num(key, shown(a) - shown(b), digits, description, source)
}

# --------------------------------------------------------------------------- #
# Introduction
# --------------------------------------------------------------------------- #
add_num("intro_power_indep", 0.8 ^ 2, 2,
        "Overall power of two independent endpoints with power 0.80 each",
        "0.8^2")

# --------------------------------------------------------------------------- #
# Section 2: anchor weight, maximum reduction rate, threshold
# --------------------------------------------------------------------------- #
w_one <- anchor_obj$w_one
for (pw in c(0.80, 0.85, 0.90)) {
  add_num(sprintf("w_one_kappa_%02.0f", 100 * pw),
          w_one$kappa[abs(w_one$power - pw) < 1e-9], 2,
          sprintf("kappa at which w* reaches one, 1 - beta = %.2f", pw),
          "fig1_anchor_weight.rds: w_one")
}
add_num("w_gap_max", max(anchor_obj$gap$gap), 2,
        "Largest excess of w* over the exact root of (S2), 1 < kappa <= 3",
        "fig1_anchor_weight.rds: gap")

rmax_curve <- reduction_obj$data
rmax_k1 <- function(pw) {
  rmax_curve$R_max[rmax_curve$kappa == 1 & abs(rmax_curve$power - pw) < 1e-9]
}
for (pw in c(0.80, 0.85, 0.90)) {
  add_num(sprintf("rmax_k1_%02.0f", 100 * pw), rmax_k1(pw), 3,
          sprintf("R_max(1, beta), 1 - beta = %.2f", pw),
          "fig2_reduction_curve.rds: data")
}
add_num("rmax_k1_80_r0", rmax_k1(0.80), 2,
        "R_max(1, 0.20), rounded to a whole percent",
        "fig2_reduction_curve.rds: data")

kst <- kappa_star_obj$data
kstar_of <- function(eps, pw) {
  kst$kappa_star[abs(kst$epsilon - eps) < 1e-9 & abs(kst$power - pw) < 1e-9]
}
add_num("kstar_05_80", kstar_of(0.05, 0.80), 2,
        "kappa_star(0.05, 0.20)", "tab3_kappa_star.rds: data")
add_num("kstar_15_80", kstar_of(0.15, 0.80), 2,
        "kappa_star(0.15, 0.20)", "tab3_kappa_star.rds: data")
add_num("kstar_05_90", kstar_of(0.05, 0.90), 2,
        "kappa_star(0.05, 0.10)", "tab3_kappa_star.rds: data")
add_num("kstar_15_90", kstar_of(0.15, 0.90), 2,
        "kappa_star(0.15, 0.10)", "tab3_kappa_star.rds: data")
add_num("kstar_05_max",
        max(kst$kappa_star[abs(kst$epsilon - 0.05) < 1e-9]), 1,
        "Largest kappa_star(0.05, beta) over the three target powers",
        "tab3_kappa_star.rds: data")

# --------------------------------------------------------------------------- #
# Section 3.1 and 3.2: sample size comparison (Table 2)
# --------------------------------------------------------------------------- #
d2 <- table2_obj$data
kappa_pattern <- function(dd1, dd2, s1, s2) {
  unique(d2$kappa[d2$delta1 == dd1 & d2$delta2 == dd2 & d2$sd1 == s1 &
                    d2$sd2 == s2])
}
add_num("kappa_p1", kappa_pattern(0.3, 0.3, 1, 1), 1,
        "kappa of (delta1, delta2, sigma1, sigma2) = (0.3, 0.3, 1, 1)",
        "tab2_sample_size.rds: data")
add_num("kappa_p2", kappa_pattern(0.7, 0.5, 1, 1), 1,
        "kappa of (0.7, 0.5, 1, 1)", "tab2_sample_size.rds: data")
add_num("kappa_p3", kappa_pattern(0.3, 0.3, 1, 1.5), 1,
        "kappa of (0.3, 0.3, 1, 1.5)", "tab2_sample_size.rds: data")
add_num("kappa_p4", kappa_pattern(0.7, 0.5, 1, 1.5), 1,
        "kappa of (0.7, 0.5, 1, 1.5)", "tab2_sample_size.rds: data")

add_num("n_cells", nrow(d2), 0, "Number of design combinations",
        "tab2_sample_size.rds: data")
add_num("n_cells_per_power", sum(abs(d2$power - 0.80) < 1e-9), 0,
        "Number of combinations at each target power",
        "tab2_sample_size.rds: data")
same <- d2$N_prop == d2$N_conv
n_match <- function(pw) sum(same[abs(d2$power - pw) < 1e-9])
# The text states one number for the two target powers 0.80 and 0.85
# ("all 32 combinations at 0.80 and 0.85"), so the smaller count is cited
add_num("n_match_80_85", min(n_match(0.80), n_match(0.85)), 0,
        "Combinations with N_prop = N_conv, smaller of 1 - beta = 0.80, 0.85",
        "tab2_sample_size.rds: data")
add_num("n_match_90", n_match(0.90), 0,
        "Combinations with N_prop = N_conv, 1 - beta = 0.90",
        "tab2_sample_size.rds: data")
add_num("n_match_total", sum(same), 0,
        "Combinations with N_prop = N_conv, all target powers",
        "tab2_sample_size.rds: data")

mis <- d2[!same, , drop = FALSE][1L, ]
add_num("mismatch_N_prop", mis$N_prop, 0,
        "N_prop in the combination with N_prop != N_conv",
        "tab2_sample_size.rds: data")
add_num("mismatch_N_conv", mis$N_conv, 0,
        "N_conv in the combination with N_prop != N_conv",
        "tab2_sample_size.rds: data")
add_num("mismatch_n2_exact", mis$n2_exact, 4,
        "Exact continuous n2 in that combination",
        "tab2_sample_size.rds: data$n2_exact")
add_num("mismatch_n2_cf", mis$n2_cont, 4,
        "Closed-form continuous n2 in that combination",
        "tab2_sample_size.rds: data$n2_cont")
add_num("mismatch_n2_diff", mis$n2_exact - mis$n2_cont, 3,
        "Exact minus closed-form continuous n2 in that combination",
        "tab2_sample_size.rds: data")
add_num("mismatch_power", mis$pow, 7,
        "Achieved power at N_prop in that combination",
        "tab2_sample_size.rds: data$pow")

for (pw in c(0.80, 0.85, 0.90)) {
  pw_ok <- d2$pow[same & abs(d2$power - pw) < 1e-9]
  add_num(sprintf("pow_min_%02.0f", 100 * pw), min(pw_ok), 4,
          sprintf("Smallest achieved power, N_prop = N_conv, 1 - beta = %.2f",
                  pw), "tab2_sample_size.rds: data$pow")
  add_num(sprintf("pow_max_%02.0f", 100 * pw), max(pw_ok), 4,
          sprintf("Largest achieved power, N_prop = N_conv, 1 - beta = %.2f",
                  pw), "tab2_sample_size.rds: data$pow")
}

# N_prop over rho for one design pattern (rows of Table 2)
N_over_rho <- function(dd1, dd2, s1, s2, rr, pw) {
  d <- d2[d2$delta1 == dd1 & d2$delta2 == dd2 & d2$sd1 == s1 &
            d2$sd2 == s2 & d2$r == rr & abs(d2$power - pw) < 1e-9, ]
  d$N_prop[order(d$rho)]
}
red_frac <- function(N_vec) (N_vec[1] - N_vec[length(N_vec)]) / N_vec[1]

N_k1_80 <- N_over_rho(0.3, 0.3, 1, 1, 1, 0.80)
N_k1_90 <- N_over_rho(0.3, 0.3, 1, 1, 1, 0.90)
add_num("N_k1_80_rho0", N_k1_80[1], 0, "N_prop, kappa = 1, r = 1, 0.80, rho = 0",
        "tab2_sample_size.rds: data")
add_num("N_k1_80_rho8", N_k1_80[4], 0,
        "N_prop, kappa = 1, r = 1, 0.80, rho = 0.8",
        "tab2_sample_size.rds: data")
add_num("red_k1_80", red_frac(N_k1_80), 2,
        "Reduction of N_prop from rho = 0 to 0.8, kappa = 1, r = 1, 0.80",
        "tab2_sample_size.rds: data")
add_num("N_k1_90_rho0", N_k1_90[1], 0, "N_prop, kappa = 1, r = 1, 0.90, rho = 0",
        "tab2_sample_size.rds: data")
add_num("N_k1_90_rho8", N_k1_90[4], 0,
        "N_prop, kappa = 1, r = 1, 0.90, rho = 0.8",
        "tab2_sample_size.rds: data")
add_num("red_k1_90", red_frac(N_k1_90), 2,
        "Reduction of N_prop from rho = 0 to 0.8, kappa = 1, r = 1, 0.90",
        "tab2_sample_size.rds: data")

red_mid <- unlist(lapply(c(0.80, 0.85, 0.90), function(pw) {
  c(red_frac(N_over_rho(0.7, 0.5, 1, 1, 1, pw)),
    red_frac(N_over_rho(0.7, 0.5, 1, 1, 2, pw)),
    red_frac(N_over_rho(0.3, 0.3, 1, 1.5, 1, pw)),
    red_frac(N_over_rho(0.3, 0.3, 1, 1.5, 2, pw)))
}))
add_num("red_mid_max", max(red_mid), 2,
        "Largest reduction from rho = 0 to 0.8 for kappa = 1.4 and 1.5",
        "tab2_sample_size.rds: data")
add_num("n_r2_match", sum(same & d2$r == 2), 0,
        "Combinations with r = 2 and N_prop = N_conv",
        "tab2_sample_size.rds: data")

# --------------------------------------------------------------------------- #
# Section 3.3: published tables and computational cost
# --------------------------------------------------------------------------- #
pub <- published_obj$sozu
lo  <- pub$rho <= 0.8
add_num("pub_n", nrow(pub), 0, "Tabulated values of Sozu et al. (2015)",
        "published_tables.rds: sozu")
add_num("pub_exact_match", sum(abs(pub$dev_exact) < 1e-9), 0,
        "Tabulated values reproduced by the exact root to three decimals",
        "published_tables.rds: sozu")
add_num("pub_lo_n", sum(lo), 0, "Tabulated values with rho <= 0.8",
        "published_tables.rds: sozu")
add_num("pub_lo_match", sum(abs(pub$dev_cf[lo]) < 1e-9), 0,
        "Values with rho <= 0.8 reproduced by the closed form",
        "published_tables.rds: sozu")
add_num("pub_hi_n", sum(!lo), 0, "Tabulated values with rho = 0.95",
        "published_tables.rds: sozu")
add_num("pub_hi_match", sum(abs(pub$dev_cf[!lo]) < 1e-9), 0,
        "Values with rho = 0.95 reproduced by the closed form",
        "published_tables.rds: sozu")
add_num("pub_hi_maxdev", max(abs(pub$dev_cf[!lo])), 3,
        "Largest deviation of the closed form at rho = 0.95",
        "published_tables.rds: sozu")

ev <- phi2_obj$data
add_num("phi2_median", stats::median(ev$n_phi2_conv), 1,
        "Median number of Phi2 evaluations of the conventional search",
        "phi2_evaluations.rds: data")
add_num("phi2_min", min(ev$n_phi2_conv), 0,
        "Smallest number of Phi2 evaluations of the conventional search",
        "phi2_evaluations.rds: data")
add_num("phi2_max", max(ev$n_phi2_conv), 0,
        "Largest number of Phi2 evaluations of the conventional search",
        "phi2_evaluations.rds: data")
add_num("phi2_argmax_kappa", ev$kappa[which.max(ev$n_phi2_conv)], 1,
        "kappa of the combination with the largest number of evaluations",
        "phi2_evaluations.rds: data")

# --------------------------------------------------------------------------- #
# Section 3.4: maximum reduction rate and threshold
# --------------------------------------------------------------------------- #
add_num("cs_increase_90", 1 / (1 - rmax_k1(0.90)) - 1, 3,
        "1 / {1 - R_max(1, 0.10)} - 1 (Chuang-Stein and Li)",
        "fig2_reduction_curve.rds: data")
hw <- published_obj$hung_wang
for (pw in c(0.80, 0.85, 0.90)) {
  add_num(sprintf("hw_ratio_%02.0f", 100 * pw),
          hw$ratio_rmax[abs(hw$power - pw) < 1e-9], 2,
          sprintf("1 / {1 - R_max(1, beta)} (Hung and Wang), 1 - beta = %.2f",
                  pw), "published_tables.rds: hung_wang")
}
add_num("kappa_varga", reduction_obj$varga$kappa, 2,
        "kappa of the example of Varga et al. (effect sizes 0.4 and 0.5)",
        "fig2_reduction_curve.rds: varga")
add_num("rmax_varga", reduction_obj$varga$R_max, 3,
        "R_max(1.25, 0.20) for the example of Varga et al.",
        "fig2_reduction_curve.rds: varga")

# --------------------------------------------------------------------------- #
# Section 4: EXPEDITION 1 (Table 4)
# --------------------------------------------------------------------------- #
t4 <- table4_obj$data
add_num("trial_es1", table4_obj$delta1 / table4_obj$sd1, 3,
        "Standardized effect size of ADAS-Cog11", "tab4_real_example.rds")
add_num("trial_es2", table4_obj$delta2 / table4_obj$sd2, 3,
        "Standardized effect size of ADCS-ADL", "tab4_real_example.rds")
add_num("trial_kappa", table4_obj$kappa, 3, "Effect size ratio of EXPEDITION 1",
        "tab4_real_example.rds: kappa")
add_num("trial_power_product", 0.89 * 0.98, 2,
        "Product of the marginal powers stated in the protocol",
        "0.89 * 0.98 (protocol)")
add_num("trial_N_rho0", t4$N_prop[t4$rho == 0], 0, "N_prop at rho = 0",
        "tab4_real_example.rds: data")
add_num("trial_N_rho8", t4$N_prop[t4$rho == 0.8], 0, "N_prop at rho = 0.8",
        "tab4_real_example.rds: data")
add_diff("trial_N_red", t4$N_prop[t4$rho == 0], t4$N_prop[t4$rho == 0.8], 0,
         0, "Reduction of N_prop from rho = 0 to 0.8 (Table 4)",
         "tab4_real_example.rds: data")
trial_red <- t4$reduction_pct[t4$rho == 0.8] / 100
add_num("trial_red", trial_red, 2,
        "Relative reduction of N_prop from rho = 0 to 0.8",
        "tab4_real_example.rds: data$reduction_pct")
add_num("trial_red_1", trial_red, 3,
        "Relative reduction of N_prop from rho = 0 to 0.8 (one decimal in %)",
        "tab4_real_example.rds: data$reduction_pct")
add_num("trial_rmax", table4_obj$R_max_pct / 100, 3,
        "R_max(1.292, 0.13)", "tab4_real_example.rds: R_max_pct")
add_num("trial_rmax_r0", table4_obj$R_max_pct / 100, 2,
        "R_max(1.292, 0.13), rounded to a whole percent",
        "tab4_real_example.rds: R_max_pct")
add_num("trial_recover", trial_red / (table4_obj$R_max_pct / 100), 2,
        "Integer reduction at rho = 0.8 relative to R_max",
        "tab4_real_example.rds")
lam4 <- table4_obj$lambda
add_num("trial_red_cont", lam4$reduction_cont_pct[lam4$rho == 0.8] / 100, 3,
        "Reduction at rho = 0.8 on the continuous scale of (eq:R_def)",
        "tab4_real_example.rds: lambda")
add_num("trial_kstar_05", table4_obj$kappa_star_eps05, 3,
        "kappa_star(0.05, 0.13)", "tab4_real_example.rds: kappa_star_eps05")
add_num("trial_kstar_01", table4_obj$kappa_star_eps01, 3,
        "kappa_star(0.01, 0.13)", "tab4_real_example.rds: kappa_star_eps01")
add_num("trial_kappa_round", table4_obj$kappa_round, 2,
        "kappa with the effect size 0.25 stated in the protocol",
        "tab4_real_example.rds: kappa_round")
add_num("trial_rmax_round", table4_obj$R_max_round_pct / 100, 3,
        "R_max(1.25, 0.13)", "tab4_real_example.rds: R_max_round_pct")
add_num("trial_rmax_round_r0", table4_obj$R_max_round_pct / 100, 2,
        "R_max(1.25, 0.13), rounded to a whole percent",
        "tab4_real_example.rds: R_max_round_pct")
add_num("trial_gain_bound",
        table4_obj$R_max_pct / 100 * t4$N_prop[t4$rho == 0], 0,
        "R_max(1.292, 0.13) times N_prop at rho = 0 (patients)",
        "tab4_real_example.rds")

# --------------------------------------------------------------------------- #
# Figure 3: reduction on the continuous scale at kappa = 1, 1 - beta = 0.80
# --------------------------------------------------------------------------- #
d3 <- reduction_rho_obj$data
red_at <- function(k, rr, pw) {
  d3$reduction[abs(d3$kappa - k) < 1e-9 & abs(d3$rho - rr) < 1e-9 &
                 abs(d3$power - pw) < 1e-9]
}
for (rr in c(0.3, 0.5, 0.8)) {
  add_num(sprintf("red_cont_k1_rho%02.0f_80", 10 * rr), red_at(1, rr, 0.80),
          3, sprintf(paste0("Reduction on the continuous scale, kappa = 1, ",
                            "rho = %.1f, 1 - beta = 0.80"), rr),
          "fig3_reduction_rho.rds: data")
}

ledger_df <- if (length(ledger) > 0L) do.call(rbind, ledger) else
  data.frame(key = character(0), value = numeric(0), digits = integer(0),
             text = character(0), description = character(0),
             source = character(0))
rownames(ledger_df) <- NULL
utils::write.csv(ledger_df, file.path(results_dir, "manuscript_numbers.csv"),
                 row.names = FALSE)
message("Ledger saved: ", file.path(results_dir, "manuscript_numbers.csv"),
        " (", nrow(ledger_df), " numbers)")
