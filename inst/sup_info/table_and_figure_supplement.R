# =============================================================================
# table_and_figure_supplement.R
#
# Rendering stage for the Supporting Information of the twocpcont article.
# It reads inst/sup_info/data/*.rds and produces the figures, the LaTeX files
# and the ledger of the numbers cited in the Supporting Information, without
# any computation. The helpers are those of
# inst/scripts/table_and_figure_manuscript.R.
#
# Usage
#   Open twocpcont.Rproj and source this file after
#   data_generate_supplement.R:
#     source("inst/sup_info/table_and_figure_supplement.R")
#
# Required packages: ggplot2
#
# Outputs (relative to the package root)
#   inst/sup_info/results/appendix_software_example.tex   Section S5 (code
#                                   and output in a verbatim environment;
#                                   \input{appendix_software_example})
#   inst/sup_info/results/tableS1_cost_comparison.tex     Table S1
#   inst/sup_info/results/tableS2_accuracy_grid.tex       Table S2
#   inst/sup_info/results/tableS3_allocation.tex          Table S3
#   inst/sup_info/results/tableS4_planning_rho.tex        Table S4
#   inst/sup_info/results/FigureS1_error_decomposition.eps     Figure S1
#   inst/sup_info/results/FigureS2_accuracy_grid.eps           Figure S2
#   inst/sup_info/results/FigureS3_expedition_sensitivity.eps  Figure S3
#   inst/sup_info/results/supplement_numbers.csv   Ledger of cited numbers
#
# Each table file is a complete table environment (caption, label, tabular
# and table notes), included with \input{tableS1_cost_comparison}; the
# preamble of the Supporting Information must load threeparttable. Each
# figure is included with
# \includegraphics[width=\textwidth]{FigureS1_error_decomposition.eps}, and
# the captions of the figures are written in the Supporting Information.
# Captions refer to the main article by number (e.g., Table~2 of the main
# article), not by \ref.
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "twocpcont")) {
  stop("Run this script from the root of the twocpcont package ",
       "(open twocpcont.Rproj first).", call. = FALSE)
}

library(ggplot2)

data_dir    <- file.path("inst", "sup_info", "data")
results_dir <- file.path("inst", "sup_info", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
out_path <- function(name) file.path(results_dir, name)

BASE_SIZE <- 10
LW        <- 0.6
DPI       <- 800
FIG_WIDTH <- 160 / 25.4          # inches (160 mm)

# Save a figure as an 800 dpi EPS file (see table_and_figure_manuscript.R)
use_cairo <- isTRUE(capabilities("cairo"))
save_figure <- function(plot, name, width, height) {
  f_eps <- file.path(results_dir, paste0(name, ".eps"))
  if (use_cairo) {
    ggsave(f_eps, plot, device = grDevices::cairo_ps, width = width,
           height = height, units = "in", dpi = DPI,
           fallback_resolution = DPI, onefile = FALSE)
  } else {
    ggsave(f_eps, plot, device = "eps", width = width, height = height,
           units = "in", dpi = DPI)
  }
  message("Figure saved: ", f_eps)
}

theme_paper <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(size = 9))
}

# Okabe-Ito colours (colour-blind friendly)
OKABE_ITO <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2",
               "#D55E00", "#CC79A7", "#999999")

# The .tex files are written with LF line endings on every operating
# system (see table_and_figure_manuscript.R)
write_tex <- function(lines, filename) {
  con <- file(filename, open = "wb")
  on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
  message("LaTeX file saved: ", filename)
}

read_data <- function(name) {
  path <- file.path(data_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run data_generate_supplement.R first.",
                 path), call. = FALSE)
  }
  readRDS(path)
}

# Fixed-point formatting; NA becomes "---"
fmt <- function(x, digits = 3, na = "---") {
  ifelse(is.na(x), na, sprintf(paste0("%.", digits, "f"), x))
}

# Scientific notation in math mode, e.g. "$3.1 \times 10^{-4}$"; NA becomes
# "---" and zero "$0$"
fmt_sci <- function(x, digits = 1, na = "---") {
  vapply(x, function(v) {
    if (is.na(v)) return(na)
    if (v == 0) return("$0$")
    e <- floor(log10(abs(v)))
    m <- v / 10 ^ e
    if (round(abs(m), digits) >= 10) {
      e <- e + 1
      m <- v / 10 ^ e
    }
    sprintf(paste0("$%.", digits, "f \\times 10^{%d}$"), m, e)
  }, character(1L))
}

# A number with 'digits' significant digits (trailing zeros dropped)
fmt_sig <- function(x, digits = 3) {
  formatC(signif(x, digits), digits = digits, format = "fg")
}

# Format a vector of numbers as "$a$, $b$, and $c$" for captions
fmt_list <- function(x, fmt) {
  s <- paste0("$", sprintf(fmt, x), "$")
  if (length(s) <= 2L) return(paste(s, collapse = " and "))
  paste0(paste(s[-length(s)], collapse = ", "), ", and ", s[length(s)])
}

# Wrap a tabular in a complete table environment with caption, label and
# optional table notes (see table_and_figure_manuscript.R)
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

near <- function(x, y) abs(x - y) < 1e-9

# =============================================================================
# Section S5: console transcript of the software example
# =============================================================================
software_example <- read_data("software_example")
transcript <- software_example$transcript

# Guard against lines that would overflow the verbatim block
too_long <- nchar(transcript) > 76L
if (any(too_long)) {
  stop("Transcript lines longer than 76 characters:\n",
       paste(transcript[too_long], collapse = "\n"), call. = FALSE)
}

# Single spacing (setspace package), so that the transcript is compact
# in the double-spaced manuscript
tex <- c(paste0("% Generated by inst/sup_info/table_and_figure_supplement.R; ",
                "do not edit by hand."),
         "\\begin{singlespace}",
         "\\small",
         "\\begin{verbatim}",
         transcript,
         "\\end{verbatim}",
         "\\end{singlespace}")
write_tex(tex, file.path(results_dir, "appendix_software_example.tex"))

# =============================================================================
# Read the data of the other parts
# =============================================================================
cost_obj   <- read_data("cost_comparison")
decomp_obj <- read_data("error_decomposition")
grid_obj   <- read_data("accuracy_grid")
alloc_obj  <- read_data("allocation")
sens_obj   <- read_data("expedition_sensitivity")

# =============================================================================
# Table S1: comparison of the methods on the 96 combinations of Table 2 of
# the main article
# =============================================================================
pc <- cost_obj$data
n_conv <- pc[pc$method == "bivariate", c("cell_id", "N")]
names(n_conv)[2L] <- "N_conv"
pc <- merge(pc, n_conv, by = "cell_id", sort = FALSE)
pc <- pc[order(pc$cell_id, match(pc$method, unique(cost_obj$data$method))), ]

METHOD_ORDER  <- c("bivariate", "newton", "uniroot", "univariate")
METHOD_LABELS <- c(bivariate  = "Sequential search",
                   newton     = "Newton--Raphson algorithm",
                   uniroot    = "Brent's method (\\texttt{uniroot})",
                   univariate = "Proposed closed form")
timing <- cost_obj$timing
s1 <- do.call(rbind, lapply(METHOD_ORDER, function(m) {
  d <- pc[pc$method == m, ]
  data.frame(method      = m,
             n           = nrow(d),
             eval_median = stats::median(d$n_phi2),
             eval_min    = min(d$n_phi2),
             eval_max    = max(d$n_phi2),
             ms          = timing$ms_per_design[timing$method == m],
             reps        = timing$reps[timing$method == m],
             n_match     = sum(d$N == d$N_conv),
             rel_max     = if (all(is.na(d$rel_n2))) NA_real_ else
               max(abs(d$rel_n2)),
             stringsAsFactors = FALSE)
}))
n_cost <- unique(s1$n)
if (length(n_cost) != 1L) stop("Table S1: unequal numbers of designs.")
nr_iter <- pc$newton_iterations[pc$method == "newton"]

tabS1 <- c(
  "\\begin{tabular}{lrrrrr}",
  "\\hline",
  paste0(" & \\multicolumn{2}{c}{Evaluations of $\\Phi_2$} & ",
         "Time per design & & Largest relative \\\\"),
  "\\cline{2-3}",
  paste0("Method & Median & Range & (ms) & $N = N_{\\rm conv}$ & ",
         "error of $n_2$ \\\\"),
  "\\hline"
)
for (i in seq_len(nrow(s1))) {
  m <- s1$method[i]
  tabS1 <- c(tabS1, sprintf(
    "%s & %s & %d--%d & %s & %s & %s \\\\",
    METHOD_LABELS[[m]], fmt_sig(s1$eval_median[i], 3),
    as.integer(s1$eval_min[i]), as.integer(s1$eval_max[i]),
    fmt_sig(s1$ms[i], 3),
    if (m == "bivariate") "---" else sprintf("%d", s1$n_match[i]),
    fmt_sci(s1$rel_max[i])))
}
tabS1 <- c(tabS1, "\\hline", "\\end{tabular}")

ext <- cost_obj$extreme
fmt_set <- function(x) paste(sprintf("%g", sort(unique(x))), collapse = ", ")
capS1 <- sprintf(paste0(
  "Comparison of the sequential search, the Newton--Raphson algorithm of ",
  "Sugimoto et al.~\\cite{Sugimoto2012}, Brent's method and the proposed ",
  "closed form on the %d combinations of Table~2 of the main article. The ",
  "columns give the number of evaluations of the bivariate normal ",
  "distribution function $\\Phi_2$ per design, the computing time per ",
  "design, the number of combinations in which the total sample size $N$ ",
  "equals that of the sequential search ($N_{\\rm conv}$), and the largest ",
  "relative error of the continuous $n_2$ over the %d combinations. The ",
  "significance level is $\\alpha = %s$ (one-sided)."),
  n_cost, n_cost, format(cost_obj$alpha))
notesS1 <- c(
  sprintf(paste0(
    "\\item[] Evaluations of $\\Phi_2$ were counted by tracing ",
    "\\texttt{pbivnorm::pbivnorm()}, which the other three methods call. ",
    "The closed form evaluates no $\\Phi_2$; it evaluates the bivariate ",
    "normal density at the five Gauss--Legendre nodes. The Newton--Raphson ",
    "algorithm made ",
    "%d to %d updates; its stopping rule and step of the finite ",
    "differences are those of Sugimoto et al. Brent's method solves the ",
    "exact power equation for $\\lambda_w$ with tolerance $10^{%d}$."),
    as.integer(min(nr_iter)), as.integer(max(nr_iter)),
    as.integer(round(log10(cost_obj$uniroot_tol)))),
  sprintf(paste0(
    "\\item[] Times are elapsed times on one computer with \\texttt{R} ",
    "version %s, averaged over %s runs of the %d combinations (closed ",
    "form, Newton--Raphson algorithm, Brent's method, and sequential ",
    "search), and include the construction of the output of each method."),
    sub("^R version ([0-9.]+).*$", "\\1", cost_obj$meta$R_version),
    fmt_list(s1$reps[match(c("univariate", "newton", "uniroot", "bivariate"),
                           s1$method)], "%d"), n_cost),
  sprintf(paste0(
    "\\item[] The relative error is $|\\tilde n_2 / n_2 - 1|$, where ",
    "$\\tilde n_2$ is the value of $n_2$ before the ceiling operation and ",
    "$n_2 = c\\,\\lambda_w^2$ with $\\lambda_w$ the exact root of the power ",
    "equation. The sequential search has no continuous $n_2$.")),
  sprintf(paste0(
    "\\item[] At the %d combinations of $\\kappa \\in \\{%s\\}$, ",
    "$\\rho \\in \\{%s\\}$, and $1 - \\beta \\in \\{%s\\}$, the ",
    "Newton--Raphson algorithm met its stopping rule in %d, with %d to %d ",
    "updates, and the largest relative error of its continuous $n_2$ was ",
    "%s; that of the closed form was %s."),
    nrow(ext), fmt_set(ext$kappa), fmt_set(ext$rho), fmt_set(ext$power),
    sum(ext$newton_converged == 1), as.integer(min(ext$newton_iterations)),
    as.integer(max(ext$newton_iterations)),
    fmt_sci(max(abs(ext$newton_rel_n2))),
    fmt_sci(max(abs(ext$cf_rel_n2))))
)
write_table_env(tabS1, out_path("tableS1_cost_comparison.tex"),
                caption = capS1, label = "tab:cost_comparison",
                notes = notesS1, size = "\\footnotesize")

# =============================================================================
# Figure S1: error of the continuous n2 by approximation step
# =============================================================================
dec <- decomp_obj$data
STEPS <- c(anchor     = "Anchor",
           first      = "First-order correction",
           proposed   = "Proposed (quadratic correction)",
           gl_root    = "Quadrature error only",
           exact_quad = "Taylor error only")
ERR_FLOOR <- 1e-7      # |relative error| below this is drawn at this value
dl <- do.call(rbind, lapply(names(STEPS), function(v) {
  data.frame(kappa = dec$kappa, rho = dec$rho, power = dec$power, step = v,
             err = abs(dec[[paste0("rel_n2_", v)]]))
}))
dl$err_plot  <- pmax(dl$err, ERR_FLOOR)
dl$step      <- factor(dl$step, levels = names(STEPS))
dl$rho_lab   <- factor(sprintf("%g", dl$rho),
                       levels = sprintf("%g", sort(unique(dl$rho))))
dl$power_lab <- sprintf("%.2f", dl$power)

pS1 <- ggplot(dl, aes(x = kappa, y = err_plot, colour = step,
                      linetype = step)) +
  geom_line(linewidth = 0.5) +
  facet_grid(power_lab ~ rho_lab,
             labeller = label_bquote(rows = 1 - beta == .(power_lab),
                                     cols = rho == .(as.character(rho_lab)))) +
  scale_colour_manual(values = c(anchor = "#999999", first = "#E69F00",
                                 proposed = "#000000", gl_root = "#56B4E9",
                                 exact_quad = "#D55E00"),
                      labels = STEPS) +
  scale_linetype_manual(values = c(anchor = "solid", first = "dotdash",
                                   proposed = "solid", gl_root = "dashed",
                                   exact_quad = "longdash"),
                        labels = STEPS) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  scale_y_log10(limits = c(ERR_FLOOR, 1), breaks = 10 ^ c(-6, -4, -2, 0),
                labels = parse(text = c("10^-6", "10^-4", "10^-2", "1"))) +
  labs(x = expression(kappa),
       y = expression("|Relative error of continuous " * italic(n)[2] * "|")) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_paper() +
  theme(legend.text = element_text(size = 9),
        legend.key.width = unit(1.0, "cm"))
save_figure(pS1, "FigureS1_error_decomposition", width = FIG_WIDTH,
            height = 110 / 25.4)

# =============================================================================
# Table S2 and Figure S2: accuracy of the closed form over a wide grid
# =============================================================================
acc <- grid_obj$data
rho_set   <- sort(unique(acc$rho))
alpha_set <- sort(unique(acc$alpha))
power_set_grid <- sort(unique(acc$power))
kappa_grid <- sort(unique(acc$kappa))

tabS2 <- c(
  sprintf("\\begin{tabular}{r%s}",
          paste(rep("r", 3L + length(alpha_set)), collapse = "")),
  "\\hline",
  sprintf(paste0(" & & \\multicolumn{%d}{c}{Largest relative error of ",
                 "$n_2$} & Largest power & Smallest \\\\"),
          length(alpha_set)),
  sprintf("\\cline{3-%d}", 2L + length(alpha_set)),
  paste0("$\\rho$ & $\\max|\\lambda_w^{*} - \\lambda_w|$ & ",
         paste(sprintf("$\\alpha = %g$", alpha_set), collapse = " & "),
         " & difference & discriminant \\\\"),
  "\\hline"
)
for (rr in rho_set) {
  d <- acc[near(acc$rho, rr), ]
  rel_by_alpha <- vapply(alpha_set, function(a) {
    max(abs(d$rel_n2[near(d$alpha, a)]))
  }, numeric(1L))
  tabS2 <- c(tabS2, sprintf(
    "%s & %s & %s & %s & %s \\\\",
    sprintf("%g", rr), fmt_sci(max(abs(d$err_lambda))),
    paste(fmt_sci(rel_by_alpha), collapse = " & "),
    fmt_sci(max(abs(d$power_diff))), fmt_sig(min(d$disc), 2)))
}
tabS2 <- c(tabS2, "\\hline", "\\end{tabular}")

n_per_rho <- nrow(acc) / length(rho_set)
capS2 <- sprintf(paste0(
  "Accuracy of the closed form over $\\kappa$ from $%g$ to $%g$ in steps of ",
  "$%g$, target powers $1 - \\beta$ from $%.2f$ to $%.2f$ in steps of ",
  "$%.2f$, and one-sided significance levels $\\alpha$ of %s (%d designs ",
  "for each $\\rho$). The columns give the largest absolute error of the ",
  "noncentrality parameter $\\lambda_w^{*}$, the largest relative error ",
  "$|\\tilde n_2 / n_2 - 1|$ of the continuous $n_2$ for each $\\alpha$, ",
  "the largest absolute difference between the power at $\\tilde n_2$ and ",
  "the target power, and the smallest discriminant of the correction term."),
  min(kappa_grid), max(kappa_grid), diff(kappa_grid[1:2]),
  min(power_set_grid), max(power_set_grid), diff(power_set_grid[1:2]),
  fmt_list(alpha_set, "%g"), as.integer(n_per_rho))
notesS2 <- c(paste0(
  "\\item[] $\\lambda_w$ and $n_2 = c\\,\\lambda_w^2$ are the exact values ",
  "from the power equation with the bivariate normal distribution ",
  "function, and $\\tilde n_2 = c\\,(\\lambda_w^{*})^2$. The relative ",
  "error does not depend on the mean differences, the standard deviations ",
  "or the allocation ratio. The power at $\\tilde n_2$ is the bivariate ",
  "normal power at the noncentrality parameter $\\lambda_w^{*}$."))
write_table_env(tabS2, out_path("tableS2_accuracy_grid.tex"),
                caption = capS2, label = "tab:accuracy_grid",
                notes = notesS2, size = "\\footnotesize", tabcolsep = "4pt")

# Figure S2: signed relative error at alpha = 0.025, one panel per rho
ALPHA_FIG_S2 <- 0.025
dS2 <- acc[near(acc$alpha, ALPHA_FIG_S2), ]
dS2$rho_lab   <- factor(sprintf("%g", dS2$rho),
                        levels = sprintf("%g", rho_set))
dS2$power_lab <- factor(sprintf("%.2f", dS2$power),
                        levels = sprintf("%.2f", power_set_grid))
power_labels_S2 <- lapply(power_set_grid, function(p) {
  bquote(1 - beta == .(sprintf("%.2f", p)))
})
pS2 <- ggplot(dS2, aes(x = kappa, y = 100 * rel_n2, colour = power_lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "#888888",
             linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ rho_lab, ncol = 3, scales = "free_y",
             labeller = label_bquote(rho == .(as.character(rho_lab)))) +
  scale_colour_manual(values = OKABE_ITO[seq_along(power_set_grid)],
                      labels = power_labels_S2) +
  scale_x_continuous(breaks = seq(1, 3, by = 0.5)) +
  labs(x = expression(kappa),
       y = expression("Relative error of continuous " * italic(n)[2] *
                        " (%)")) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_paper() +
  theme(legend.text = element_text(size = 9))
save_figure(pS2, "FigureS2_accuracy_grid", width = FIG_WIDTH,
            height = 150 / 25.4)

# =============================================================================
# Table S3: allocation ratios other than one and two
# =============================================================================
al <- alloc_obj$data
r_set_al     <- sort(unique(al$r))
rho_set_al   <- sort(unique(al$rho))
power_set_al <- sort(unique(al$power))
al$round_diff <- al$N_prop != al$N_round_exact

tabS3 <- c(
  sprintf("\\begin{tabular}{ccr%s}",
          paste(rep("r", length(rho_set_al) * length(power_set_al)),
                collapse = "")),
  "\\hline",
  paste0(" & & & ", paste(sprintf(
    "\\multicolumn{%d}{c}{$1 - \\beta = %.2f$}", length(rho_set_al),
    power_set_al), collapse = " & "), " \\\\"),
  paste(sprintf("\\cline{%d-%d}",
                4L + length(rho_set_al) * (seq_along(power_set_al) - 1L),
                3L + length(rho_set_al) * seq_along(power_set_al)),
        collapse = " "),
  paste0("$(\\delta_1, \\delta_2)$ & $(\\sigma_1, \\sigma_2)$ & $r$ & ",
         paste(rep(sprintf("$\\rho = %s$", sprintf("%g", rho_set_al)),
                   length(power_set_al)), collapse = " & "), " \\\\"),
  "\\hline"
)
for (pt in sort(unique(al$pattern))) {
  for (j in seq_along(r_set_al)) {
    rr <- r_set_al[j]
    cells <- unlist(lapply(power_set_al, function(pw) {
      vapply(rho_set_al, function(rh) {
        d <- al[al$pattern == pt & near(al$r, rr) & near(al$rho, rh) &
                  near(al$power, pw), ]
        if (nrow(d) != 1L) stop("Table S3: cell not found or not unique.")
        sprintf("%d/%d%s", as.integer(d$N_conv), as.integer(d$N_prop),
                if (d$round_diff) "$^{\\dagger}$" else "")
      }, character(1L))
    }))
    p1 <- al[al$pattern == pt, ][1L, ]
    lead <- if (j == 1L) {
      sprintf("(%.1f, %.1f) & (%.1f, %.1f)", p1$delta1, p1$delta2, p1$sd1,
              p1$sd2)
    } else {
      " & "
    }
    tabS3 <- c(tabS3, sprintf("%s & %s & %s \\\\", lead, sprintf("%g", rr),
                              paste(cells, collapse = " & ")))
  }
}
tabS3 <- c(tabS3, "\\hline", "\\end{tabular}")

capS3 <- sprintf(paste0(
  "Total sample size of the sequential search and of the proposed closed ",
  "form, $N_{\\rm conv}/N_{\\rm prop}$, for allocation ratios $r$ of %s, ",
  "with $n_2 = \\lceil c\\,(\\lambda_w^{*})^2 \\rceil$ and ",
  "$n_1 = \\lceil r\\,n_2 \\rceil$. The significance level is ",
  "$\\alpha = %s$ (one-sided)."),
  fmt_list(r_set_al, "%g"), format(alloc_obj$alpha))
n_al <- nrow(al)
n_ne <- sum(al$N_prop != al$N_conv)
notesS3 <- c(
  sprintf(paste0(
    "\\item[] $N_{\\rm prop}$ differs from $N_{\\rm conv}$ in %d of the %d ",
    "combinations. Rounding the exact continuous $n_2 = c\\,\\lambda_w^2$ ",
    "in the same way gives $N_{\\rm prop}$ in %d of the %d combinations. ",
    "This rounded total differs from $N_{\\rm conv}$ in %d combinations, ",
    "because $n_1 = \\lceil r\\,n_2 \\rceil$ can exceed $r\\,n_2$ enough for ",
    "a smaller $n_2$ to attain the target power."),
    n_ne, n_al, sum(!al$round_diff), n_al,
    sum(al$N_round_exact != al$N_conv)),
  sprintf(paste0(
    "\\item[] The power at $N_{\\rm prop}$ exceeds the target by $%s$ to ",
    "$%s$. The Newton--Raphson algorithm gives $N_{\\rm conv}$ in %d of ",
    "the %d combinations."),
    fmt(min(al$power_prop - al$power), 5),
    fmt(max(al$power_prop - al$power), 5),
    sum(al$N_newton == al$N_conv), n_al),
  if (any(al$round_diff)) {
    paste0("\\item[$\\dagger$] $N_{\\rm prop}$ differs from the total ",
           "obtained by rounding the exact continuous $n_2$.")
  }
)
write_table_env(tabS3, out_path("tableS3_allocation.tex"),
                caption = capS3, label = "tab:allocation", notes = notesS3,
                size = "\\footnotesize", tabcolsep = "4pt")

# =============================================================================
# Table S4: power at the true correlation of the sample size planned at
# another correlation
# =============================================================================
pl <- sens_obj$planning
rho_true_set <- sort(unique(pl$rho_true))
rho_plan_set <- sort(unique(pl$rho_plan))
DESIGN_LABELS <- c(kappa1 = "A", expedition1 = "EXPEDITION~1")
designs_S4 <- unique(pl[, c("design", "delta1", "delta2", "sd1", "sd2",
                            "power")])
designs_S4 <- designs_S4[match(names(DESIGN_LABELS), designs_S4$design), ]

tabS4 <- c(
  sprintf("\\begin{tabular}{lrr%s}",
          paste(rep("r", length(rho_true_set)), collapse = "")),
  "\\hline",
  sprintf(paste0(" & & & \\multicolumn{%d}{c}{Power at the true ",
                 "correlation} \\\\"), length(rho_true_set)),
  sprintf("\\cline{4-%d}", 3L + length(rho_true_set)),
  paste0("Design & $\\rho_{\\rm plan}$ & $N$ & ",
         paste(sprintf("$\\rho = %s$", sprintf("%g", rho_true_set)),
               collapse = " & "), " \\\\"),
  "\\hline"
)
for (dsg in designs_S4$design) {
  for (j in seq_along(rho_plan_set)) {
    rp <- rho_plan_set[j]
    d <- pl[pl$design == dsg & near(pl$rho_plan, rp), ]
    d <- d[order(d$rho_true), ]
    if (nrow(d) != length(rho_true_set)) {
      stop("Table S4: rows not found or not unique.")
    }
    pw_cells <- sprintf("%.3f%s", d$power_true,
                        ifelse(d$power_true < d$power, "$^{\\dagger}$", ""))
    tabS4 <- c(tabS4, sprintf(
      "%s & %s & %d & %s \\\\",
      if (j == 1L) DESIGN_LABELS[[dsg]] else "", sprintf("%g", rp),
      as.integer(unique(d$N_plan)), paste(pw_cells, collapse = " & ")))
  }
  if (dsg != designs_S4$design[nrow(designs_S4)]) {
    tabS4 <- c(tabS4, "\\hline")
  }
}
tabS4 <- c(tabS4, "\\hline", "\\end{tabular}")

dA <- designs_S4[designs_S4$design == "kappa1", ]
dE <- designs_S4[designs_S4$design == "expedition1", ]
capS4 <- sprintf(paste0(
  "Co-primary power at the true correlation $\\rho$ of the total sample ",
  "size $N$ of the proposed closed form planned at the correlation ",
  "$\\rho_{\\rm plan}$. Design A has $(\\delta_1, \\delta_2) = ",
  "(%.1f, %.1f)$, $(\\sigma_1, \\sigma_2) = (%g, %g)$, and target power ",
  "$1 - \\beta = %.2f$; EXPEDITION~1 is the design of Table~4 of the main ",
  "article, with $1 - \\beta = %.2f$. Both have $r = 1$ and $\\alpha = %s$ ",
  "(one-sided)."),
  dA$delta1, dA$delta2, dA$sd1, dA$sd2, dA$power, dE$power,
  format(sens_obj$alpha))
notesS4 <- paste0("\\item[$\\dagger$] Power below the target power ",
                  "$1 - \\beta$ of the design.")
write_table_env(tabS4, out_path("tableS4_planning_rho.tex"),
                caption = capS4, label = "tab:planning_rho",
                notes = notesS4, size = "\\footnotesize")

# =============================================================================
# Figure S3: EXPEDITION 1, sample size and reduction over the standardized
# effect size of ADCS-ADL
# =============================================================================
sg   <- sens_obj$data
rmx  <- sens_obj$r_max
thr  <- sens_obj$thresholds
es1  <- sens_obj$es1
rho_set_S3 <- sort(unique(sg$rho))
SERIES_LEVELS <- c(sprintf("%g", rho_set_S3), "limit")
series_labels <- c(lapply(rho_set_S3, function(rr) {
  bquote(rho == .(sprintf("%g", rr)))
}), list(bquote(rho %->% 1 ~ (italic(R)[max]))))

# Panel N: total sample size; panel R: reduction relative to rho = 0 on the
# continuous scale, with its limit R_max as rho -> 1
dN <- data.frame(es2 = sg$es2, y = sg$N,
                 series = sprintf("%g", sg$rho), panel = "N")
dR <- sg[sg$rho > 0, ]
dR <- data.frame(es2 = dR$es2, y = 100 * dR$reduction_cont,
                 series = sprintf("%g", dR$rho), panel = "R")
dL <- data.frame(es2 = rmx$es2, y = 100 * rmx$R_max, series = "limit",
                 panel = "R")
dS3 <- rbind(dN, dR, dL)
dS3$series <- factor(dS3$series, levels = SERIES_LEVELS)
dS3$panel  <- factor(dS3$panel, levels = c("N", "R"))

# Effect size of ADCS-ADL at kappa_star(0.05, 0.13), drawn in panel R if it
# lies in the plotted range; the effect sizes of the protocol and of Table 4
# of the main article are drawn in both panels
es2_kstar <- es1 * thr$kappa_star[near(thr$epsilon, 0.05)]
es2_kstar <- es2_kstar[es2_kstar >= min(sg$es2) & es2_kstar <= max(sg$es2)]
vlines <- data.frame(panel = factor(rep("R", length(es2_kstar)),
                                    levels = c("N", "R")),
                     x = es2_kstar)
es2_design <- data.frame(x = c(sens_obj$es2_protocol, sens_obj$es2_design))

pS3 <- ggplot(dS3, aes(x = es2, y = y, colour = series, linetype = series)) +
  geom_vline(data = es2_design, aes(xintercept = x), colour = "#BBBBBB",
             linewidth = 0.4, linetype = "solid") +
  geom_vline(data = vlines, aes(xintercept = x), colour = "#888888",
             linewidth = 0.4, linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, scales = "free_y", strip.position = "left",
             labeller = as_labeller(c(N = "Total~sample~size~italic(N)",
                                      R = "Reduction~relative~to~rho==0~('%')"),
                                    label_parsed)) +
  scale_colour_manual(values = c(OKABE_ITO[seq_along(rho_set_S3)],
                                 "#D55E00"),
                      breaks = SERIES_LEVELS, labels = series_labels) +
  scale_linetype_manual(values = c(rep("solid", length(rho_set_S3)),
                                   "dashed"),
                        breaks = SERIES_LEVELS, labels = series_labels) +
  scale_x_continuous(breaks = seq(0.20, 0.30, by = 0.02)) +
  labs(x = "Standardized effect size of ADCS-ADL", y = NULL) +
  theme_paper() +
  theme(strip.placement  = "outside",
        strip.background = element_blank(),
        strip.text       = element_text(size = BASE_SIZE),
        legend.text      = element_text(size = 9))
save_figure(pS3, "FigureS3_expedition_sensitivity", width = FIG_WIDTH,
            height = 80 / 25.4)

# =============================================================================
# Ledger of the numbers cited in the Supporting Information:
# results/supplement_numbers.csv (same columns as manuscript_numbers.csv;
# keys must differ from those of the manuscript). Percentages are entered
# as fractions; number_locations.csv gives scale 100 where the text writes
# them in percent. The numbers of Section S5 are printed in the transcript
# itself, and the claims made about them in the text are checked by
# inst/validation/validate_claims_example.R.
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

# --------------------------------------------------------------------------- #
# Table S1: comparison of the methods
# --------------------------------------------------------------------------- #
s1_of <- function(m, v) s1[[v]][s1$method == m]
add_num("cost_n", n_cost, 0, "Combinations compared in Table S1",
        "cost_comparison.rds: data")
for (m in c("bivariate", "newton", "uniroot")) {
  nm <- c(bivariate = "conv", newton = "newton", uniroot = "uniroot")[[m]]
  add_num(sprintf("cost_eval_%s_median", nm), s1_of(m, "eval_median"), 1,
          sprintf("Median evaluations of Phi2 per design, %s", m),
          "cost_comparison.rds: data$n_phi2")
  add_num(sprintf("cost_eval_%s_min", nm), s1_of(m, "eval_min"), 0,
          sprintf("Smallest evaluations of Phi2 per design, %s", m),
          "cost_comparison.rds: data$n_phi2")
  add_num(sprintf("cost_eval_%s_max", nm), s1_of(m, "eval_max"), 0,
          sprintf("Largest evaluations of Phi2 per design, %s", m),
          "cost_comparison.rds: data$n_phi2")
}
add_num("cost_iter_newton_min", min(nr_iter), 0,
        "Fewest Newton-Raphson updates", "cost_comparison.rds: data")
add_num("cost_iter_newton_max", max(nr_iter), 0,
        "Most Newton-Raphson updates", "cost_comparison.rds: data")
for (m in METHOD_ORDER) {
  nm <- c(bivariate = "conv", newton = "newton", uniroot = "uniroot",
          univariate = "closed")[[m]]
  add_num(sprintf("cost_ms_%s", nm), s1_of(m, "ms"), 2,
          sprintf("Time per design (ms), %s", m),
          "cost_comparison.rds: timing")
}
add_num("cost_ratio_newton", s1_of("newton", "ms") / s1_of("univariate", "ms"),
        0, "Time of the Newton-Raphson algorithm relative to the closed form",
        "cost_comparison.rds: timing")
add_num("cost_ratio_uniroot",
        s1_of("uniroot", "ms") / s1_of("univariate", "ms"), 0,
        "Time of Brent's method relative to the closed form",
        "cost_comparison.rds: timing")
for (m in c("newton", "uniroot", "univariate")) {
  nm <- c(newton = "newton", uniroot = "uniroot", univariate = "closed")[[m]]
  add_num(sprintf("cost_match_%s", nm), s1_of(m, "n_match"), 0,
          sprintf("Combinations with N = N_conv, %s", m),
          "cost_comparison.rds: data")
}
add_num("cost_rel_closed", s1_of("univariate", "rel_max"), 4,
        "Largest relative error of the continuous n2, closed form",
        "cost_comparison.rds: data$rel_n2")
add_num("cost_rel_newton", s1_of("newton", "rel_max"), 6,
        "Largest relative error of the continuous n2, Newton-Raphson",
        "cost_comparison.rds: data$rel_n2")
add_num("cost_ext_n", nrow(ext), 0,
        "Combinations with extreme correlations",
        "cost_comparison.rds: extreme")
add_num("cost_ext_converged", sum(ext$newton_converged == 1), 0,
        "Extreme combinations in which the Newton-Raphson algorithm converged",
        "cost_comparison.rds: extreme")
add_num("cost_ext_iter_max", max(ext$newton_iterations), 0,
        "Most Newton-Raphson updates at the extreme correlations",
        "cost_comparison.rds: extreme")
add_num("cost_ext_rel_closed", max(abs(ext$cf_rel_n2)), 3,
        "Largest relative error of the closed form at the extreme correlations",
        "cost_comparison.rds: extreme")

# --------------------------------------------------------------------------- #
# Figure S1: error decomposition
# --------------------------------------------------------------------------- #
dec_digits <- c(anchor = 2, first = 3, proposed = 4, gl_root = 4,
                exact_quad = 4)
for (v in names(STEPS)) {
  add_num(sprintf("dec_rel_%s_max", v), max(abs(dec[[paste0("rel_n2_", v)]])),
          dec_digits[[v]],
          sprintf("Largest |relative error| of the continuous n2, %s",
                  STEPS[[v]]),
          sprintf("error_decomposition.rds: data$rel_n2_%s", v))
}
dec_lo <- dec[dec$rho <= 0.8, ]
add_num("dec_rel_proposed_max_lo", max(abs(dec_lo$rel_n2_proposed)), 4,
        "Largest |relative error| of the proposed n2 for rho <= 0.8",
        "error_decomposition.rds: data$rel_n2_proposed")
add_num("dec_rel_first_max_lo", max(abs(dec_lo$rel_n2_first)), 3,
        "Largest |relative error| of the first-order n2 for rho <= 0.8",
        "error_decomposition.rds: data$rel_n2_first")

# --------------------------------------------------------------------------- #
# Table S2: accuracy over a wide grid
# --------------------------------------------------------------------------- #
add_num("grid_n", nrow(acc), 0, "Designs of the accuracy grid",
        "accuracy_grid.rds: data")
acc_lo <- acc[acc$rho <= 0.9, ]
add_num("grid_rel_max_lo", max(abs(acc_lo$rel_n2)), 4,
        "Largest |relative error| of the continuous n2 for rho <= 0.9",
        "accuracy_grid.rds: data$rel_n2")
for (rr in c(0.95, 0.99)) {
  add_num(sprintf("grid_rel_max_rho%02.0f", 100 * rr),
          max(abs(acc$rel_n2[near(acc$rho, rr)])), 3,
          sprintf("Largest |relative error| of the continuous n2, rho = %g",
                  rr), "accuracy_grid.rds: data$rel_n2")
}
add_num("grid_pow_max_lo", max(abs(acc_lo$power_diff)), 5,
        "Largest |power at the closed-form n2 - target| for rho <= 0.9",
        "accuracy_grid.rds: data$power_diff")
add_num("grid_disc_min", min(acc$disc), 3,
        "Smallest discriminant of the correction term over the grid",
        "accuracy_grid.rds: data$disc")
add_num("grid_disc_neg", sum(acc$disc < 0), 0,
        "Designs of the grid with a negative discriminant",
        "accuracy_grid.rds: data$disc")

# --------------------------------------------------------------------------- #
# Table S3: allocation ratios
# --------------------------------------------------------------------------- #
add_num("alloc_n", n_al, 0, "Combinations of Table S3",
        "allocation.rds: data")
add_num("alloc_ne", n_ne, 0, "Combinations with N_prop != N_conv",
        "allocation.rds: data")
add_num("alloc_round_same", sum(!al$round_diff), 0,
        "Combinations with N_prop = N from the rounded exact continuous n2",
        "allocation.rds: data")
add_num("alloc_round_ne", sum(al$N_round_exact != al$N_conv), 0,
        "Combinations with N from the rounded exact continuous n2 != N_conv",
        "allocation.rds: data")
add_num("alloc_excess_max", max(al$N_prop - al$N_conv), 0,
        "Largest N_prop - N_conv", "allocation.rds: data")
add_num("alloc_pow_excess_min", min(al$power_prop - al$power), 5,
        "Smallest excess of the power at N_prop over the target",
        "allocation.rds: data")
add_num("alloc_match_newton", sum(al$N_newton == al$N_conv), 0,
        "Combinations with N_newton = N_conv", "allocation.rds: data")

# --------------------------------------------------------------------------- #
# Table S4 and Figure S3: EXPEDITION 1 and the planning correlation
# --------------------------------------------------------------------------- #
plan_of <- function(dsg, rp, rt) {
  pl$power_true[pl$design == dsg & near(pl$rho_plan, rp) &
                  near(pl$rho_true, rt)]
}
for (dsg in c("kappa1", "expedition1")) {
  nm <- c(kappa1 = "A", expedition1 = "E1")[[dsg]]
  add_num(sprintf("plan_%s_pow_p08_t0", nm), plan_of(dsg, 0.8, 0), 3,
          sprintf("Power at rho = 0 of N planned at rho = 0.8, design %s", nm),
          "expedition_sensitivity.rds: planning")
  add_num(sprintf("plan_%s_pow_p05_t0", nm), plan_of(dsg, 0.5, 0), 3,
          sprintf("Power at rho = 0 of N planned at rho = 0.5, design %s", nm),
          "expedition_sensitivity.rds: planning")
  add_num(sprintf("plan_%s_pow_p0_tm03", nm), plan_of(dsg, 0, -0.3), 3,
          sprintf("Power at rho = -0.3 of N planned at rho = 0, design %s",
                  nm), "expedition_sensitivity.rds: planning")
}
rmax_at <- function(e2) rmx$R_max[near(rmx$es2, e2)]
add_num("sens_rmax_es020", rmax_at(min(rmx$es2)), 3,
        "R_max at the smallest effect size of ADCS-ADL in Figure S3",
        "expedition_sensitivity.rds: r_max")
add_num("sens_rmax_es030", rmax_at(max(rmx$es2)), 3,
        "R_max at the largest effect size of ADCS-ADL in Figure S3",
        "expedition_sensitivity.rds: r_max")
add_num("sens_es2_kstar05", es1 * thr$kappa_star[near(thr$epsilon, 0.05)], 3,
        "Effect size of ADCS-ADL at kappa_star(0.05, 0.13)",
        "expedition_sensitivity.rds: thresholds")
add_num("sens_red_rho05_max",
        max(sg$reduction_cont[near(sg$rho, 0.5)]), 3,
        "Largest reduction at rho = 0.5 over the effect sizes of Figure S3",
        "expedition_sensitivity.rds: data$reduction_cont")

ledger_df <- do.call(rbind, ledger)
rownames(ledger_df) <- NULL
utils::write.csv(ledger_df, file.path(results_dir, "supplement_numbers.csv"),
                 row.names = FALSE)
message("Ledger saved: ", file.path(results_dir, "supplement_numbers.csv"),
        " (", nrow(ledger_df), " numbers)")
