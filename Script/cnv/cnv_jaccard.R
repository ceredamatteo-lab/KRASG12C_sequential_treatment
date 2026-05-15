#!/usr/bin/env Rscript
# =============================================================================
# cnv_jaccard.R
# Triangular Jaccard-similarity heatmaps of CNV gain (Amplification) and
# loss (Loss) gene sets across conditions {naive, sotorasib, adagrasib} for
# each cell line (H23, CALU, MEP). One PDF per cell line:
#
#   Figures/CNV/jaccard_<line>.pdf
#     - Two 3 × 3 lower-triangular heatmaps side-by-side: Amp | Loss
#     - Cells labelled with J(A,B) = |A ∩ B| / |A ∪ B|
#     - Axis labels show set size n in parentheses
#     - Upper triangle masked (Jaccard is symmetric)
#
# Sample → condition mapping uses the CNV table's `condition` column with
# the same naive/sotorasib/adagrasib labels as cnv_upset_heatmaps.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scico)
})

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
INPUT_DIR <- file.path(.project_root, "input")
OUT_DIR   <- file.path(.project_root, "Figures", "CNV")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

CONDS      <- c("000", "010", "020")
COND_LABEL <- c("000" = "naive",
                "010" = "sotorasib",
                "020" = "adagrasib")
COND_LEVELS <- unname(COND_LABEL[CONDS])
EVENTS     <- c(Amp = "Amplification", Loss = "Loss")

LINES <- list(
  H23  = file.path(INPUT_DIR, "h23_cnv_snv_sv.rds"),
  CALU = file.path(INPUT_DIR, "calu_cnv_snv_sv.rds"),
  MEP  = file.path(INPUT_DIR, "mep_cnv_snv_sv.rds")
)

# ---- helpers ----------------------------------------------------------------
gene_sets_for_line <- function(rds_path) {
  cnv <- as.data.frame(readRDS(rds_path)$CNV)
  cnv$condition <- as.character(cnv$condition)
  out <- list()
  for (ev_name in names(EVENTS)) {
    cat_str <- EVENTS[[ev_name]]
    sub <- cnv[cnv$category == cat_str & cnv$condition %in% CONDS, , drop = FALSE]
    by_cond <- split(sub$gene_name, sub$condition)
    by_cond <- lapply(by_cond, function(g) unique(g[!is.na(g) & nzchar(g)]))
    out[[ev_name]] <- setNames(by_cond[CONDS], COND_LEVELS)
  }
  out
}

jaccard_matrix <- function(set_list) {
  n <- length(set_list)
  m <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(names(set_list), names(set_list)))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    a <- set_list[[i]]; b <- set_list[[j]]
    u <- length(union(a, b))
    m[i, j] <- if (u == 0) NA_real_ else length(intersect(a, b)) / u
  }
  m
}

triangular_long <- function(mat, set_sizes) {
  df <- as.data.frame(mat) %>%
    tibble::rownames_to_column("row") %>%
    pivot_longer(-row, names_to = "col", values_to = "jaccard")
  ord <- COND_LEVELS
  df %>%
    mutate(row = factor(row, levels = ord),
           col = factor(col, levels = ord),
           # keep lower triangle + diagonal
           keep = as.integer(row) >= as.integer(col)) %>%
    filter(keep) %>%
    mutate(row_lab = sprintf("%s\n(n=%d)", row, set_sizes[as.character(row)]),
           col_lab = sprintf("%s\n(n=%d)", col, set_sizes[as.character(col)]),
           row_lab = factor(row_lab,
                            levels = sprintf("%s\n(n=%d)", ord, set_sizes[ord])),
           col_lab = factor(col_lab,
                            levels = sprintf("%s\n(n=%d)", ord, set_sizes[ord])))
}

plot_one <- function(mat, set_sizes, title) {
  d <- triangular_long(mat, set_sizes)
  ggplot(d, aes(x = col_lab, y = row_lab, fill = jaccard)) +
    geom_tile(colour = "white", size = 0.6) +
    geom_text(aes(label = sprintf("%.2f", jaccard)),
              size = 3, colour = "black") +
    scale_fill_gradientn(
      colours = scico::scico(11, palette = "lajolla"),
      limits = c(0, 1), na.value = "transparent",
      name = "Jaccard") +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(limits = rev) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid       = element_blank(),
          axis.text.x      = element_text(size = 8, lineheight = 0.85),
          axis.text.y      = element_text(size = 8, lineheight = 0.85),
          plot.title       = element_text(size = 10, face = "bold"),
          plot.title.position = "plot",
          legend.position  = "right",
          legend.key.height = unit(8, "mm"),
          legend.key.width  = unit(3, "mm"))
}

# ---- driver -----------------------------------------------------------------
for (line in names(LINES)) {
  message(sprintf("[%s] computing Jaccard matrices ...", line))
  sets <- gene_sets_for_line(LINES[[line]])

  panels <- list()
  for (ev_name in names(EVENTS)) {
    sl <- sets[[ev_name]]
    sizes <- vapply(sl, length, integer(1))
    if (sum(sizes) == 0) {
      message(sprintf("  [%s/%s] all sets empty — skipping", line, ev_name))
      next
    }
    m <- jaccard_matrix(sl)
    panels[[ev_name]] <- plot_one(
      m, set_sizes = sizes,
      title = sprintf("%s — %s",
                      line,
                      ifelse(ev_name == "Amp", "Amplification (gain)", "Loss")))
  }

  # Combined view: per condition, union of Amp ∪ Loss (any CNV-affected gene).
  combined <- mapply(union, sets$Amp, sets$Loss, SIMPLIFY = FALSE)
  combined <- lapply(combined, unique)
  combined_sizes <- vapply(combined, length, integer(1))
  if (sum(combined_sizes) > 0) {
    panels[["All"]] <- plot_one(
      jaccard_matrix(combined),
      set_sizes = combined_sizes,
      title = sprintf("%s — Amp ∪ Loss (any CNV)", line))
  }

  if (!length(panels)) next

  panel <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      title = sprintf("%s — Jaccard similarity of CNV gene sets across conditions",
                      line),
      theme = theme(plot.title = element_text(size = 11, face = "bold")))

  out_pdf <- file.path(OUT_DIR, sprintf("jaccard_%s.pdf", line))
  panel_w <- 3.6 * length(panels) - 0.5
  ggsave(out_pdf, panel, width = panel_w, height = 3.6, useDingbats = FALSE)
  message(sprintf("Wrote %s   (%.1f × 3.6 in)", out_pdf, panel_w))
}
