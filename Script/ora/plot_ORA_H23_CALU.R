#!/usr/bin/env Rscript
# =============================================================================
# plot_ORA_H23_CALU.R
# Per cell line, per ontology group, two ORA dot-plots:
#
#   ontology groups
#     kh : kegg + hallmarks   filter: p.adjust ≤ 0.1
#     go : go-bp              filter: p.adjust < 0.001 AND FoldEnrichment ≥ 2
#
#   plot variants
#     all     : de.status == "all"           (one column)
#     updown  : de.status %in% {"up","down"} (two columns: up | down)
#
#   x-axis variants
#     FoldEnrichment  (default)
#     richFactor      (extra "_richFactor" PDFs)
#
# Total per cell line: 2 groups × 2 variants × 2 x-axes = 8 PDFs.
# All significant pathways are shown (no top-N cap). Each facet panel is
# 4 cm wide; height is proportional to pathway count (≈ 0.4 cm per row,
# clamped to [MIN_PANEL_CM, MAX_PANEL_CM]). All theme text ≥ 8 pt.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(ggh4x)
  library(grid)
})

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
FIGURES_DIR   <- file.path(.project_root, "Figures")
OUT_DIR       <- file.path(FIGURES_DIR, "ORA")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Sweep stale PDFs from previous layouts so the directory holds only the
# current set of outputs.
old <- list.files(OUT_DIR, pattern = "\\.pdf$", full.names = TRUE)
if (length(old)) file.remove(old)

DESC_MAX_CHAR <- 55
PANEL_W_CM    <- 4         # forced panel width
PER_PATH_CM   <- 0.4       # height per pathway row at 8 pt font
MIN_PANEL_CM  <- 0.8       # don't go thinner than this
MAX_PANEL_CM  <- 30        # safety cap on per-panel height
MARGIN_H_CM   <- 6         # title + subtitle + legend + strip + axis title
MARGIN_W_CM   <- 14        # y-axis text + axis title + strip + legend slack
CM_PER_IN     <- 2.54

BASE_FONT   <- 9
AXIS_FONT   <- 8
STRIP_FONT  <- 8
LEGEND_FONT <- 8

CONTRAST_COL <- c("sotorasib" = "#0072B2",
                  "adagrasib" = "#D55E00")
DESTATUS_SHAPE <- c("up" = 24, "down" = 25, "all" = 21)
X_LABEL <- c(FoldEnrichment = "Fold enrichment (sqrt scale)",
             richFactor     = "Rich factor (sqrt scale; k / K = DEG ∩ set / set)")
CELL_LINE_LEVELS <- c("H23", "CALU")

# Sqrt x-axis everywhere: compresses occasional FE outliers (e.g. Ribosome)
# and expands the small-value range for richFactor (which is bounded [0, 1]
# and otherwise clusters near zero on a linear axis).
x_scale <- function(x_var) {
  scale_x_sqrt(name = X_LABEL[[x_var]])
}

ONT_GROUPS <- list(
  kh = list(
    ontologies = c("kegg", "hallmarks"),
    filter_expr = quote(p.adjust <= 0.1),
    label       = "kegg + hallmarks  |  p.adjust ≤ 0.1"
  ),
  go = list(
    ontologies  = c("go-bp"),
    filter_expr = quote(p.adjust < 0.001 & FoldEnrichment >= 2),
    label       = "go-bp  |  p.adjust < 0.001 & FE ≥ 2"
  )
)

# ---------- helpers ----------------------------------------------------------

relabel_contrast <- function(x) {
  x <- sub("^[^_]+__", "", x)
  x <- sub("AMG510", "sotorasib", x)
  x <- sub("MRTX",   "adagrasib", x)
  sub("_vs_NAIVE$", "", x)
}

prepare_data <- function(df, group, de_keep) {
  fexpr <- group$filter_expr
  df %>%
    filter(ontology %in% group$ontologies,
           !is.na(p.adjust)) %>%
    filter(!!fexpr) %>%
    filter(de.status %in% de_keep) %>%
    mutate(contrast    = relabel_contrast(contrast),
           contrast    = factor(contrast, levels = names(CONTRAST_COL)),
           ontology    = factor(ontology, levels = group$ontologies),
           de.status   = factor(de.status, levels = de_keep),
           Description = str_trunc(as.character(Description), DESC_MAX_CHAR))
}

# Within each ontology, sort pathways by best (smallest) p.adjust ascending.
order_by_padjust <- function(df) {
  df %>%
    group_by(ontology, Description) %>%
    mutate(.best_padj = min(p.adjust, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(ontology, .best_padj) %>%
    mutate(Description = fct_inorder(Description)) %>%
    select(-.best_padj)
}

# Per-ontology panel height (cm), bounded; returns a unit vector aligned to
# the ontology factor levels.
panel_heights_unit <- function(d) {
  per_ont <- d %>% distinct(ontology, Description) %>%
    count(ontology, .drop = FALSE) %>%
    mutate(h_cm = pmin(MAX_PANEL_CM,
                      pmax(MIN_PANEL_CM, n * PER_PATH_CM)))
  ord <- match(levels(d$ontology), as.character(per_ont$ontology))
  unit(per_ont$h_cm[ord], "cm")
}

# Total page (in) inferred from the data; just generous enough that the
# forced panels + axis labels + legend always fit.
page_dims <- function(d, n_cols) {
  per_ont <- d %>% distinct(ontology, Description) %>%
    count(ontology, .drop = FALSE) %>%
    mutate(h_cm = pmin(MAX_PANEL_CM,
                      pmax(MIN_PANEL_CM, n * PER_PATH_CM)))
  total_h_cm <- sum(per_ont$h_cm) + MARGIN_H_CM
  total_w_cm <- PANEL_W_CM * n_cols + MARGIN_W_CM
  c(width = total_w_cm / CM_PER_IN, height = total_h_cm / CM_PER_IN)
}

write_pdf <- function(p, path, w_in, h_in) {
  grDevices::pdf(path, width = w_in, height = h_in,
                 useDingbats = FALSE, bg = "white")
  print(p)
  invisible(dev.off())
  message(sprintf("Wrote %s   (%.1f × %.1f in)", path, w_in, h_in))
}

# ---------- plot builders ----------------------------------------------------

plot_all <- function(df, group, x_var, title_prefix) {
  d <- prepare_data(df, group, de_keep = "all")
  if (!nrow(d)) {
    message(sprintf("  [%s, all]  no rows pass the filter — skipping",
                    title_prefix))
    return(NULL)
  }
  d <- order_by_padjust(d)
  list(
    plot = ggplot(d, aes(x = .data[[x_var]], y = fct_rev(Description),
                         colour = contrast, size = Count)) +
      geom_point(shape = 21, fill = NA, stroke = 1.1, alpha = 0.95) +
      scale_size(range = c(2, 5), name = "Genes overlap") +
      scale_colour_manual(values = CONTRAST_COL, drop = FALSE,
                          name = "Contrast") +
      x_scale(x_var) +
      facet_grid(ontology ~ ., scales = "free_y") +
      force_panelsizes(rows = panel_heights_unit(d),
                       cols = unit(PANEL_W_CM, "cm")) +
      labs(title    = sprintf("%s — %s", title_prefix, group$label),
           subtitle = "de.status = all   |   all significant pathways",
           y = NULL) +
      theme_bw(base_size = BASE_FONT) +
      theme(panel.grid.minor   = element_blank(),
            strip.text         = element_text(face = "bold", size = STRIP_FONT),
            axis.text.y        = element_text(size = AXIS_FONT),
            axis.text.x        = element_text(size = AXIS_FONT),
            legend.text        = element_text(size = LEGEND_FONT),
            legend.title       = element_text(size = LEGEND_FONT),
            legend.position    = "bottom",
            plot.title.position = "plot"),
    dims = page_dims(d, n_cols = 1)
  )
}

plot_updown <- function(df, group, x_var, title_prefix) {
  d <- prepare_data(df, group, de_keep = c("up","down"))
  if (!nrow(d)) {
    message(sprintf("  [%s, updown]  no rows pass the filter — skipping",
                    title_prefix))
    return(NULL)
  }
  d <- order_by_padjust(d)
  list(
    plot = ggplot(d, aes(x = .data[[x_var]], y = fct_rev(Description),
                         colour = contrast, shape = de.status, size = Count)) +
      geom_point(stroke = 1.0, alpha = 0.95, fill = NA) +
      scale_size(range = c(2, 5), name = "Genes overlap") +
      scale_colour_manual(values = CONTRAST_COL, drop = FALSE,
                          name = "Contrast") +
      scale_shape_manual(values = DESTATUS_SHAPE,
                         labels = c("up" = "up (▲)", "down" = "down (▼)"),
                         name   = "DE direction") +
      x_scale(x_var) +
      facet_grid(ontology ~ de.status, scales = "free_y") +
      force_panelsizes(rows = panel_heights_unit(d),
                       cols = unit(PANEL_W_CM, "cm")) +
      labs(title    = sprintf("%s — %s", title_prefix, group$label),
           subtitle = "de.status: up (▲, left col) | down (▼, right col)   |   all significant pathways",
           y = NULL) +
      theme_bw(base_size = BASE_FONT) +
      theme(panel.grid.minor   = element_blank(),
            strip.text         = element_text(face = "bold", size = STRIP_FONT),
            axis.text.y        = element_text(size = AXIS_FONT),
            axis.text.x        = element_text(size = AXIS_FONT),
            legend.text        = element_text(size = LEGEND_FONT),
            legend.title       = element_text(size = LEGEND_FONT),
            legend.position    = "bottom",
            plot.title.position = "plot"),
    dims = page_dims(d, n_cols = 2)
  )
}

# ---------- cross-line comparison plot --------------------------------------
# kh group, de.status == "all", H23 vs CALU as columns. Pathway set per
# ontology = union of significant pathways across both cell lines, with the
# same y-ordering applied to both columns so a missing dot in one column
# means "not significant in that line" — making cross-line similarity /
# uniqueness visible at a glance.
plot_cross_line <- function(combined_df, group, x_var) {
  fexpr <- group$filter_expr
  d <- combined_df %>%
    filter(ontology %in% group$ontologies,
           de.status == "all",
           !is.na(p.adjust)) %>%
    filter(!!fexpr) %>%
    mutate(contrast    = relabel_contrast(contrast),
           contrast    = factor(contrast,  levels = names(CONTRAST_COL)),
           ontology    = factor(ontology,  levels = group$ontologies),
           cell_line   = factor(cell_line, levels = CELL_LINE_LEVELS),
           Description = str_trunc(as.character(Description), DESC_MAX_CHAR))
  if (!nrow(d)) return(NULL)
  d <- order_by_padjust(d)

  list(
    plot = ggplot(d, aes(x = .data[[x_var]], y = fct_rev(Description),
                         colour = contrast, size = Count)) +
      geom_point(shape = 21, fill = NA, stroke = 1.1, alpha = 0.95) +
      scale_size(range = c(2, 5), name = "Genes overlap") +
      scale_colour_manual(values = CONTRAST_COL, drop = FALSE,
                          name = "Contrast") +
      x_scale(x_var) +
      facet_grid(ontology ~ cell_line, scales = "free_y") +
      force_panelsizes(rows = panel_heights_unit(d),
                       cols = unit(PANEL_W_CM, "cm")) +
      labs(title    = "ORA: H23 vs CALU — pathway similarity",
           subtitle = sprintf("%s   |   de.status = all   |   one column per cell line; missing dot = not significant in that line",
                              group$label),
           y = NULL) +
      theme_bw(base_size = BASE_FONT) +
      theme(panel.grid.minor   = element_blank(),
            strip.text         = element_text(face = "bold", size = STRIP_FONT),
            axis.text.y        = element_text(size = AXIS_FONT),
            axis.text.x        = element_text(size = AXIS_FONT),
            legend.text        = element_text(size = LEGEND_FONT),
            legend.title       = element_text(size = LEGEND_FONT),
            legend.position    = "bottom",
            plot.title.position = "plot"),
    dims = page_dims(d, n_cols = length(CELL_LINE_LEVELS))
  )
}

# Variant of the cross-line plot for "Nature-style" enrichment dotplots:
# colour-fill encodes -log10(p.adjust) (so darker = more significant), shape
# encodes the contrast (sotorasib ●, adagrasib ▲). Same row/col facetting and
# panel-size logic as plot_cross_line().
plot_cross_line_padj <- function(combined_df, group, x_var) {
  fexpr <- group$filter_expr
  d <- combined_df %>%
    filter(ontology %in% group$ontologies,
           de.status == "all",
           !is.na(p.adjust)) %>%
    filter(!!fexpr) %>%
    mutate(contrast      = relabel_contrast(contrast),
           contrast      = factor(contrast,  levels = names(CONTRAST_COL)),
           ontology      = factor(ontology,  levels = group$ontologies),
           cell_line     = factor(cell_line, levels = CELL_LINE_LEVELS),
           Description   = str_trunc(as.character(Description), DESC_MAX_CHAR),
           neglog10_padj = -log10(p.adjust))
  if (!nrow(d)) return(NULL)
  d <- order_by_padjust(d)

  list(
    plot = ggplot(d, aes(x = .data[[x_var]], y = fct_rev(Description),
                         fill = neglog10_padj, size = Count, shape = contrast)) +
      geom_point(stroke = 0.4, colour = "grey25", alpha = 0.95) +
      scale_size(range = c(2, 5), name = "Genes overlap") +
      scale_fill_viridis_c(option = "plasma", direction = -1,
                           name = "-log10(p.adjust)") +
      scale_shape_manual(values = c("sotorasib" = 21, "adagrasib" = 24),
                         name = "Contrast") +
      x_scale(x_var) +
      facet_grid(ontology ~ cell_line, scales = "free_y") +
      force_panelsizes(rows = panel_heights_unit(d),
                       cols = unit(PANEL_W_CM, "cm")) +
      labs(title    = "ORA: H23 vs CALU — pathway similarity (p.adjust shown)",
           subtitle = sprintf("%s   |   de.status = all   |   fill: -log10(p.adjust)   |   shape: contrast",
                              group$label),
           y = NULL) +
      theme_bw(base_size = BASE_FONT) +
      theme(panel.grid.minor   = element_blank(),
            strip.text         = element_text(face = "bold", size = STRIP_FONT),
            axis.text.y        = element_text(size = AXIS_FONT),
            axis.text.x        = element_text(size = AXIS_FONT),
            legend.text        = element_text(size = LEGEND_FONT),
            legend.title       = element_text(size = LEGEND_FONT),
            legend.position    = "bottom",
            plot.title.position = "plot") +
      guides(shape = guide_legend(override.aes = list(fill = "grey60", size = 4)),
             size  = guide_legend(override.aes = list(shape = 21, fill = "grey60"))),
    dims = page_dims(d, n_cols = length(CELL_LINE_LEVELS))
  )
}

# ---------- driver -----------------------------------------------------------

ora_h23  <- readRDS(file.path(.project_root, "ora_H23.rds"))
ora_calu <- readRDS(file.path(.project_root, "ora_CALU.rds"))

cell_lines <- list(H23 = ora_h23, CALU = ora_calu)

for (line in names(cell_lines)) {
  df <- cell_lines[[line]]
  for (gname in names(ONT_GROUPS)) {
    grp <- ONT_GROUPS[[gname]]
    title_prefix <- sprintf("%s — ORA: sotorasib vs adagrasib", line)
    for (x_var in c("FoldEnrichment","richFactor")) {
      x_suffix <- if (x_var == "richFactor") "_richFactor" else ""
      pa <- plot_all   (df, grp, x_var = x_var, title_prefix = title_prefix)
      pu <- plot_updown(df, grp, x_var = x_var, title_prefix = title_prefix)
      if (!is.null(pa)) write_pdf(pa$plot,
        file.path(OUT_DIR, sprintf("ora_%s_%s_all%s.pdf",    line, gname, x_suffix)),
        pa$dims["width"], pa$dims["height"])
      if (!is.null(pu)) write_pdf(pu$plot,
        file.path(OUT_DIR, sprintf("ora_%s_%s_updown%s.pdf", line, gname, x_suffix)),
        pu$dims["width"], pu$dims["height"])
    }
  }
}

# Cross-cell-line pathway similarity, kh group, de.status = "all"
combined <- bind_rows(
  ora_h23  %>% mutate(cell_line = "H23"),
  ora_calu %>% mutate(cell_line = "CALU")
)
for (x_var in c("FoldEnrichment","richFactor")) {
  x_suffix <- if (x_var == "richFactor") "_richFactor" else ""
  pc <- plot_cross_line(combined, ONT_GROUPS$kh, x_var = x_var)
  if (!is.null(pc)) write_pdf(pc$plot,
    file.path(OUT_DIR, sprintf("ora_kh_all_byCellLine%s.pdf", x_suffix)),
    pc$dims["width"], pc$dims["height"])

  pcp <- plot_cross_line_padj(combined, ONT_GROUPS$kh, x_var = x_var)
  if (!is.null(pcp)) write_pdf(pcp$plot,
    file.path(OUT_DIR, sprintf("ora_kh_all_byCellLine_padj%s.pdf", x_suffix)),
    pcp$dims["width"], pcp$dims["height"])
}
