#!/usr/bin/env Rscript
# =============================================================================
# summary_MAPK_ERK_MEM_MEP.R
# Mouse counterpart of summary_MAPK_ERK.R, for the MEF lines MEM and MEP.
#
# A4-friendly summary panels, one per gene list (MAPK / ERK / RAS84):
#
#   panel_<sig>_summary.pdf  —  2 cell-line heatmaps side-by-side  + count bar
#     - heatmap (one per cell line)  : top-N up + top-N down DEGs from the
#                                       signature in that cell line, sample-
#                                       level log2(CPM+1), row-z-scored.
#     - count bar (single)            : signature DEG counts per
#                                       (cell line × drug × direction),
#                                       coord-flipped, ≤ 4 cm × < 4 cm panels.
#
# Per-signature outputs land in Figures/RNA/MEM_MEP/summary/top{10,20,25}/.
# Sample → treatment from second digit (0=NAIVE, 1=sotorasib, 2=adagrasib).
# Mouse symbol mapping: RAS84 via ras84$CCLE_feature_id (already mouse);
# MAPK / ERK via title-case conversion of human symbols (ABCA7 → Abca7).
# Palette: scico::vik (11-stop, clipped at ±2.5 z-scores).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
  library(scico)
  library(ggh4x)
})

HEATMAP_RAMP <- colorRamp2(seq(-2.5, 2.5, length.out = 11),
                           scico::scico(11, palette = "vik"))

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
INPUT_DIR <- file.path(.project_root, "input")
OUT_DIR   <- file.path(.project_root, "Figures", "RNA_min", "MEM_MEP", "summary")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

old <- list.files(OUT_DIR, pattern = "\\.pdf$",
                  full.names = TRUE, recursive = TRUE)
if (length(old)) file.remove(old)

PADJ_MAX <- 0.01
LFC_MIN  <- 0.5
TOP_N_LIST <- c(20)        # minimal: top-20 only

PER_ROW_CM   <- 0.2
BODY_W_CM    <- 1.8
ROW_LABEL_PT <- 6

BAR_FACET_W_CM <- 3.8
BAR_FACET_H_CM <- 3.5
BAR_BASE_FONT  <- 8

DRUG_LABEL <- c(AMG510 = "sotorasib", MRTX = "adagrasib")
TREATMENT_LEVELS <- c("NAIVE","sotorasib","adagrasib")
TREATMENT_COL    <- c(NAIVE     = "grey60",
                      sotorasib = "#0072B2",
                      adagrasib = "#D55E00")
DIGIT_TO_TREAT   <- c("0" = "NAIVE", "1" = "sotorasib", "2" = "adagrasib")

strip_version <- function(x) sub("\\..*$", "", x)
to_mouse_sym  <- function(s) {
  s <- as.character(s); s <- s[!is.na(s) & nzchar(s)]
  paste0(substr(s, 1, 1), tolower(substr(s, 2, nchar(s))))
}

samples_for_line <- function(counts_df, line) {
  tag <- if (line == "MEM") "mem01" else if (line == "MEP") "mep01" else stop(line)
  cols <- grep(tag, colnames(counts_df), value = TRUE)
  code <- sub(".*-", "", cols)
  treat <- unname(DIGIT_TO_TREAT[substr(code, 2, 2)])
  data.frame(sample = cols, treatment = treat, cell_line = line,
             stringsAsFactors = FALSE) %>%
    mutate(treatment = factor(treatment, levels = TREATMENT_LEVELS)) %>%
    arrange(treatment, sample)
}

normalise_counts <- function(mat) {
  lib <- colSums(mat)
  cpm <- t(t(mat) / lib) * 1e6
  log2(cpm + 1)
}

# ---- inputs -----------------------------------------------------------------
counts_df <- readRDS(file.path(INPUT_DIR, "MEM_MEP_counts.rds"))
deg       <- readRDS(file.path(INPUT_DIR, "MEM_MEP_deg.rds"))

ras84    <- read.csv(file.path(INPUT_DIR,
                               "RAS_84_PhilipEast_et_al_NatCom_2022.csv"),
                     stringsAsFactors = FALSE)
mapk_erk <- readRDS(file.path(INPUT_DIR,
                              "mapk_erk_signature_East_NatCom_2022.rds"))

# ---- mouse symbol → ENSMUSG lookup -----------------------------------------
sym2ensg_df <- counts_df %>%
  transmute(gene_id_nv = strip_version(Geneid), gene_name) %>%
  filter(!is.na(gene_name), nzchar(gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE)
sym2ensg <- setNames(sym2ensg_df$gene_id_nv, sym2ensg_df$gene_name)

map_to_ensg <- function(syms) {
  ensg <- unname(sym2ensg[syms])
  ensg[!is.na(ensg)]
}

ras84_ensg <- map_to_ensg(unique(ras84$CCLE_feature_id))

mapk_human <- unique(mapk_erk$gene_symbol[
  mapk_erk$gs_name == "GOBP_MAPK_CASCADE"])
erk_human  <- unique(mapk_erk$gene_symbol[
  mapk_erk$gs_name == "GOBP_ERK1_AND_ERK2_CASCADE"])

mapk_ensg <- map_to_ensg(to_mouse_sym(mapk_human))
erk_ensg  <- map_to_ensg(to_mouse_sym(erk_human))

message(sprintf("RAS84    : %d input → %d ENSMUSG (lost %d)",
                length(unique(ras84$CCLE_feature_id)),
                length(ras84_ensg),
                length(unique(ras84$CCLE_feature_id)) - length(ras84_ensg)))
message(sprintf("MAPK     : %d input → %d ENSMUSG (lost %d)",
                length(mapk_human), length(mapk_ensg),
                length(mapk_human) - length(mapk_ensg)))
message(sprintf("ERK      : %d input → %d ENSMUSG (lost %d)",
                length(erk_human),  length(erk_ensg),
                length(erk_human)  - length(erk_ensg)))
message(sprintf("MAPK ∩ ERK: %d", length(intersect(mapk_ensg, erk_ensg))))

# Minimal variant: MAPK and ERK only — RAS84 dropped per project request.
SIG_LISTS <- list(
  MAPK  = list(label = "MAPK cascade", ensg = mapk_ensg,  tag = "MAPK"),
  ERK   = list(label = "ERK cascade",  ensg = erk_ensg,   tag = "ERK")
)

# ---- helpers (cloned from summary_MAPK_ERK.R) ------------------------------

top_per_line <- function(line, sig_ensg, n) {
  contrasts_keep <- grep(paste0("^", line, "__"),
                         unique(deg$contrast), value = TRUE)
  d <- deg %>%
    filter(contrast %in% contrasts_keep,
           gene_id_nv %in% sig_ensg,
           !is.na(padj),
           !is.na(log2FoldChange),
           padj < PADJ_MAX,
           abs(log2FoldChange) > LFC_MIN)
  rng <- d %>% group_by(gene_id_nv) %>%
    summarise(max_lfc = max(log2FoldChange),
              min_lfc = min(log2FoldChange),
              .groups = "drop")
  up   <- rng %>% filter(max_lfc >  0) %>% slice_max(max_lfc, n = n) %>%
    pull(gene_id_nv)
  down <- rng %>% filter(min_lfc <  0) %>% slice_min(min_lfc, n = n) %>%
    pull(gene_id_nv)
  list(up = up, down = down)
}

build_heatmap_one_line <- function(line, sig_ensg, sig_label, gene_name_lookup, n,
                                   shared_syms = NULL,
                                   shared_only = FALSE,
                                   shared_ensg = NULL) {
  picks <- top_per_line(line, sig_ensg, n)
  smpl  <- samples_for_line(counts_df, line)
  raw   <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  norm  <- normalise_counts(raw)

  keep <- intersect(rownames(norm), union(picks$up, picks$down))
  if (shared_only && !is.null(shared_ensg))
    keep <- intersect(keep, shared_ensg)
  if (length(keep) < 2) return(NULL)

  m <- norm[keep, smpl$sample, drop = FALSE]
  sym <- gene_name_lookup[rownames(m)]
  sym[is.na(sym) | !nzchar(sym)] <- rownames(m)[is.na(sym) | !nzchar(sym)]
  dup <- duplicated(sym)
  if (any(dup)) sym[dup] <- paste0(sym[dup], "_", seq_len(sum(dup)))
  rownames(m) <- sym

  m_scaled <- t(scale(t(m)))
  ramp     <- HEATMAP_RAMP

  top_anno <- HeatmapAnnotation(
    treatment = smpl$treatment,
    col = list(treatment = TREATMENT_COL),
    annotation_legend_param = list(
      title_gp  = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8)),
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_side = "left",
    height               = unit(4, "mm"),
    show_legend          = TRUE
  )

  # Left annotation: flags genes that also appear in the OTHER cell line's
  # top-N heatmap for the same signature (cross-line "shared" picks).
  # Suppressed when shared_only = TRUE (every row is shared by construction).
  left_anno <- NULL
  if (!shared_only && !is.null(shared_syms) && length(shared_syms) > 0) {
    in_both <- ifelse(rownames(m_scaled) %in% shared_syms, "shared", "unique")
    left_anno <- rowAnnotation(
      "in both" = factor(in_both, levels = c("shared", "unique")),
      col = list("in both" = c(shared = "#117733", unique = "grey88")),
      annotation_legend_param = list(
        "in both" = list(title    = "in both heatmaps",
                          title_gp = gpar(fontsize = 8),
                          labels_gp = gpar(fontsize = 8))),
      annotation_name_gp   = gpar(fontsize = 7),
      annotation_name_rot  = 90,
      annotation_name_side = "top",
      simple_anno_size     = unit(2.5, "mm"),
      show_legend          = TRUE
    )
  }

  body_h <- unit(nrow(m_scaled) * PER_ROW_CM, "cm")
  Heatmap(
    m_scaled,
    name        = "z-score",
    col         = ramp,
    height      = body_h,
    width       = unit(BODY_W_CM, "cm"),
    cluster_rows         = TRUE,
    cluster_columns      = FALSE,
    show_row_dend        = TRUE,
    show_column_dend     = FALSE,
    show_row_names       = TRUE,
    row_names_side       = "right",
    row_names_gp         = gpar(fontsize = ROW_LABEL_PT),
    show_column_names    = FALSE,
    column_title         = sprintf("%s — %s top%d ↑ + ↓", line, sig_label, n),
    column_title_gp      = gpar(fontsize = 9, fontface = "bold"),
    top_annotation       = top_anno,
    left_annotation      = if (!is.null(left_anno)) left_anno else NULL,
    heatmap_legend_param = list(
      title          = "row z-score",
      title_gp       = gpar(fontsize = 8),
      labels_gp      = gpar(fontsize = 8),
      title_position = "topcenter")
  )
}

build_count_bar <- function(sig_ensg, sig_label) {
  d <- deg %>%
    filter(gene_id_nv %in% sig_ensg,
           !is.na(padj),
           padj < PADJ_MAX,
           abs(log2FoldChange) > LFC_MIN) %>%
    separate(contrast, into = c("cell_line","drug_str"), sep = "__") %>%
    mutate(drug      = sub("_vs_NAIVE$", "", drug_str),
           drug      = unname(DRUG_LABEL[drug]),
           drug      = factor(drug,      levels = c("sotorasib","adagrasib")),
           cell_line = factor(cell_line, levels = c("MEM","MEP")),
           direction = factor(ifelse(log2FoldChange > 0, "up", "down"),
                              levels = c("up","down"))) %>%
    count(cell_line, drug, direction, .drop = FALSE) %>%
    group_by(cell_line, drug) %>%
    mutate(total = sum(n),
           pct   = n / total * 100) %>%
    ungroup()

  totals <- d %>% distinct(cell_line, drug, total)

  ggplot(d, aes(x = drug, y = pct, fill = direction)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    geom_text(aes(label = sprintf("%d (%.0f%%)", n, pct)),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 2.8, colour = "white", lineheight = 0.85) +
    geom_text(data = totals,
              aes(x = drug, y = 103, label = sprintf("n = %d", total)),
              inherit.aes = FALSE, size = 2.8, hjust = 0) +
    scale_fill_manual(values = c(up = "#CC0000", down = "#3333FF"),
                      name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 115),
                       breaks = seq(0, 100, by = 25),
                       expand = expansion(add = c(0, 0))) +
    coord_flip() +
    facet_wrap(~ cell_line, nrow = 1) +
    force_panelsizes(cols = unit(BAR_FACET_W_CM, "cm"),
                     rows = unit(BAR_FACET_H_CM, "cm")) +
    labs(y = sprintf("%% of %s signature DEGs", sig_label),
         x = NULL,
         title    = sprintf("%s pathway DEG composition", sig_label),
         subtitle = sprintf("padj < %g and |log2FC| > %g   |   denominator = signature DEGs; n at right = total, n inside = up / down",
                            PADJ_MAX, LFC_MIN)) +
    theme_bw(base_size = BAR_BASE_FONT) +
    theme(strip.text          = element_text(face = "bold",
                                              size = BAR_BASE_FONT),
          axis.text           = element_text(size = BAR_BASE_FONT),
          axis.title          = element_text(size = BAR_BASE_FONT),
          legend.text         = element_text(size = BAR_BASE_FONT),
          legend.title        = element_text(size = BAR_BASE_FONT),
          plot.title          = element_text(size = BAR_BASE_FONT + 1,
                                              face = "bold"),
          plot.subtitle       = element_text(size = BAR_BASE_FONT - 1),
          plot.title.position = "plot",
          legend.position     = "right",
          panel.grid.major.y  = element_blank())
}

# ---- driver -----------------------------------------------------------------

raw_full <- as.matrix(counts_df[, grep("caks-", colnames(counts_df))])
rownames(raw_full) <- strip_version(counts_df$Geneid)
gene_name_lookup <- setNames(counts_df$gene_name, rownames(raw_full))

for (n in TOP_N_LIST) {
  thr_dir <- file.path(OUT_DIR, sprintf("top%d", n))
  if (!dir.exists(thr_dir)) dir.create(thr_dir, recursive = TRUE)

  for (sig_name in names(SIG_LISTS)) {
    sig <- SIG_LISTS[[sig_name]]
    message(sprintf("\n=== %s — top%d ===", sig$label, n))

    # Cross-line shared genes (intersection of MEM's and MEP's top-N picks).
    picks_mem <- top_per_line("MEM", sig$ensg, n)
    picks_mep <- top_per_line("MEP", sig$ensg, n)
    shared_ensg <- intersect(union(picks_mem$up, picks_mem$down),
                             union(picks_mep$up, picks_mep$down))
    shared_syms <- unname(gene_name_lookup[shared_ensg])
    shared_syms <- shared_syms[!is.na(shared_syms) & nzchar(shared_syms)]
    message(sprintf("    shared (in both top-%d heatmaps): %d genes",
                    n, length(shared_syms)))

    hm_mem  <- build_heatmap_one_line("MEM", sig$ensg, sig$label, gene_name_lookup, n,
                                       shared_syms = shared_syms)
    hm_mep  <- build_heatmap_one_line("MEP", sig$ensg, sig$label, gene_name_lookup, n,
                                       shared_syms = shared_syms)
    p_bar   <- build_count_bar(sig$ensg, sig$label)

    body_to_in <- function(n_rows) (n_rows * PER_ROW_CM) / 2.54
    for (pair in list(list(hm = hm_mem, line = "MEM"),
                      list(hm = hm_mep, line = "MEP"))) {
      if (is.null(pair$hm)) next
      n_rows <- nrow(pair$hm@matrix)
      h <- body_to_in(n_rows) + 2.2
      out <- file.path(thr_dir,
                       sprintf("heatmap_%s_%s_compact.pdf",
                               sig$tag, pair$line))
      grDevices::pdf(out, width = 4.5, height = h, useDingbats = FALSE, bg = "white")
      draw(pair$hm, heatmap_legend_side = "right",
           annotation_legend_side = "right", merge_legend = TRUE)
      invisible(dev.off())
      message(sprintf("Wrote %s   (%.1f × %.1f in)", out, 4.5, h))
    }

    bar_out <- file.path(thr_dir, sprintf("barplot_%s_DEG_counts.pdf", sig$tag))
    ggsave(bar_out, p_bar, width = 5, height = 2.6, useDingbats = FALSE)
    message(sprintf("Wrote %s", bar_out))

    hm_mem_grob <- grid.grabExpr(draw(hm_mem, heatmap_legend_side = "right",
                                      annotation_legend_side = "right",
                                      merge_legend = TRUE))
    hm_mep_grob <- grid.grabExpr(draw(hm_mep, heatmap_legend_side = "right",
                                      annotation_legend_side = "right",
                                      merge_legend = TRUE))
    p_bar_grob  <- grid.grabExpr(print(p_bar))

    panel <- (patchwork::wrap_elements(full = hm_mem_grob) |
              patchwork::wrap_elements(full = hm_mep_grob)) /
             patchwork::wrap_elements(full = p_bar_grob) +
             patchwork::plot_layout(heights = c(3, 1)) +
             patchwork::plot_annotation(
               title = sprintf("%s pathway under sotorasib / adagrasib (top%d ↑ + ↓ per cell line)  —  MEM / MEP",
                               sig$label, n),
               theme = theme(plot.title = element_text(size = 11, face = "bold")))

    panel_out <- file.path(thr_dir, sprintf("panel_%s_summary.pdf", sig$tag))
    panel_w <- 9
    max_body_in <- (max(nrow(hm_mem@matrix), nrow(hm_mep@matrix)) *
                    PER_ROW_CM) / 2.54
    panel_h <- max_body_in + 4.6
    ggsave(panel_out, panel, width = panel_w, height = panel_h,
           useDingbats = FALSE, limitsize = FALSE)
    message(sprintf("Wrote %s   (%.1f × %.1f in)", panel_out, panel_w, panel_h))

    # ------------------------------------------------------------------
    # Shared-only variant: same heatmaps + bar but rows restricted to the
    # cross-line intersection.
    # ------------------------------------------------------------------
    if (length(shared_ensg) >= 2) {
      hm_mem_sh <- build_heatmap_one_line(
        "MEM", sig$ensg, sig$label, gene_name_lookup, n,
        shared_only = TRUE, shared_ensg = shared_ensg)
      hm_mep_sh <- build_heatmap_one_line(
        "MEP", sig$ensg, sig$label, gene_name_lookup, n,
        shared_only = TRUE, shared_ensg = shared_ensg)

      for (pair in list(list(hm = hm_mem_sh, line = "MEM"),
                        list(hm = hm_mep_sh, line = "MEP"))) {
        if (is.null(pair$hm)) next
        n_rows <- nrow(pair$hm@matrix)
        h <- body_to_in(n_rows) + 2.2
        out <- file.path(thr_dir,
                         sprintf("heatmap_%s_%s_compact_shared.pdf",
                                 sig$tag, pair$line))
        grDevices::pdf(out, width = 4.5, height = h, useDingbats = FALSE, bg = "white")
        draw(pair$hm, heatmap_legend_side = "right",
             annotation_legend_side = "right", merge_legend = TRUE)
        invisible(dev.off())
        message(sprintf("Wrote %s   (%.1f × %.1f in)", out, 4.5, h))
      }

      if (!is.null(hm_mem_sh) && !is.null(hm_mep_sh)) {
        hm_mem_sh_grob <- grid.grabExpr(draw(hm_mem_sh, heatmap_legend_side = "right",
                                              annotation_legend_side = "right",
                                              merge_legend = TRUE))
        hm_mep_sh_grob <- grid.grabExpr(draw(hm_mep_sh, heatmap_legend_side = "right",
                                              annotation_legend_side = "right",
                                              merge_legend = TRUE))
        p_bar_grob_sh <- grid.grabExpr(print(p_bar))
        panel_sh <- (patchwork::wrap_elements(full = hm_mem_sh_grob) |
                     patchwork::wrap_elements(full = hm_mep_sh_grob)) /
                    patchwork::wrap_elements(full = p_bar_grob_sh) +
                    patchwork::plot_layout(heights = c(3, 1)) +
                    patchwork::plot_annotation(
                      title = sprintf("%s — shared genes only (top%d ↑+↓ ∩ across MEM / MEP; %d genes)",
                                      sig$label, n, length(shared_ensg)),
                      theme = theme(plot.title = element_text(size = 11, face = "bold")))
        panel_sh_out <- file.path(thr_dir,
                                  sprintf("panel_%s_summary_shared.pdf", sig$tag))
        max_body_sh_in <- (max(nrow(hm_mem_sh@matrix),
                               nrow(hm_mep_sh@matrix)) * PER_ROW_CM) / 2.54
        panel_sh_h <- max_body_sh_in + 4.6
        ggsave(panel_sh_out, panel_sh, width = panel_w, height = panel_sh_h,
               useDingbats = FALSE, limitsize = FALSE)
        message(sprintf("Wrote %s   (%.1f × %.1f in)",
                        panel_sh_out, panel_w, panel_sh_h))
      }
    } else {
      message(sprintf("    shared set < 2 — skipping shared-only %s top%d",
                      sig$tag, n))
    }
  }
}

# ---- Standalone RAS-family heatmap (Kras / Hras / Nras / Mras, mouse) ------
# Single 4-row heatmap: both cell lines (MEM, MEP) side-by-side; left strip
# = log2FC vs naive per (cell line × drug) with "*" where the contrast
# clears padj<0.01 & |log2FC|>0.5.
CELL_LINE_COL <- c(MEM = "#117733", MEP = "#882255")
LFC_RAMP_RAS  <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

build_ras_family_heatmap_mouse <- function() {
  ras_syms <- c("Kras", "Hras", "Nras", "Mras")

  ras_lookup <- counts_df %>%
    filter(gene_name %in% ras_syms) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    transmute(gene_name, gene_id_nv = strip_version(Geneid))
  if (!nrow(ras_lookup)) return(NULL)
  ras_ensg <- setNames(ras_lookup$gene_id_nv, ras_lookup$gene_name)
  ras_ensg <- ras_ensg[ras_syms]
  ras_ensg <- ras_ensg[!is.na(ras_ensg)]

  smpl <- bind_rows(samples_for_line(counts_df, "MEM"),
                    samples_for_line(counts_df, "MEP")) %>%
    mutate(cell_line = factor(cell_line, levels = c("MEM","MEP"))) %>%
    arrange(cell_line, treatment, sample)

  raw <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  norm <- normalise_counts(raw)
  ras_ensg <- ras_ensg[ras_ensg %in% rownames(norm)]
  if (length(ras_ensg) < 1) return(NULL)

  m <- norm[ras_ensg, smpl$sample, drop = FALSE]
  rownames(m) <- names(ras_ensg)
  m_scaled <- t(scale(t(m)))

  contrasts_v <- c("MEM__AMG510_vs_NAIVE", "MEM__MRTX_vs_NAIVE",
                   "MEP__AMG510_vs_NAIVE", "MEP__MRTX_vs_NAIVE")
  contrast_labels <- c("MEM\nsoto", "MEM\nada", "MEP\nsoto", "MEP\nada")

  lfc_mat <- matrix(NA_real_, nrow = length(ras_ensg), ncol = length(contrasts_v),
                    dimnames = list(names(ras_ensg), contrast_labels))
  sig_mat <- matrix(FALSE, nrow = length(ras_ensg), ncol = length(contrasts_v),
                    dimnames = list(names(ras_ensg), contrast_labels))
  for (j in seq_along(contrasts_v)) {
    sub <- deg %>% filter(contrast == contrasts_v[j], gene_id_nv %in% ras_ensg)
    for (i in seq_along(ras_ensg)) {
      r <- sub %>% filter(gene_id_nv == ras_ensg[i])
      if (nrow(r) >= 1) {
        lfc_mat[i, j] <- r$log2FoldChange[1]
        sig_mat[i, j] <- !is.na(r$padj[1]) & r$padj[1] < 0.01 &
                        !is.na(r$log2FoldChange[1]) &
                        abs(r$log2FoldChange[1]) > 0.5
      }
    }
  }

  top_anno <- HeatmapAnnotation(
    cell_line = smpl$cell_line,
    treatment = smpl$treatment,
    col = list(cell_line = CELL_LINE_COL, treatment = TREATMENT_COL),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_side = "left",
    height = unit(7, "mm"),
    show_legend = TRUE
  )

  star_pt <- function(sig_col) ifelse(sig_col, "*", "")
  left_anno <- rowAnnotation(
    "MEM soto" = anno_simple(lfc_mat[, 1], col = LFC_RAMP_RAS,
                             pch = star_pt(sig_mat[, 1]),
                             pt_size = unit(3, "mm"),
                             gp = gpar(col = "grey60", lwd = 0.4)),
    "MEM ada"  = anno_simple(lfc_mat[, 2], col = LFC_RAMP_RAS,
                             pch = star_pt(sig_mat[, 2]),
                             pt_size = unit(3, "mm"),
                             gp = gpar(col = "grey60", lwd = 0.4)),
    "MEP soto" = anno_simple(lfc_mat[, 3], col = LFC_RAMP_RAS,
                             pch = star_pt(sig_mat[, 3]),
                             pt_size = unit(3, "mm"),
                             gp = gpar(col = "grey60", lwd = 0.4)),
    "MEP ada"  = anno_simple(lfc_mat[, 4], col = LFC_RAMP_RAS,
                             pch = star_pt(sig_mat[, 4]),
                             pt_size = unit(3, "mm"),
                             gp = gpar(col = "grey60", lwd = 0.4)),
    annotation_name_gp   = gpar(fontsize = 7),
    annotation_name_rot  = 60,
    annotation_name_side = "bottom",
    simple_anno_size     = unit(4, "mm"),
    gap                  = unit(0.6, "mm")
  )

  body_h <- unit(nrow(m_scaled) * 0.5, "cm")
  body_w <- unit(0.3 * ncol(m_scaled), "cm")

  hm <- Heatmap(
    m_scaled,
    name        = "z-score",
    col         = HEATMAP_RAMP,
    height      = body_h,
    width       = body_w,
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    show_row_dend        = FALSE,
    show_column_dend     = FALSE,
    show_row_names       = TRUE,
    row_names_side       = "left",
    row_names_gp         = gpar(fontsize = 9),
    show_column_names    = FALSE,
    column_split         = smpl$cell_line,
    column_title         = "RAS family — MEM / MEP expression and drug log2FC",
    column_title_gp      = gpar(fontsize = 10, fontface = "bold"),
    top_annotation       = top_anno,
    right_annotation     = left_anno,
    heatmap_legend_param = list(
      title          = "row z-score",
      title_gp       = gpar(fontsize = 8),
      labels_gp      = gpar(fontsize = 8),
      title_position = "topcenter")
  )

  lfc_legend <- Legend(col_fun = LFC_RAMP_RAS, title = "log2FC vs naive",
                       at = c(-3, 0, 3),
                       title_gp = gpar(fontsize = 8),
                       labels_gp = gpar(fontsize = 8))
  sig_legend <- Legend(labels = "padj<0.01 & |log2FC|>0.5",
                       title  = "significance",
                       type   = "points", pch = "*",
                       title_gp = gpar(fontsize = 8),
                       labels_gp = gpar(fontsize = 8),
                       background = "white")
  list(hm = hm, extra_legends = list(lfc_legend, sig_legend))
}

# Standalone RAS-family heatmap intentionally skipped in the _min variant.
