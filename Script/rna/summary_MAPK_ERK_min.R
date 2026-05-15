#!/usr/bin/env Rscript
# =============================================================================
# summary_MAPK_ERK.R
# A4-friendly summary panels of pathway regulation under sotorasib /
# adagrasib in H23 and CALU. One panel per gene list (MAPK/ERK and RAS84):
#
#   panel_<sig>_summary.pdf  —  2 cell-line heatmaps side-by-side  + count bar
#     - heatmap (one per cell line)  : top-10 up + top-10 down DEGs from the
#                                       signature in that cell line, sample-
#                                       level log2(CPM+1), row-z-scored.
#     - count bar (single)            : MAPK/ERK or RAS84 DEG counts per
#                                       (cell line × drug × direction).
#
# Per-signature outputs in Figures/RNA/summary/:
#   panel_MAPK_ERK_summary.pdf
#   panel_RAS84_summary.pdf
#   heatmap_MAPK_ERK_H23_compact.pdf, heatmap_MAPK_ERK_CALU_compact.pdf
#   heatmap_RAS84_H23_compact.pdf,    heatmap_RAS84_CALU_compact.pdf
#   barplot_MAPK_ERK_DEG_counts.pdf,  barplot_RAS84_DEG_counts.pdf
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

# Bar-plot panel size — coord-flipped so bars are horizontal. Each facet
# locked to ≤4 cm wide × <4 cm tall.
BAR_FACET_W_CM <- 3.8
BAR_FACET_H_CM <- 3.5
BAR_BASE_FONT  <- 8

# 11-stop perceptually-uniform diverging palette (scico "vik"), clipped at
# ±2.5 z-scores so the midpoint sits cleanly at 0.
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
OUT_DIR   <- file.path(.project_root, "Figures", "RNA_min", "summary")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Sweep stale outputs from earlier layouts (root + any old top* subfolders)
old <- list.files(OUT_DIR, pattern = "\\.pdf$",
                  full.names = TRUE, recursive = TRUE)
if (length(old)) file.remove(old)

PADJ_MAX <- 0.01
LFC_MIN  <- 0.5
TOP_N_LIST <- c(20)        # minimal: top-20 only

# Pin per-row height so top20 (≤ 40 rows) caps at 8 cm. The other top-N
# variants scale proportionally — top10 ≤ 4 cm, top25 ≤ 10 cm.
PER_ROW_CM <- 0.2
BODY_W_CM    <- 1.8  # heatmap body width — kept under 2 cm
ROW_LABEL_PT <- 6   # small font is required to fit 0.2 cm rows without overlap

DRUG_LABEL <- c(AMG510 = "sotorasib", MRTX = "adagrasib")
TREATMENT_LEVELS <- c("NAIVE", "sotorasib", "adagrasib")
TREATMENT_COL    <- c(NAIVE     = "grey60",
                      sotorasib = "#0072B2",
                      adagrasib = "#D55E00")
DIGIT_TO_TREAT   <- c("0" = "NAIVE", "1" = "sotorasib", "2" = "adagrasib")

# CNV right-side annotation: Amp/Loss event per (gene × condition) per cell line.
CNV_FILES <- list(
  H23  = file.path(INPUT_DIR, "h23_cnv_snv_sv.rds"),
  CALU = file.path(INPUT_DIR, "calu_cnv_snv_sv.rds")
)
CNV_CONDS      <- c("000", "010", "020")
CNV_COND_LABEL <- c("000" = "naive", "010" = "sotorasib", "020" = "adagrasib")
# Continuous CNt colour ramp for the right-side annotation: blue at 0
# (loss), white at 2 (neutral), red at 6 (gain) — capped at 6 to match the
# drug-only CNV heatmaps. NA (no CNV row) renders white.
CNV_CNT_RAMP <- colorRamp2(c(0, 2, 6), c("#2166AC", "white", "#B2182B"))

strip_version <- function(x) sub("\\..*$", "", x)

# Per-line CNV CNt table: gene × condition → CNt (integer, max-magnitude
# from {Amp, Loss} rows in the CNV table). Cached across topN × signature
# loops so the rds is read once per cell line.
.cnv_cnt_cache <- list()
cnv_cnt_table_cached <- function(line) {
  if (!is.null(.cnv_cnt_cache[[line]])) return(.cnv_cnt_cache[[line]])
  if (!line %in% names(CNV_FILES)) {
    .cnv_cnt_cache[[line]] <<- NULL; return(NULL)
  }
  cnv <- as.data.frame(readRDS(CNV_FILES[[line]])$CNV)
  ev <- cnv %>%
    transmute(gene_name = as.character(gene_name),
              condition = as.character(condition),
              category  = as.character(category),
              CNt       = as.integer(CNt)) %>%
    filter(!is.na(gene_name), nzchar(gene_name),
           condition %in% CNV_CONDS,
           category %in% c("Amplification", "Loss")) %>%
    group_by(gene_name, condition) %>%
    summarise(
      # CNt with largest deviation from 2 (same rule as cnv_upset_heatmaps.R)
      CNt = CNt[which.max(abs(CNt - 2))],
      .groups = "drop")
  .cnv_cnt_cache[[line]] <<- ev
  ev
}

# Build a (genes × 3 conds) numeric CNt matrix; NA where no Amp/Loss row
# exists for that gene × condition. Returns NULL when no CNV file exists
# for the cell line (e.g. mouse MEM).
cnv_cnt_matrix <- function(line, genes) {
  ev <- cnv_cnt_table_cached(line)
  if (is.null(ev)) return(NULL)
  m <- matrix(NA_real_, nrow = length(genes), ncol = length(CNV_CONDS),
              dimnames = list(genes, unname(CNV_COND_LABEL[CNV_CONDS])))
  for (k in seq_along(CNV_CONDS)) {
    sub <- ev %>% filter(condition == CNV_CONDS[k], gene_name %in% genes)
    m[sub$gene_name, k] <- as.numeric(sub$CNt)
  }
  m
}

samples_for_line <- function(counts_df, line) {
  tag <- if (line == "H23") "h2301" else if (line == "CALU") "cal01" else stop(line)
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
counts_df <- readRDS(file.path(INPUT_DIR, "H23_CALU_counts.rds"))
deg       <- readRDS(file.path(INPUT_DIR, "H23_CALU_deg.rds"))

ras84    <- read.csv(file.path(INPUT_DIR,
                               "RAS_84_PhilipEast_et_al_NatCom_2022.csv"),
                     stringsAsFactors = FALSE)
mapk_erk <- readRDS(file.path(INPUT_DIR,
                              "mapk_erk_signature_East_NatCom_2022.rds"))

ras84_ensg <- unique(ras84$Uppsala_feature_ids)
ras84_ensg <- ras84_ensg[!is.na(ras84_ensg) & nzchar(ras84_ensg)]

# Split the East 2022 MAPK/ERK signature back into its two source GO sets.
mapk_ensg <- mapk_erk %>%
  filter(gs_name == "GOBP_MAPK_CASCADE") %>%
  pull(human_ensembl_gene) %>% unique()
mapk_ensg <- mapk_ensg[!is.na(mapk_ensg) & nzchar(mapk_ensg)]

erk_ensg <- mapk_erk %>%
  filter(gs_name == "GOBP_ERK1_AND_ERK2_CASCADE") %>%
  pull(human_ensembl_gene) %>% unique()
erk_ensg <- erk_ensg[!is.na(erk_ensg) & nzchar(erk_ensg)]

message(sprintf("MAPK cascade genes  : %d", length(mapk_ensg)))
message(sprintf("ERK cascade genes   : %d", length(erk_ensg)))
message(sprintf("MAPK ∩ ERK overlap  : %d",
                length(intersect(mapk_ensg, erk_ensg))))
message(sprintf("RAS84 genes         : %d", length(ras84_ensg)))

# Each pathway gets its own heatmap; overlapping genes (ERK ⊂ MAPK in
# practice) appear in both heatmaps so each panel shows its full pathway.
# Minimal variant: MAPK and ERK only — RAS84 dropped per project request.
SIG_LISTS <- list(
  MAPK  = list(label = "MAPK cascade", ensg = mapk_ensg,  tag = "MAPK"),
  ERK   = list(label = "ERK cascade",  ensg = erk_ensg,   tag = "ERK")
)

# ---- shared helpers ---------------------------------------------------------

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

CELL_LINE_COL <- c(H23 = "#117733", CALU = "#882255")
LFC_RAMP_RAS  <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

# Build a unified RAS-family heatmap: rows = {KRAS, HRAS, NRAS, MRAS},
# columns = all 18 samples (9 H23 + 9 CALU), z-scored row-wise across all
# samples. Left annotation = log2FC under each (cell line × drug vs naive)
# with an asterisk when the contrast clears padj<0.01 & |log2FC|>0.5.
build_ras_family_heatmap <- function(deg, counts_df, gene_name_lookup) {
  ras_syms <- c("KRAS", "HRAS", "NRAS", "MRAS")

  # Resolve symbols → ENSG via the counts gene_name column
  ras_lookup <- counts_df %>%
    filter(gene_name %in% ras_syms) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    transmute(gene_name, gene_id_nv = strip_version(Geneid))
  if (!nrow(ras_lookup)) return(NULL)
  ras_ensg <- setNames(ras_lookup$gene_id_nv, ras_lookup$gene_name)
  ras_ensg <- ras_ensg[ras_syms]   # preserve KRAS, HRAS, NRAS, MRAS order
  ras_ensg <- ras_ensg[!is.na(ras_ensg)]

  smpl <- bind_rows(samples_for_line(counts_df, "H23"),
                    samples_for_line(counts_df, "CALU")) %>%
    mutate(cell_line = factor(cell_line, levels = c("H23","CALU"))) %>%
    arrange(cell_line, treatment, sample)

  raw <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  norm <- normalise_counts(raw)
  if (!all(ras_ensg %in% rownames(norm))) {
    miss <- names(ras_ensg)[!ras_ensg %in% rownames(norm)]
    message("RAS-family genes missing from counts: ",
            paste(miss, collapse = ", "))
    ras_ensg <- ras_ensg[ras_ensg %in% rownames(norm)]
  }
  if (length(ras_ensg) < 1) return(NULL)

  m <- norm[ras_ensg, smpl$sample, drop = FALSE]
  rownames(m) <- names(ras_ensg)
  m_scaled <- t(scale(t(m)))

  # log2FC + significance per gene × contrast
  contrasts_v <- c("H23__AMG510_vs_NAIVE",  "H23__MRTX_vs_NAIVE",
                   "CALU__AMG510_vs_NAIVE", "CALU__MRTX_vs_NAIVE")
  contrast_labels <- c("H23\nsoto", "H23\nada",
                       "CALU\nsoto", "CALU\nada")

  lfc_mat <- matrix(NA_real_, nrow = length(ras_ensg), ncol = length(contrasts_v),
                    dimnames = list(names(ras_ensg), contrast_labels))
  sig_mat <- matrix(FALSE, nrow = length(ras_ensg), ncol = length(contrasts_v),
                    dimnames = list(names(ras_ensg), contrast_labels))
  for (j in seq_along(contrasts_v)) {
    sub <- deg %>%
      filter(contrast == contrasts_v[j], gene_id_nv %in% ras_ensg)
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
    annotation_name_gp = gpar(fontsize = 8),
    annotation_name_side = "left",
    height = unit(7, "mm"),
    show_legend = TRUE
  )

  # Left annotation: 4 anno_simple tracks for log2FC, with "*" overlay when sig
  star_pt <- function(sig_col) ifelse(sig_col, "*", "")
  left_anno <- rowAnnotation(
    "H23 soto"  = anno_simple(lfc_mat[, 1], col = LFC_RAMP_RAS,
                              pch = star_pt(sig_mat[, 1]),
                              pt_size = unit(3, "mm"),
                              gp = gpar(col = "grey60", lwd = 0.4)),
    "H23 ada"   = anno_simple(lfc_mat[, 2], col = LFC_RAMP_RAS,
                              pch = star_pt(sig_mat[, 2]),
                              pt_size = unit(3, "mm"),
                              gp = gpar(col = "grey60", lwd = 0.4)),
    "CALU soto" = anno_simple(lfc_mat[, 3], col = LFC_RAMP_RAS,
                              pch = star_pt(sig_mat[, 3]),
                              pt_size = unit(3, "mm"),
                              gp = gpar(col = "grey60", lwd = 0.4)),
    "CALU ada"  = anno_simple(lfc_mat[, 4], col = LFC_RAMP_RAS,
                              pch = star_pt(sig_mat[, 4]),
                              pt_size = unit(3, "mm"),
                              gp = gpar(col = "grey60", lwd = 0.4)),
    annotation_name_gp   = gpar(fontsize = 7),
    annotation_name_rot  = 60,
    annotation_name_side = "bottom",
    simple_anno_size     = unit(4, "mm"),
    gap                  = unit(0.6, "mm")
  )

  body_h <- unit(nrow(m_scaled) * 0.5, "cm")  # 4 rows × 0.5 cm = 2 cm body
  body_w <- unit(0.3 * ncol(m_scaled), "cm")  # ~5.4 cm body width

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
    column_title         = "RAS family — H23 / CALU expression and drug log2FC",
    column_title_gp      = gpar(fontsize = 10, fontface = "bold"),
    top_annotation       = top_anno,
    right_annotation     = left_anno,
    heatmap_legend_param = list(
      title          = "row z-score",
      title_gp       = gpar(fontsize = 8),
      labels_gp      = gpar(fontsize = 8),
      title_position = "topcenter")
  )

  # Stand-alone log2FC legend (asterisk = padj<0.01 & |log2FC|>0.5)
  lfc_legend <- Legend(
    col_fun  = LFC_RAMP_RAS,
    title    = "log2FC vs naive",
    at       = c(-3, 0, 3),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8))
  sig_legend <- Legend(
    labels   = "padj<0.01 & |log2FC|>0.5",
    title    = "significance",
    type     = "points",
    pch      = "*",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    background = "white")

  list(hm = hm,
       extra_legends = list(lfc_legend, sig_legend))
}

# ---- unified heatmap: both cell lines × shared genes -----------------------
# One Heatmap with 18 columns (9 H23 + 9 CALU). Rows = intersection of the
# two top-N selections. Z-score per row across all 18 samples so values are
# directly comparable across the two cell-line blocks. Column-split on
# cell_line so the blocks are visually separated.
build_unified_shared_heatmap <- function(sig_ensg, sig_label, n,
                                         gene_name_lookup, shared_ensg) {
  if (length(shared_ensg) < 2) return(NULL)

  smpl <- bind_rows(samples_for_line(counts_df, "H23"),
                    samples_for_line(counts_df, "CALU")) %>%
    mutate(cell_line = factor(cell_line, levels = c("H23","CALU"))) %>%
    arrange(cell_line, treatment, sample)

  raw <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  norm <- normalise_counts(raw)

  keep <- intersect(rownames(norm), shared_ensg)
  if (length(keep) < 2) return(NULL)
  m <- norm[keep, smpl$sample, drop = FALSE]

  sym <- gene_name_lookup[rownames(m)]
  sym[is.na(sym) | !nzchar(sym)] <- rownames(m)[is.na(sym) | !nzchar(sym)]
  dup <- duplicated(sym)
  if (any(dup)) sym[dup] <- paste0(sym[dup], "_", seq_len(sum(dup)))
  rownames(m) <- sym

  m_scaled <- t(scale(t(m)))   # row z-score across all 18 samples

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

  # Right CNV strip: 6 columns (H23 × {naive, soto, ada} + CALU × {naive, soto, ada}).
  cnv_h23  <- cnv_cnt_matrix("H23",  rownames(m_scaled))
  cnv_calu <- cnv_cnt_matrix("CALU", rownames(m_scaled))
  right_anno <- NULL
  if (!is.null(cnv_h23) && !is.null(cnv_calu)) {
    cnv_combo <- cbind(cnv_h23, cnv_calu)
    colnames(cnv_combo) <- c(paste0("H23_",  colnames(cnv_h23)),
                             paste0("CALU_", colnames(cnv_calu)))
    args <- c(
      lapply(seq_len(ncol(cnv_combo)),
             function(j) cnv_combo[, j]),
      list(col = stats::setNames(
              rep(list(CNV_CNT_RAMP), ncol(cnv_combo)),
              colnames(cnv_combo)),
           na_col = "white",
           show_legend = c(TRUE, rep(FALSE, ncol(cnv_combo) - 1)),
           annotation_legend_param = stats::setNames(
              list(list(title    = "CNt",
                        title_gp = gpar(fontsize = 8),
                        labels_gp = gpar(fontsize = 8),
                        at       = c(0, 2, 4, 6))),
              colnames(cnv_combo)[1]),
           annotation_name_gp   = gpar(fontsize = 6),
           annotation_name_rot  = 90,
           annotation_name_side = "top",
           simple_anno_size     = unit(2.5, "mm"),
           gap                  = unit(0.4, "mm")))
    names(args)[seq_len(ncol(cnv_combo))] <- colnames(cnv_combo)
    right_anno <- do.call(rowAnnotation, args)
  }

  body_h <- unit(nrow(m_scaled) * PER_ROW_CM, "cm")
  body_w <- unit(0.3 * ncol(m_scaled), "cm")     # ~0.3 cm per sample column
  Heatmap(
    m_scaled,
    name        = "z-score",
    col         = HEATMAP_RAMP,
    height      = body_h,
    width       = body_w,
    cluster_rows         = TRUE,
    cluster_columns      = FALSE,
    show_row_dend        = TRUE,
    show_column_dend     = FALSE,
    show_row_names       = TRUE,
    row_names_side       = "right",
    row_names_gp         = gpar(fontsize = ROW_LABEL_PT),
    show_column_names    = FALSE,
    column_split         = smpl$cell_line,
    column_title         = sprintf("%s — top%d ↑+↓ ∩ across H23 / CALU (%d genes)",
                                    sig_label, n, nrow(m_scaled)),
    column_title_gp      = gpar(fontsize = 9, fontface = "bold"),
    top_annotation       = top_anno,
    right_annotation     = if (!is.null(right_anno)) right_anno else NULL,
    heatmap_legend_param = list(
      title          = "row z-score",
      title_gp       = gpar(fontsize = 8),
      labels_gp      = gpar(fontsize = 8),
      title_position = "topcenter")
  )
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
    annotation_name_gp = gpar(fontsize = 8),
    annotation_name_side = "left",
    height = unit(4, "mm"),
    show_legend = TRUE
  )

  # Left annotation: flags genes that also appear in the OTHER cell line's
  # top-N heatmap for the same signature (= cross-line "shared" picks).
  # Suppressed when shared_only = TRUE (every row is shared by construction).
  left_anno <- NULL
  if (!shared_only && !is.null(shared_syms) && length(shared_syms) > 0) {
    in_both <- ifelse(rownames(m_scaled) %in% shared_syms,
                      "shared", "unique")
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

  # Right annotation: per-gene CNt at each of the 3 conditions, mapped via
  # the continuous CNV_CNT_RAMP. Built from the CNV table for this cell
  # line (NULL for lines without CNV — e.g. mouse MEM).
  cnv_m <- cnv_cnt_matrix(line, rownames(m_scaled))
  right_anno <- NULL
  if (!is.null(cnv_m)) {
    right_anno <- rowAnnotation(
      naive     = cnv_m[, 1],
      sotorasib = cnv_m[, 2],
      adagrasib = cnv_m[, 3],
      col = list(naive     = CNV_CNT_RAMP,
                 sotorasib = CNV_CNT_RAMP,
                 adagrasib = CNV_CNT_RAMP),
      na_col = "white",
      show_legend = c(naive = TRUE, sotorasib = FALSE, adagrasib = FALSE),
      annotation_legend_param = list(
        naive = list(title       = "CNt",
                     title_gp    = gpar(fontsize = 8),
                     labels_gp   = gpar(fontsize = 8),
                     at          = c(0, 2, 4, 6))),
      annotation_name_gp   = gpar(fontsize = 7),
      annotation_name_rot  = 90,
      annotation_name_side = "top",
      simple_anno_size     = unit(3.5, "mm"),
      gap                  = unit(0.5, "mm")
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
    left_annotation      = if (!is.null(left_anno))  left_anno  else NULL,
    right_annotation     = if (!is.null(right_anno)) right_anno else NULL,
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
           cell_line = factor(cell_line, levels = c("H23","CALU")),
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

write_pdf_grob <- function(grob, path, w, h) {
  grDevices::pdf(path, width = w, height = h, useDingbats = FALSE, bg = "white")
  grid::grid.newpage()
  grid::grid.draw(grob)
  invisible(dev.off())
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

    # Cross-line shared genes (intersection of the two top-N selections,
    # mapped to symbols for matching against heatmap row names).
    picks_h23  <- top_per_line("H23",  sig$ensg, n)
    picks_calu <- top_per_line("CALU", sig$ensg, n)
    shared_ensg <- intersect(union(picks_h23$up,  picks_h23$down),
                             union(picks_calu$up, picks_calu$down))
    shared_syms <- unname(gene_name_lookup[shared_ensg])
    shared_syms <- shared_syms[!is.na(shared_syms) & nzchar(shared_syms)]
    message(sprintf("    shared (in both top-%d heatmaps): %d genes",
                    n, length(shared_syms)))

    hm_h23   <- build_heatmap_one_line("H23",  sig$ensg, sig$label, gene_name_lookup, n,
                                        shared_syms = shared_syms)
    hm_calu  <- build_heatmap_one_line("CALU", sig$ensg, sig$label, gene_name_lookup, n,
                                        shared_syms = shared_syms)
    p_bar    <- build_count_bar(sig$ensg, sig$label)

    # Per-cell-line heatmap PDFs. Body is fixed at nrow × PER_ROW_CM so the
    # cm spec is honoured; total page = body + ~2 in for title / annotation /
    # legend / margins.
    body_to_in <- function(n_rows) (n_rows * PER_ROW_CM) / 2.54
    for (pair in list(list(hm = hm_h23,  line = "H23"),
                      list(hm = hm_calu, line = "CALU"))) {
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

    # Bar plot — duplicated per threshold subfolder so each top-N folder is
    # self-contained. The DEG-count summary itself doesn't depend on N.
    bar_out <- file.path(thr_dir, sprintf("barplot_%s_DEG_counts.pdf", sig$tag))
    ggsave(bar_out, p_bar, width = 5, height = 2.6, useDingbats = FALSE)
    message(sprintf("Wrote %s", bar_out))

    # Combined panel: 2 heatmaps side-by-side on top, bar below
    hm_h23_grob  <- grid.grabExpr(draw(hm_h23,  heatmap_legend_side = "right",
                                       annotation_legend_side = "right",
                                       merge_legend = TRUE))
    hm_calu_grob <- grid.grabExpr(draw(hm_calu, heatmap_legend_side = "right",
                                       annotation_legend_side = "right",
                                       merge_legend = TRUE))

    p_bar_grob <- grid.grabExpr(print(p_bar))
    panel <- (patchwork::wrap_elements(full = hm_h23_grob) |
              patchwork::wrap_elements(full = hm_calu_grob)) /
             patchwork::wrap_elements(full = p_bar_grob) +
             patchwork::plot_layout(heights = c(3, 1)) +
             patchwork::plot_annotation(
               title = sprintf("%s pathway under sotorasib / adagrasib (top%d ↑ + ↓ per cell line)",
                               sig$label, n),
               theme = theme(plot.title = element_text(size = 11, face = "bold")))

    panel_out <- file.path(thr_dir, sprintf("panel_%s_summary.pdf", sig$tag))
    panel_w <- 9
    max_body_in <- (max(nrow(hm_h23@matrix), nrow(hm_calu@matrix)) *
                    PER_ROW_CM) / 2.54
    panel_h <- max_body_in + 4.6   # heatmap body + bar plot + title margins
    ggsave(panel_out, panel, width = panel_w, height = panel_h,
           useDingbats = FALSE, limitsize = FALSE)
    message(sprintf("Wrote %s   (%.1f × %.1f in)", panel_out, panel_w, panel_h))

    # ------------------------------------------------------------------
    # Shared-only variant: same heatmaps + bar but rows restricted to the
    # cross-line intersection. Suppresses the "in both" left strip since
    # every row is shared. Skipped if the intersection is too small.
    # ------------------------------------------------------------------
    if (length(shared_ensg) >= 2) {
      hm_h23_sh  <- build_heatmap_one_line(
        "H23",  sig$ensg, sig$label, gene_name_lookup, n,
        shared_only = TRUE, shared_ensg = shared_ensg)
      hm_calu_sh <- build_heatmap_one_line(
        "CALU", sig$ensg, sig$label, gene_name_lookup, n,
        shared_only = TRUE, shared_ensg = shared_ensg)

      for (pair in list(list(hm = hm_h23_sh,  line = "H23"),
                        list(hm = hm_calu_sh, line = "CALU"))) {
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

      if (!is.null(hm_h23_sh) && !is.null(hm_calu_sh)) {
        hm_h23_sh_grob  <- grid.grabExpr(draw(hm_h23_sh,  heatmap_legend_side = "right",
                                              annotation_legend_side = "right",
                                              merge_legend = TRUE))
        hm_calu_sh_grob <- grid.grabExpr(draw(hm_calu_sh, heatmap_legend_side = "right",
                                              annotation_legend_side = "right",
                                              merge_legend = TRUE))
        p_bar_grob_sh <- grid.grabExpr(print(p_bar))
        panel_sh <- (patchwork::wrap_elements(full = hm_h23_sh_grob) |
                     patchwork::wrap_elements(full = hm_calu_sh_grob)) /
                    patchwork::wrap_elements(full = p_bar_grob_sh) +
                    patchwork::plot_layout(heights = c(3, 1)) +
                    patchwork::plot_annotation(
                      title = sprintf("%s — shared genes only (top%d ↑+↓ ∩ across H23 / CALU; %d genes)",
                                      sig$label, n, length(shared_ensg)),
                      theme = theme(plot.title = element_text(size = 11, face = "bold")))
        panel_sh_out <- file.path(thr_dir,
                                  sprintf("panel_%s_summary_shared.pdf", sig$tag))
        max_body_sh_in <- (max(nrow(hm_h23_sh@matrix),
                               nrow(hm_calu_sh@matrix)) * PER_ROW_CM) / 2.54
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

    # ------------------------------------------------------------------
    # Unified compact heatmap: ONE Heatmap with both cell lines (18 cols),
    # rows = the cross-line intersection. Z-scored row-wise across all
    # 18 samples so the two cell-line blocks are directly comparable.
    # ------------------------------------------------------------------
    if (length(shared_ensg) >= 2) {
      hm_unified <- build_unified_shared_heatmap(
        sig$ensg, sig$label, n,
        gene_name_lookup = gene_name_lookup,
        shared_ensg = shared_ensg)
      if (!is.null(hm_unified)) {
        n_rows <- nrow(hm_unified@matrix)
        h_unified <- (n_rows * PER_ROW_CM) / 2.54 + 2.4
        out_unified <- file.path(thr_dir,
                                  sprintf("heatmap_%s_unified_shared.pdf", sig$tag))
        grDevices::pdf(out_unified, width = 8, height = h_unified,
                       useDingbats = FALSE, bg = "white")
        draw(hm_unified, heatmap_legend_side = "right",
             annotation_legend_side = "right", merge_legend = TRUE)
        invisible(dev.off())
        message(sprintf("Wrote %s   (8 × %.1f in)", out_unified, h_unified))
      }
    }
  }
}

# Standalone RAS-family heatmap intentionally skipped in the _min variant.
