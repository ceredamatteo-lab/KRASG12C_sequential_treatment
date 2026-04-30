#!/usr/bin/env Rscript
# =============================================================================
# heatmaps_MAPK_ERK_MEM_MEP.R
# Mouse counterpart of heatmaps_MAPK_ERK.R, for the MEF lines MEM and MEP.
#
# Per cell line × per gene list, Heatmap3-style heatmap of DEG genes that
# fall in the (mouse-mapped) RAS84 or MAPK/ERK signature. Input is the raw
# count matrix; values are CPM-normalised log2(counts+1) and z-scored per
# row inside get_heatmap3(). scico::vik palette (3-anchor approximation).
#
# Two gene lists, two cell lines  →  4 PDFs:
#   Figures/RNA/MEM_MEP/heatmap_RAS84_MEM.pdf
#   Figures/RNA/MEM_MEP/heatmap_RAS84_MEP.pdf
#   Figures/RNA/MEM_MEP/heatmap_MAPK_ERK_MEM.pdf
#   Figures/RNA/MEM_MEP/heatmap_MAPK_ERK_MEP.pdf
# Plus the additional top-50 / top-25 ↑+↓ MAPK/ERK heatmaps per cell line.
#
# Sample → treatment mapping (verified by Spearman ≈ +1 between inferred
# log2FC and the deg log2FC for both contrasts and both cell lines):
#   second digit of trailing condition code
#     0 → NAIVE   1 → sotorasib (AMG510)   2 → adagrasib (MRTX)
#
# Mouse gene-symbol mapping:
#   RAS84    : ras84$CCLE_feature_id (already mouse symbols)
#   MAPK/ERK : human symbols → mouse via title-case (ABCA7 → Abca7).
#              Matches >95% of orthologues; non-trivial cases like
#              TP53→Trp53 are missed and logged.
#
# Inputs : input/MEM_MEP_counts.rds  (21701 × 23; 4 metadata + 18 samples + gene_id_nv col)
#          input/MEM_MEP_deg.rds
#          input/RAS_84_PhilipEast_et_al_NatCom_2022.csv
#          input/mapk_erk_signature_East_NatCom_2022.rds
#          Script/utils_RNA.R          (provides get_heatmap3)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(scico)
})

VIK_3 <- scico::scico(3, palette = "vik")   # 3-anchor sample for get_heatmap3

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
INPUT_DIR <- file.path(.project_root, "input")
OUT_DIR   <- file.path(.project_root, "Figures", "RNA_min", "MEM_MEP")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Sweep stale outputs
old <- list.files(OUT_DIR, pattern = "^heatmap_.*\\.pdf$",
                  full.names = TRUE)
if (length(old)) file.remove(old)

source(file.path(.script_dir, "..", "utils", "utils_RNA.R"))

PADJ_MAX <- 0.01
LFC_MIN  <- 0.5

DRUG_LABEL <- c(AMG510 = "sotorasib", MRTX = "adagrasib")
DRUG_ORDER <- c("sotorasib", "adagrasib")
DRUG_COL   <- c(sotorasib = "#0072B2", adagrasib = "#D55E00")

TREATMENT_LEVELS <- c("NAIVE","sotorasib","adagrasib")
TREATMENT_COL    <- c(NAIVE = "grey60", sotorasib = "#0072B2",
                      adagrasib = "#D55E00")
DIGIT_TO_TREAT   <- c("0" = "NAIVE", "1" = "sotorasib", "2" = "adagrasib")

# Minimal variant: MAPK/ERK only — RAS84 dropped per project request.
GENE_LISTS <- list(
  MAPK_ERK = list(
    label = "MAPK/ERK (GOBP_MAPK_CASCADE ∪ GOBP_ERK1_AND_ERK2_CASCADE; East 2022)",
    file  = "mapk_erk_signature_East_NatCom_2022.rds"
  )
)

# ---- inputs -----------------------------------------------------------------
counts_df <- readRDS(file.path(INPUT_DIR, "MEM_MEP_counts.rds"))
deg       <- readRDS(file.path(INPUT_DIR, "MEM_MEP_deg.rds"))

mapk_erk <- readRDS(file.path(INPUT_DIR, GENE_LISTS$MAPK_ERK$file))

# ---- gene-symbol → mouse ENSMUSG mapping -----------------------------------
strip_version <- function(x) sub("\\..*$", "", x)

# Build a name → ENSMUSG lookup from the counts table (which has mouse symbols).
# When a symbol maps to multiple ENSMUSGs, keep the first.
mouse_sym_to_ensg <- counts_df %>%
  transmute(gene_id_nv = strip_version(Geneid), gene_name) %>%
  filter(!is.na(gene_name), nzchar(gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE)
sym2ensg <- setNames(mouse_sym_to_ensg$gene_id_nv, mouse_sym_to_ensg$gene_name)

# Title-case helper for mouse orthologue convention.
to_mouse <- function(s) {
  s <- as.character(s)
  s <- s[!is.na(s) & nzchar(s)]
  paste0(substr(s, 1, 1),
         tolower(substr(s, 2, nchar(s))))
}

# MAPK/ERK: human symbols → mouse via title-case → ENSMUSG.
mapk_erk_sym_human <- unique(mapk_erk$gene_symbol)
mapk_erk_sym_mouse <- to_mouse(mapk_erk_sym_human)
mapk_erk_ensg <- unname(sym2ensg[mapk_erk_sym_mouse])
mapk_erk_ensg <- mapk_erk_ensg[!is.na(mapk_erk_ensg)]
mapk_erk_missing <- setdiff(mapk_erk_sym_mouse, names(sym2ensg))

message(sprintf(
  "MAPK_ERK: %d input symbols  -> %d mapped to ENSMUSG  (%d unmatched)",
  length(mapk_erk_sym_human), length(mapk_erk_ensg), length(mapk_erk_missing)))

SIG_GENES <- list(MAPK_ERK = mapk_erk_ensg)

# ---- helpers ----------------------------------------------------------------

samples_for_line <- function(line) {
  tag <- if (line == "MEM") "mem01" else if (line == "MEP") "mep01" else stop(line)
  cols <- grep(tag, colnames(counts_df), value = TRUE)
  code <- sub(".*-", "", cols)
  treat <- unname(DIGIT_TO_TREAT[substr(code, 2, 2)])
  data.frame(sample = cols, treatment = treat,
             stringsAsFactors = FALSE) %>%
    mutate(treatment = factor(treatment, levels = TREATMENT_LEVELS)) %>%
    arrange(treatment, sample)
}

normalise_counts <- function(mat) {
  lib <- colSums(mat)
  cpm <- t(t(mat) / lib) * 1e6
  log2(cpm + 1)
}

cell_line_deg_genes <- function(line) {
  contrasts_keep <- grep(paste0("^", line, "__"),
                         unique(deg$contrast), value = TRUE)
  deg %>%
    filter(contrast %in% contrasts_keep,
           !is.na(padj),
           padj < PADJ_MAX,
           abs(log2FoldChange) > LFC_MIN) %>%
    pull(gene_id_nv) %>% unique()
}

top_up_down_for_line <- function(line, sig_ensg, n = 50) {
  contrasts_keep <- grep(paste0("^", line, "__"),
                         unique(deg$contrast), value = TRUE)
  d <- deg %>%
    filter(contrast %in% contrasts_keep,
           gene_id_nv %in% sig_ensg,
           !is.na(padj),
           !is.na(log2FoldChange))
  deg_genes <- d %>%
    filter(padj < PADJ_MAX, abs(log2FoldChange) > LFC_MIN) %>%
    pull(gene_id_nv) %>% unique()
  d <- d %>% filter(gene_id_nv %in% deg_genes)

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

# ---- driver: per (cell line × gene list) heatmaps --------------------------
for (line in c("MEM", "MEP")) {

  smpl <- samples_for_line(line)
  raw  <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  gene_name_lookup <- setNames(counts_df$gene_name, rownames(raw))

  cl_deg <- cell_line_deg_genes(line)
  norm_full <- normalise_counts(raw)

  for (gname in names(GENE_LISTS)) {
    keep_ensg <- intersect(intersect(rownames(norm_full), cl_deg),
                           SIG_GENES[[gname]])
    message(sprintf("[%s / %s] sig DEGs ∩ signature : %d genes",
                    line, gname, length(keep_ensg)))

    if (length(keep_ensg) < 2) {
      message(sprintf("  too few genes — skipping plot"))
      next
    }

    m <- norm_full[keep_ensg, smpl$sample, drop = FALSE]
    sym <- gene_name_lookup[rownames(m)]
    sym[is.na(sym) | !nzchar(sym)] <- rownames(m)[is.na(sym) | !nzchar(sym)]
    dup <- duplicated(sym)
    if (any(dup)) sym[dup] <- paste0(sym[dup], "_", seq_len(sum(dup)))
    rownames(m) <- sym

    annotDF  <- data.frame(treatment = smpl$treatment, row.names = smpl$sample)
    annotCol <- list(treatment = TREATMENT_COL)

    show_rn <- nrow(m) <= 80

    out_pdf <- file.path(OUT_DIR,
                         sprintf("heatmap_%s_%s.pdf", gname, line))
    page_h <- max(5, min(22, 0.18 * nrow(m) + 3))
    page_w <- if (show_rn) 7 else 5

    grDevices::pdf(out_pdf, width = page_w, height = page_h,
                   useDingbats = FALSE, bg = "white")
    get_heatmap3(m,
                 annotDF       = annotDF,
                 annotCol      = annotCol,
                 show_rownames = show_rn,
                 myPalette     = VIK_3,
                 cluster_columns = FALSE,
                 column_title  = sprintf("%s — %s  (n=%d)",
                                         line, GENE_LISTS[[gname]]$label,
                                         nrow(m)))
    invisible(dev.off())
    message(sprintf("  Wrote %s   (%.1f × %.1f in)", out_pdf, page_w, page_h))
  }
}

# ---- additional: top-50 / top-25 up + down MAPK/ERK heatmaps ---------------
for (TOPN in c(20))    # minimal: top-20 only
for (line in c("MEM", "MEP")) {
  smpl <- samples_for_line(line)
  raw  <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  gene_name_lookup <- setNames(counts_df$gene_name, rownames(raw))
  norm_full <- normalise_counts(raw)

  picks <- top_up_down_for_line(line, SIG_GENES$MAPK_ERK, n = TOPN)
  keep_ensg <- intersect(rownames(norm_full),
                         union(picks$up, picks$down))
  if (length(keep_ensg) < 2) {
    message(sprintf("[%s / MAPK_ERK top%d±]  too few genes (%d) — skipping",
                    line, TOPN, length(keep_ensg)))
    next
  }
  message(sprintf("[%s / MAPK_ERK top%d±]  up=%d  down=%d  union=%d  (overlap=%d)",
                  line, TOPN,
                  length(picks$up), length(picks$down),
                  length(union(picks$up, picks$down)),
                  length(intersect(picks$up, picks$down))))

  ordered <- c(intersect(picks$up,  keep_ensg),
               setdiff(picks$down, picks$up))
  m <- norm_full[ordered, smpl$sample, drop = FALSE]

  sym <- gene_name_lookup[rownames(m)]
  sym[is.na(sym) | !nzchar(sym)] <- rownames(m)[is.na(sym) | !nzchar(sym)]
  dup <- duplicated(sym)
  if (any(dup)) sym[dup] <- paste0(sym[dup], "_", seq_len(sum(dup)))
  rownames(m) <- sym

  direction <- c(rep("up",   length(intersect(picks$up,  keep_ensg))),
                 rep("down", length(setdiff(picks$down, picks$up))))

  annotDF  <- data.frame(treatment = smpl$treatment, row.names = smpl$sample)
  annotCol <- list(treatment = TREATMENT_COL)

  row_anno <- ComplexHeatmap::rowAnnotation(
    direction = direction,
    col = list(direction = c(up = "#CC0000", down = "#3333FF")),
    annotation_legend_param = list(
      title_gp  = grid::gpar(fontsize = 8),
      labels_gp = grid::gpar(fontsize = 8)),
    annotation_name_gp = grid::gpar(fontsize = 8),
    width = grid::unit(3, "mm")
  )

  out_pdf <- file.path(OUT_DIR,
                       sprintf("heatmap_MAPK_ERK_top%d_updwn_%s.pdf", TOPN, line))
  page_h <- max(8, min(24, 0.18 * nrow(m) + 3))
  page_w <- 7

  grDevices::pdf(out_pdf, width = page_w, height = page_h,
                 useDingbats = FALSE, bg = "white")
  get_heatmap3(m,
               annotDF       = annotDF,
               annotCol      = annotCol,
               rowAnnot      = row_anno,
               show_rownames = TRUE,
               myPalette     = VIK_3,
               cluster_rows  = TRUE,
               cluster_columns = FALSE,
               column_title  = sprintf("%s — MAPK/ERK top%d ↑ + top%d ↓ DEGs (n=%d)",
                                       line, TOPN, TOPN, nrow(m)))
  invisible(dev.off())
  message(sprintf("  Wrote %s   (%.1f × %.1f in)", out_pdf, page_w, page_h))
}
