#!/usr/bin/env Rscript
# =============================================================================
# heatmaps_MAPK_ERK.R
# Per cell line × per gene list, Heatmap3-style heatmap of DEG genes that fall
# in the signature. Input is the raw count matrix; values are CPM-normalised
# log2(counts+1) and z-scored per row inside get_heatmap3().
#
# Two gene lists, two cell lines  →  4 PDFs:
#   Figures/RNA/heatmap_RAS84_H23.pdf
#   Figures/RNA/heatmap_RAS84_CALU.pdf
#   Figures/RNA/heatmap_MAPK_ERK_H23.pdf
#   Figures/RNA/heatmap_MAPK_ERK_CALU.pdf
#
# Sample → treatment mapping (verified by Spearman ≈ ±1 between inferred
# log2FC and the deg log2FC for both contrasts and both cell lines):
#   second digit of trailing condition code = treatment
#     0 → NAIVE   1 → sotorasib (AMG510)   2 → adagrasib (MRTX)
# Example: caks-h2301-50-134-112  → 2nd digit = 1 → sotorasib.
#
# Inputs : input/H23_CALU_counts.rds   (15226 × 22; 4 metadata + 18 samples)
#          input/H23_CALU_deg.rds
#          input/RAS_84_PhilipEast_et_al_NatCom_2022.csv
#          input/mapk_erk_signature_East_NatCom_2022.rds
#          Script/utils_RNA.R          (provides get_heatmap3)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(scico)
})

# scico::vik 3-anchor sample (passed to get_heatmap3 which builds a 3-stop
# colorRamp2 between c(-2, 0, 2)). Linear interpolation between these three
# stops approximates the perceptually-uniform vik gradient.
VIK_3 <- scico::scico(3, palette = "vik")

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
INPUT_DIR <- file.path(.project_root, "input")
OUT_DIR   <- file.path(.project_root, "Figures", "RNA_min")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Sweep stale outputs from the previous (combined-signature) version
old <- list.files(OUT_DIR, pattern = "heatmap_RAS84_MAPK_ERK_.*\\.pdf$", full.names = TRUE)
if (length(old)) file.remove(old)

source(file.path(.script_dir, "..", "utils", "utils_RNA.R"))   # get_heatmap3 + ComplexHeatmap deps

PADJ_MAX <- 0.01
LFC_MIN  <- 0.5

TREATMENT_LEVELS <- c("NAIVE", "sotorasib", "adagrasib")
TREATMENT_COL    <- c(NAIVE     = "grey60",
                      sotorasib = "#0072B2",
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
counts_df <- readRDS(file.path(INPUT_DIR, "H23_CALU_counts.rds"))
deg       <- readRDS(file.path(INPUT_DIR, "H23_CALU_deg.rds"))

mapk_erk <- readRDS(file.path(INPUT_DIR, GENE_LISTS$MAPK_ERK$file))

mapk_ensg <- unique(mapk_erk$human_ensembl_gene)
mapk_ensg <- mapk_ensg[!is.na(mapk_ensg) & nzchar(mapk_ensg)]

SIG_GENES <- list(MAPK_ERK = mapk_ensg)

# ---- helpers ----------------------------------------------------------------

# Strip the version suffix (".13", ".11", ...) from rownames so we match
# deg$gene_id_nv and signature ENSGs.
strip_version <- function(x) sub("\\..*$", "", x)

samples_for_line <- function(line) {
  tag <- if (line == "H23") "h2301" else if (line == "CALU") "cal01" else stop(line)
  cols <- grep(tag, colnames(counts_df), value = TRUE)
  code <- sub(".*-", "", cols)             # last segment (e.g. "112")
  treat <- unname(DIGIT_TO_TREAT[substr(code, 2, 2)])
  data.frame(sample = cols, treatment = treat,
             stringsAsFactors = FALSE) %>%
    mutate(treatment = factor(treatment, levels = TREATMENT_LEVELS)) %>%
    arrange(treatment, sample)
}

normalise_counts <- function(mat) {
  # Library-size normalisation → CPM, then log2(CPM + 1) for visualisation
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

# Per cell line × MAPK/ERK signature, return up and down DEG sets defined as
# "top-N in sotorasib OR top-N in adagrasib" — i.e. for each contrast take
# the top-N by log2FC and union the two sets per direction. Gene must clear
# the DEG filter (padj < PADJ_MAX, |log2FC| > LFC_MIN) in the contrast it
# was selected from.
top_up_down_for_line <- function(line, sig_ensg, n = 50) {
  contrasts_keep <- grep(paste0("^", line, "__"),
                         unique(deg$contrast), value = TRUE)
  d <- deg %>%
    filter(contrast %in% contrasts_keep,
           gene_id_nv %in% sig_ensg,
           !is.na(padj),
           !is.na(log2FoldChange),
           padj < PADJ_MAX,
           abs(log2FoldChange) > LFC_MIN)

  pick <- function(per_contrast_filter, slicer) {
    d %>% filter(per_contrast_filter(log2FoldChange)) %>%
      group_by(contrast) %>%
      slicer(log2FoldChange, n = n) %>%
      ungroup() %>%
      pull(gene_id_nv) %>% unique()
  }
  up   <- pick(function(x) x > 0,  slice_max)
  down <- pick(function(x) x < 0,  slice_min)
  list(up = up, down = down)
}

# ---- driver -----------------------------------------------------------------
for (line in c("H23", "CALU")) {

  smpl <- samples_for_line(line)
  raw  <- as.matrix(counts_df[, smpl$sample])
  rownames(raw) <- strip_version(counts_df$Geneid)
  gene_name_lookup <- setNames(counts_df$gene_name, rownames(raw))

  cl_deg <- cell_line_deg_genes(line)
  norm_full <- normalise_counts(raw)        # genes × samples, log2(CPM+1)

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

    # Replace ENSG rownames with gene symbols (fall back to ENSG if missing)
    sym <- gene_name_lookup[rownames(m)]
    sym[is.na(sym) | !nzchar(sym)] <- rownames(m)[is.na(sym) | !nzchar(sym)]
    # If duplicated symbols slip through, append an index to keep names unique
    dup <- duplicated(sym)
    if (any(dup)) sym[dup] <- paste0(sym[dup], "_", seq_len(sum(dup)))
    rownames(m) <- sym

    # Column order: NAIVE first, then sotorasib, then adagrasib
    m <- m[, smpl$sample, drop = FALSE]

    # Annotation for the column-side strip
    annotDF <- data.frame(treatment = smpl$treatment,
                          row.names = smpl$sample)
    annotCol <- list(treatment = TREATMENT_COL)

    # Show row labels only when the row count is small enough to read at 8 pt
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

# ---- additional heatmap: top-50 up + top-50 down MAPK/ERK DEGs --------------
# Per cell line, restrict to MAPK/ERK signature DEGs and pick the top 50 most
# upregulated and 50 most downregulated by max / min log2FC across the two
# drug contrasts. Row labels stay on (8 pt, get_heatmap3 default) — typical
# row count is ≤ 100 so labels are readable.
for (TOPN in c(20))    # minimal: top-20 only
for (line in c("H23", "CALU")) {
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

  # Order rows: most upregulated first, then most downregulated. Within each
  # block keep the slice_min/slice_max ordering.
  ordered <- c(intersect(picks$up,  keep_ensg),
               setdiff(picks$down, picks$up))
  m <- norm_full[ordered, smpl$sample, drop = FALSE]

  sym <- gene_name_lookup[rownames(m)]
  sym[is.na(sym) | !nzchar(sym)] <- rownames(m)[is.na(sym) | !nzchar(sym)]
  dup <- duplicated(sym)
  if (any(dup)) sym[dup] <- paste0(sym[dup], "_", seq_len(sum(dup)))
  rownames(m) <- sym

  direction <- ifelse(rownames(m) %in% gene_name_lookup[picks$up], "up", "down")
  # The lookup above can mismap if symbols collide; safer: re-derive from order
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
