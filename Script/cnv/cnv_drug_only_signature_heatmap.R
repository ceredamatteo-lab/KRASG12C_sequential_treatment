#!/usr/bin/env Rscript
# =============================================================================
# cnv_drug_only_signature_heatmap.R
# Per cell line × per signature, two heatmaps of signature genes that gain
# or lose copies ONLY under drug treatment (Amp/Loss in sotorasib or
# adagrasib, absent in naive):
#
#   filter "drug_only"  : any Amp/Loss in 010 or 020, none in 000
#   filter "extreme"    : drug_only AND the drug-condition CNt is < 1
#                         (strong loss) or > 3 (strong gain)
#
# Outputs (in Figures/CNV/drug_only/):
#   heatmap_<sig>_<line>.pdf            — basic drug-only filter
#   heatmap_<sig>_<line>_extreme.pdf    — strict CNt<1 OR CNt>3 filter
#
# 2 signatures × 3 lines × 2 filters = 12 PDFs. Body height capped at 10 cm.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
INPUT_DIR <- file.path(.project_root, "input")
OUT_DIR   <- file.path(.project_root, "Figures", "CNV", "drug_only")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Sweep stale outputs from earlier layouts
old <- list.files(OUT_DIR, pattern = "\\.pdf$", full.names = TRUE)
if (length(old)) file.remove(old)
old_flat <- list.files(file.path(.project_root, "Figures", "CNV"),
                       pattern = "^heatmap_drug_only_.*\\.pdf$",
                       full.names = TRUE)
if (length(old_flat)) file.remove(old_flat)

CONDS      <- c("000", "010", "020")
COND_LABEL <- c("000" = "naive", "010" = "sotorasib", "020" = "adagrasib")
COND_COL   <- c(naive = "grey60", sotorasib = "#0072B2", adagrasib = "#D55E00")

CNT_COL_LOW  <- "#2166AC"
CNT_COL_MID  <- "#F7F7F7"
CNT_COL_HIGH <- "#B2182B"
CNT_COL_NA   <- "grey92"

EXTREME_LO <- 1   # CNt < 1  → strong loss
EXTREME_HI <- 3   # CNt > 3  → strong gain
# Horizontal layout: 3 rows (naive/sotorasib/adagrasib) × N gene columns.
BODY_W_MAX_CM <- 6     # cap on the gene-axis (was the height cap when vertical)
BODY_H_CM     <- 1.8   # 3 condition rows at ~0.6 cm each

CNV_FILES <- list(
  H23  = file.path(INPUT_DIR, "h23_cnv_snv_sv.rds"),
  CALU = file.path(INPUT_DIR, "calu_cnv_snv_sv.rds"),
  MEP  = file.path(INPUT_DIR, "mep_cnv_snv_sv.rds")
)

# ---- signatures (HUMAN symbols) --------------------------------------------
ras84    <- read.csv(file.path(INPUT_DIR,
                               "RAS_84_PhilipEast_et_al_NatCom_2022.csv"),
                     stringsAsFactors = FALSE)
mapk_erk <- readRDS(file.path(INPUT_DIR,
                              "mapk_erk_signature_East_NatCom_2022.rds"))

ras84_sym <- unique(ras84$HGNC.symbol); ras84_sym <- ras84_sym[!is.na(ras84_sym) & nzchar(ras84_sym)]
mapk_sym  <- unique(mapk_erk$gene_symbol[mapk_erk$gs_name == "GOBP_MAPK_CASCADE"])
mapk_sym  <- mapk_sym[!is.na(mapk_sym) & nzchar(mapk_sym)]

SIGS <- list(
  MAPK       = list(label = "GOBP_MAPK_CASCADE", symbols = mapk_sym),
  RAS84      = list(label = "RAS84",             symbols = ras84_sym),
  RAS_FAMILY = list(label = "RAS family",        symbols = c("KRAS","HRAS","NRAS","MRAS"))
)

# ---- per-line aggregation ---------------------------------------------------
aggregate_line <- function(rds_path) {
  cnv <- as.data.frame(readRDS(rds_path)$CNV)
  cnv <- cnv %>%
    transmute(gene_name = as.character(gene_name),
              condition = as.character(condition),
              category  = as.character(category),
              CNt       = as.integer(CNt)) %>%
    filter(!is.na(gene_name), nzchar(gene_name),
           condition %in% CONDS,
           category %in% c("Amplification", "Loss"))

  cnv %>%
    group_by(gene_name, condition) %>%
    summarise(
      has_amp  = any(category == "Amplification"),
      has_loss = any(category == "Loss"),
      # CNt with the largest deviation from 2 across rows of this gene × cond
      CNt = CNt[which.max(abs(CNt - 2))],
      .groups = "drop")
}

ag_all <- lapply(CNV_FILES, aggregate_line)

# ---- gene-selection helpers ------------------------------------------------
drug_only_genes <- function(ag, sig_syms, extreme = FALSE) {
  d <- ag %>% filter(gene_name %in% sig_syms)
  has_naive <- d %>% filter(condition == "000") %>% pull(gene_name)
  drug <- d %>% filter(condition %in% c("010", "020"))
  if (extreme) drug <- drug %>% filter(CNt < EXTREME_LO | CNt > EXTREME_HI)
  setdiff(unique(drug$gene_name), unique(has_naive))
}

# Build the 3-column CNt matrix (naive / sotorasib / adagrasib) for a given
# gene set and cell line.
build_matrix <- function(ag, genes) {
  m <- matrix(NA_real_, nrow = length(genes), ncol = length(CONDS),
              dimnames = list(genes, unname(COND_LABEL[CONDS])))
  for (k in seq_along(CONDS)) {
    sub <- ag %>% filter(gene_name %in% genes, condition == CONDS[k])
    m[sub$gene_name, k] <- as.numeric(sub$CNt)
  }
  m
}

# ---- heatmap renderer (horizontal layout) ----------------------------------
# Genes on columns, treatments on rows. Cluster columns; rows fixed in
# naive → sotorasib → adagrasib order.
render_heatmap <- function(m, line, sig_label, filter_label, out_pdf) {
  if (!nrow(m)) {
    message(sprintf("  [%s] empty matrix — skipping %s", line, out_pdf))
    return(invisible(NULL))
  }

  mh <- t(m)   # rows = treatments, cols = genes

  cap  <- 6
  ramp <- colorRamp2(c(0, 2, cap),
                     c(CNT_COL_LOW, CNT_COL_MID, CNT_COL_HIGH))

  per_col_cm <- min(0.33, BODY_W_MAX_CM / ncol(mh))
  body_w     <- unit(ncol(mh) * per_col_cm, "cm")
  body_h     <- unit(BODY_H_CM, "cm")
  font_pt    <- max(3, min(8, floor(per_col_cm * 28.35 * 0.7)))

  # NA-tolerant clustering on columns (the gene axis).
  mh_for_clust <- mh; mh_for_clust[is.na(mh_for_clust)] <- 2
  col_dend <- if (ncol(mh_for_clust) >= 2)
    hclust(dist(t(mh_for_clust)), method = "complete") else FALSE

  left_anno <- rowAnnotation(
    treatment = rownames(mh),
    col = list(treatment = COND_COL),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_side = "top",
    width = unit(4, "mm"),
    show_legend = FALSE
  )

  ht <- Heatmap(
    mh,
    name        = "CNt",
    col         = ramp,
    na_col      = CNT_COL_NA,
    height      = body_h,
    width       = body_w,
    cluster_rows         = FALSE,
    cluster_columns      = col_dend,
    show_row_dend        = FALSE,
    show_column_dend     = TRUE,
    show_row_names       = TRUE,
    row_names_side       = "left",
    row_names_gp         = gpar(fontsize = 8),
    show_column_names    = TRUE,
    column_names_side    = "bottom",
    column_names_rot     = 90,
    column_names_gp      = gpar(fontsize = font_pt),
    column_title         = sprintf("%s — %s — %s   (n=%d genes)",
                                    line, sig_label, filter_label, ncol(mh)),
    column_title_gp      = gpar(fontsize = 9, fontface = "bold"),
    left_annotation      = left_anno,
    heatmap_legend_param = list(
      title          = "CNt",
      at             = pretty(c(0, 2, cap), 5),
      title_gp       = gpar(fontsize = 8),
      labels_gp      = gpar(fontsize = 8),
      title_position = "topcenter")
  )

  page_w <- min(BODY_W_MAX_CM / 2.54, (ncol(mh) * per_col_cm) / 2.54) + 3.5
  page_h <- BODY_H_CM / 2.54 + 3.0   # body + treatment strip + col-name labels + title
  grDevices::pdf(out_pdf, width = page_w, height = page_h,
                 useDingbats = FALSE, bg = "white")
  draw(ht, heatmap_legend_side = "right",
       annotation_legend_side = "right", merge_legend = TRUE)
  invisible(dev.off())
  message(sprintf("  Wrote %s   (%.1f × %.1f in)", out_pdf, page_w, page_h))
}

# ---- driver -----------------------------------------------------------------
FILTERS <- list(
  list(key = "",        extreme = FALSE,
       label = "drug-only Amp/Loss (absent in naive)"),
  list(key = "_extreme", extreme = TRUE,
       label = sprintf("drug-only & CNt < %d or > %d", EXTREME_LO, EXTREME_HI))
)

for (sig_name in names(SIGS)) {
  sig <- SIGS[[sig_name]]
  message(sprintf("\n=== %s — signature size = %d ===",
                  sig_name, length(sig$symbols)))
  for (line in names(ag_all)) {
    ag <- ag_all[[line]]
    for (flt in FILTERS) {
      genes <- drug_only_genes(ag, sig$symbols, extreme = flt$extreme)
      message(sprintf("  [%s] %-44s : %d genes",
                      line, flt$label, length(genes)))
      m <- build_matrix(ag, genes)
      out_pdf <- file.path(OUT_DIR,
                           sprintf("heatmap_%s_%s%s.pdf",
                                   sig_name, line, flt$key))
      render_heatmap(m, line, sig$label, flt$label, out_pdf)
    }
  }
}
