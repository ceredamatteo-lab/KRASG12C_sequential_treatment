# =============================================================================
# MEP_mutations_heatmap.R
#
# High-impact somatic mutation heatmap across three MEP conditions:
#   000 = naive, 010 = sotorasib, 020 = adagrasib.
#
# NOTE: MEP's SNV table is the largest of the three cell lines (~493 rows);
# the HIGH-impact union can therefore be substantially larger than H23's.
# The script automatically scales the figure height with the number of rows.
#
# Pipeline:
#   (1) Load the per-cell-line .rds (named list with $SNV).
#   (2) Filter to IMPACT == "HIGH".
#   (3) Build a (mutation x condition) matrix of TUMOR_AF (max-aggregate on
#       rare within-sample duplicates; NA where absent).
#   (4) ComplexHeatmap:
#         - viridis colour ramp on [0, 1] via circlize::colorRamp2
#         - NA cells in grey92
#         - row labels: gene + p.XposY (or splice/frameshift fallback)
#         - column labels: naive / sotorasib / adagrasib
#         - top annotation: untreated vs KRAS G12Ci
#         - rows clustered (binary presence/absence + Ward.D2);
#           columns are NOT clustered (experimental arm order is meaningful)
#   (5) Save <TITLE>_mutations_heatmap.pdf + .png alongside the script.
#
# Run (locally, since the sandbox has no R runtime):
#     cd /path/to/PanRAS
#     Rscript MEP_mutations_heatmap.R
#
# Override defaults from the shell:
#     Rscript MEP_mutations_heatmap.R path/to/mep_cnv_snv_sv.rds out_dir TITLE
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(grid)
})

# -----------------------------------------------------------------------------
# Parameters (edit these at the top of the script, or override via CLI args)
# -----------------------------------------------------------------------------

# Script-dir-relative defaults (matches the circos trio's idiom)
.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) {
    normalizePath(dirname(sub("^--file=", "", f[1])), mustWork = FALSE)
  } else {
    getwd()
  }
})

args       <- commandArgs(trailingOnly = TRUE)
# Project layout: PanRAS/{Script,Figures,input}/.
# .rds inputs live in ../input, PDFs are written to ../Figures.
.project_root <- normalizePath(file.path(.script_dir, "..", ".."),
                               mustWork = FALSE)
INPUT_DIR   <- file.path(.project_root, "input")
FIGURES_DIR <- file.path(.project_root, "Figures")

rds_path   <- if (length(args) >= 1) args[1] else file.path(INPUT_DIR, "mep_cnv_snv_sv.rds")
out_dir    <- if (length(args) >= 2) args[2] else file.path(FIGURES_DIR, "mutations")
TITLE      <- if (length(args) >= 3) args[3] else "MEP"

# Named vector: condition code -> human label for column headers
COND_LABEL <- c("000" = "naive",
                "010" = "sotorasib",
                "020" = "adagrasib")

# Named vector: condition code -> experimental arm (for top annotation)
ARM <- c("000" = "untreated",
         "010" = "KRAS G12Ci",
         "020" = "KRAS G12Ci")

ARM_COL <- c("untreated"  = "#BDBDBD",
             "KRAS G12Ci" = "#4575B4")

# RAS gene family (canonical RAS + R-RAS subfamily). Used by variant 2 so that
# mutations in Hras/Kras/Nras/Rras/Rras2/Mras are kept even if they are not
# flagged by the RAS84 / MAPK-ERK / NCG_cg signatures.
# NOTE: MEP is mouse — gene symbols use mouse casing (matches SYMBOL column).
RAS_FAMILY <- c("Hras", "Kras", "Nras", "Rras", "Rras2", "Mras")

# -----------------------------------------------------------------------------
# Per-gene-set scopes (2026-04-25 extension)
# -----------------------------------------------------------------------------
# Two focused panels rendered after the existing variants:
#   "RAS"         : KRAS / HRAS / MRAS / NRAS only
#   "TP53MYCSRSF" : TP53 / MYC / SRSF1..SRSF12  (mouse TP53 is "Trp53")
# Symbols are matched case-insensitively against gene_sym (resolved earlier
# from SYMBOL -> gene_name -> gene_id_nv). For MEP this catches mouse-cased
# Kras, Hras, Mras, Nras, Trp53, Myc, Srsf1..Srsf12.
GENE_SET_SCOPES <- list(
  RAS = list(
    file_tag = "RAS",
    title    = "RAS family (KRAS/HRAS/MRAS/NRAS)",
    symbols  = c("KRAS", "HRAS", "MRAS", "NRAS")
  ),
  TP53MYCSRSF = list(
    file_tag = "TP53_MYC_SRSF",
    title    = "TP53 / MYC / SRSF1-12",
    symbols  = c("TP53", "TRP53", "MYC", paste0("SRSF", 1:12))
  )
)

# Standard genetic code (used to build HGVSp labels on-the-fly from Codons +
# Protein_position when the VCF doesn't carry a pre-formatted HGVSp column)
GENETIC_CODE <- c(
  TTT="F", TTC="F", TTA="L", TTG="L",
  CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M",
  GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S",
  CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T",
  GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K",
  GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W",
  CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Translate a single codon (case-insensitive; '-' for gap-padded VEP codons)
.translate_codon <- function(codon) {
  if (is.na(codon) || !nzchar(codon)) return(NA_character_)
  codon <- toupper(gsub("-", "", codon))
  if (nchar(codon) != 3) return(NA_character_)
  aa <- GENETIC_CODE[[codon]]
  if (is.null(aa)) NA_character_ else aa
}

# Build a human-readable label per SNV row.
#   1. Protein substitution (preferred):  "<gene> p.<REF><pos><ALT>"
#   2. Splice site:                       "<gene> splice acceptor" / "donor"
#   3. Frameshift:                        "<gene> frameshift"
#   4. Fallback:                          "<gene> <consequence>" / gene alone
#
# Gene-symbol fallback cascade (MEP, mouse): SYMBOL -> gene_name -> gene_id_nv -> "?".
# SYMBOL is the VEP mouse gene symbol (0 NAs in the MEP SNV table); gene_name
# holds uppercase human orthologue symbols and is kept only as a safety net.
nice_label <- function(gene, gene_id, codons, protein_position, consequence) {
  if (is.na(gene) || !nzchar(gene)) {
    gene <- if (!is.na(gene_id) && nzchar(gene_id)) gene_id else "?"
  }

  # Tier 1: protein change from Codons ("aaa/aAa") + Protein_position
  if (!is.na(codons) && grepl("/", codons) &&
      !is.na(protein_position) && nzchar(protein_position)) {
    parts <- strsplit(codons, "/", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      ref_aa <- .translate_codon(parts[1])
      alt_aa <- .translate_codon(parts[2])
      pos    <- sub("[^0-9].*$", "", protein_position)  # keep leading digits
      if (!is.na(ref_aa) && !is.na(alt_aa) && nzchar(pos) &&
          ref_aa != alt_aa) {
        return(sprintf("%s p.%s%s%s", gene, ref_aa, pos, alt_aa))
      }
    }
  }

  # Tier 2: splice-site annotations
  if (!is.na(consequence)) {
    if (grepl("splice_acceptor", consequence))
      return(sprintf("%s splice acceptor", gene))
    if (grepl("splice_donor", consequence))
      return(sprintf("%s splice donor", gene))
    if (grepl("stop_gained", consequence))
      return(sprintf("%s stop gained", gene))
    if (grepl("stop_lost", consequence))
      return(sprintf("%s stop lost", gene))
    if (grepl("start_lost", consequence))
      return(sprintf("%s start lost", gene))
    if (grepl("frameshift", consequence))
      return(sprintf("%s frameshift", gene))
  }

  # Tier 3: fallback to first consequence term
  if (!is.na(consequence) && nzchar(consequence)) {
    first <- strsplit(consequence, "[&,]")[[1]][1]
    return(sprintf("%s %s", gene, first))
  }
  gene
}

# -----------------------------------------------------------------------------
# Load & filter
# -----------------------------------------------------------------------------

message(sprintf("[%s] loading %s", TITLE, rds_path))
dat <- readRDS(rds_path)
if (!"SNV" %in% names(dat))
  stop("rds does not contain a 'SNV' element: ", rds_path)
snv <- as.data.frame(dat$SNV, stringsAsFactors = FALSE)

required_cols <- c("ID", "IMPACT", "TUMOR_AF", "condition",
                   "SYMBOL", "gene_name",
                   "Codons", "Protein_position", "Consequence")
missing <- setdiff(required_cols, colnames(snv))
if (length(missing))
  stop("SNV table is missing required columns: ",
       paste(missing, collapse = ", "))

# Coerce types defensively (some VCF-derived frames store AF as character)
snv$TUMOR_AF  <- suppressWarnings(as.numeric(snv$TUMOR_AF))
snv$condition <- as.character(snv$condition)
snv$IMPACT    <- as.character(snv$IMPACT)

hi <- snv[!is.na(snv$IMPACT) & snv$IMPACT == "HIGH", , drop = FALSE]
message(sprintf("[%s] HIGH-impact rows: %d / %d", TITLE, nrow(hi), nrow(snv)))

if (nrow(hi) == 0)
  stop("No HIGH-impact SNVs in ", rds_path,
       " — heatmap cannot be drawn.")

# Keep only recognised conditions
hi <- hi[hi$condition %in% names(COND_LABEL), , drop = FALSE]
if (nrow(hi) == 0)
  stop("No HIGH-impact rows with condition in {000,010,020}.")

# Drop long indels: parse REF/ALT from the ID (<CHROM>_<POS>_<REF>_<ALT>_<gene>)
# and keep rows where both alleles are shorter than 5 nucleotides.
.id_parts   <- regmatches(hi$ID,
                          regexec("^(chr[^_]+)_([0-9]+)_([^_]+)_([^_]+)_(.+)$",
                                  hi$ID))
.ref_allele <- vapply(.id_parts,
                      function(x) if (length(x) >= 5L) x[4L] else NA_character_,
                      character(1))
.alt_allele <- vapply(.id_parts,
                      function(x) if (length(x) >= 5L) x[5L] else NA_character_,
                      character(1))
short_allele <- !is.na(.ref_allele) & !is.na(.alt_allele) &
                nchar(.ref_allele) < 5L & nchar(.alt_allele) < 5L
message(sprintf(
  "[%s] allele-length filter (REF & ALT < 5 nt): %d / %d kept (%d dropped)",
  TITLE, sum(short_allele), nrow(hi), sum(!short_allele)))
hi <- hi[short_allele, , drop = FALSE]
if (nrow(hi) == 0)
  stop("No HIGH-impact rows remain after the allele-length filter.")

# -----------------------------------------------------------------------------
# Build mutation x condition matrix (TUMOR_AF; NA = absent)
# -----------------------------------------------------------------------------

cond_order <- c("000", "010", "020")
cond_order <- cond_order[cond_order %in% hi$condition]

muts <- sort(unique(hi$ID))
mat  <- matrix(NA_real_,
               nrow = length(muts),
               ncol = length(cond_order),
               dimnames = list(muts, cond_order))

# max-aggregate rare intra-condition duplicates
for (i in seq_len(nrow(hi))) {
  m <- hi$ID[i]; c <- hi$condition[i]; v <- hi$TUMOR_AF[i]
  if (is.na(v)) next
  prev <- mat[m, c]
  mat[m, c] <- if (is.na(prev)) v else max(prev, v)
}

# -----------------------------------------------------------------------------
# Row labels (one per mutation ID; preserves mat row order)
# -----------------------------------------------------------------------------

# Pick the first SNV row per ID to source the label fields
first_per_id <- hi[!duplicated(hi$ID), , drop = FALSE]
rownames(first_per_id) <- first_per_id$ID

row_labels <- vapply(muts, function(id) {
  r <- first_per_id[id, ]
  # MEP is mouse: prefer SYMBOL (mouse-cased), fall back to gene_name, then gene_id_nv.
  g <- r$SYMBOL
  if (is.na(g) || !nzchar(g)) g <- r$gene_name
  if (is.na(g) || !nzchar(g)) g <- r$gene_id_nv
  nice_label(g, r$gene_id_nv,
             r$Codons, r$Protein_position, r$Consequence)
}, character(1))

# Disambiguate exact label collisions (shouldn't happen for H23, but safe)
if (anyDuplicated(row_labels)) {
  dup <- duplicated(row_labels) | duplicated(row_labels, fromLast = TRUE)
  row_labels[dup] <- paste0(row_labels[dup], " [", muts[dup], "]")
}

# -----------------------------------------------------------------------------
# Retention filter:
#   keep HIGH SNVs that either (a) carry an Existing_variation annotation
#   (dbSNP / COSMIC / etc.) or (b) are detected in more than one sample
#   (i.e. in >= 2 of the three conditions).
# -----------------------------------------------------------------------------

ex_raw       <- first_per_id[muts, "Existing_variation"]
has_existing <- !is.na(ex_raw) & nzchar(ex_raw) & ex_raw != "-"
n_samples    <- rowSums(!is.na(mat))
keep         <- has_existing | n_samples > 1L

message(sprintf(
  "[%s] retained %d / %d HIGH mutations (Existing_variation=%d, >1 sample=%d)",
  TITLE, sum(keep), length(muts),
  sum(has_existing), sum(n_samples > 1L)))

if (sum(keep) == 0L)
  stop(sprintf("[%s] no HIGH mutations pass the retention filter.", TITLE))

mat        <- mat[keep, , drop = FALSE]
row_labels <- row_labels[keep]
muts       <- muts[keep]

# -----------------------------------------------------------------------------
# Manual row order:
#   (A) rows with a value in 000, sorted by 000 TUMOR_AF (descending)
#   (B) rows absent from 000 but present in 010 and/or 020,
#       sorted by max(010, 020) TUMOR_AF (descending)
# Row clustering is therefore disabled.
# -----------------------------------------------------------------------------

has_000_col <- "000" %in% cond_order
col_010_020 <- intersect(c("010", "020"), cond_order)

if (has_000_col) {
  in_000  <- !is.na(mat[, "000"])
  idx_top <- which(in_000)[order(-mat[in_000, "000"])]
  idx_bot <- which(!in_000)
} else {
  idx_top <- integer(0)
  idx_bot <- seq_len(nrow(mat))
}

if (length(idx_bot) && length(col_010_020)) {
  bot_vals <- apply(mat[idx_bot, col_010_020, drop = FALSE], 1,
                    function(x) if (all(is.na(x))) -Inf
                                else max(x, na.rm = TRUE))
  idx_bot  <- idx_bot[order(-bot_vals)]
}

row_order <- c(idx_top, idx_bot)

# -----------------------------------------------------------------------------
# Row-level annotations (T/F per retained mutation) drawn from the SNV "is*"
# columns. A row is "rest" (= FALSE) when the column value is missing, empty,
# or the literal string "rest"; anything else is TRUE.
# -----------------------------------------------------------------------------

.is_true <- function(x) !is.na(x) & nzchar(x) & tolower(x) != "rest"

ann_ras84    <- .is_true(as.character(first_per_id[muts, "isRAS84"]))
ann_mapk_erk <- .is_true(as.character(first_per_id[muts, "isMAPK_ERK_signature"]))
ann_ncg_cg   <- .is_true(as.character(first_per_id[muts, "isNCG_cg"]))

# Per-row gene symbol (MEP is mouse; use SYMBOL -> gene_name -> gene_id_nv).
# Used to flag rows in the RAS gene family (matched case-sensitively against
# mouse-cased RAS_FAMILY above).
gene_sym <- vapply(muts, function(id) {
  r <- first_per_id[id, ]
  g <- r$SYMBOL
  if (is.na(g) || !nzchar(g)) g <- r$gene_name
  if (is.na(g) || !nzchar(g)) {
    g <- if (!is.na(r$gene_id_nv) && nzchar(r$gene_id_nv)) r$gene_id_nv
         else NA_character_
  }
  g
}, character(1))

in_ras_family <- !is.na(gene_sym) & gene_sym %in% RAS_FAMILY

# Apply the computed row order to the matrix and all row-aligned vectors
mat           <- mat[row_order, , drop = FALSE]
row_labels    <- row_labels[row_order]
ann_ras84     <- ann_ras84[row_order]
ann_mapk_erk  <- ann_mapk_erk[row_order]
ann_ncg_cg    <- ann_ncg_cg[row_order]
gene_sym      <- gene_sym[row_order]
in_ras_family <- in_ras_family[row_order]

# -----------------------------------------------------------------------------
# Heatmap rendering (helper called twice: full + without frameshift/splice)
# -----------------------------------------------------------------------------

col_fun <- colorRamp2(seq(0, 1, length.out = 256), viridis(256))

top_anno <- HeatmapAnnotation(
  treatment = ARM[cond_order],
  col = list(treatment = ARM_COL),
  annotation_legend_param = list(treatment = list(title = "Arm")),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

BOOL_COL <- c("yes" = "#D7191C", "no" = "grey92")       # RAS84
BOOL_COL2 <- c("yes" = "#2C7FB8", "no" = "grey92")      # MAPK/ERK
BOOL_COL3 <- c("yes" = "#1A9641", "no" = "grey92")      # NCG_cg

render_heatmap <- function(m, labs, a_ras84, a_mapk_erk, a_ncg_cg,
                           title_suffix, file_suffix) {

  right_anno <- rowAnnotation(
    RAS84      = ifelse(a_ras84,    "yes", "no"),
    `MAPK/ERK` = ifelse(a_mapk_erk, "yes", "no"),
    NCG_cg     = ifelse(a_ncg_cg,   "yes", "no"),
    col = list(
      RAS84      = BOOL_COL,
      `MAPK/ERK` = BOOL_COL2,
      NCG_cg     = BOOL_COL3
    ),
    annotation_name_side = "top",
    annotation_name_gp   = gpar(fontsize = 8),
    annotation_name_rot  = 90,
    simple_anno_size     = unit(3, "mm"),
    show_legend          = c(RAS84 = TRUE, `MAPK/ERK` = TRUE, NCG_cg = TRUE),
    annotation_legend_param = list(
      RAS84      = list(title = "RAS84"),
      `MAPK/ERK` = list(title = "MAPK/ERK"),
      NCG_cg     = list(title = "NCG_cg")
    )
  )

  col_title <- sprintf(
    "%s high-impact somatic mutations \u2014 KRAS G12C inhibitor response (TUMOR_AF)%s",
    TITLE,
    if (nzchar(title_suffix)) paste0("\n", title_suffix) else "")

  ht <- Heatmap(
    m,
    name               = "TUMOR_AF",
    col                = col_fun,
    na_col             = "grey92",
    rect_gp            = gpar(col = "grey80", lwd = 0.2),
    cluster_rows       = FALSE,
    cluster_columns    = FALSE,
    row_labels         = labs,
    row_names_side     = "left",
    row_names_gp       = gpar(fontsize = 8),
    column_labels      = COND_LABEL[cond_order],
    column_names_rot   = 90,
    column_names_gp    = gpar(fontsize = 10, fontface = "bold"),
    width              = unit(10 * ncol(m), "mm"),
    top_annotation     = top_anno,
    right_annotation   = right_anno,
    column_title       = col_title,
    column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title        = "TUMOR_AF",
      at           = c(0, 0.25, 0.5, 0.75, 1),
      legend_height = unit(30, "mm")
    ),
    show_row_dend      = FALSE,
    border             = TRUE
  )

  n_rows    <- nrow(m)
  width_mm  <- 180
  height_mm <- max(120, 3 * n_rows + 70)

  pdf_path <- file.path(
    out_dir, sprintf("%s_mutations_heatmap%s.pdf", TITLE, file_suffix))

  pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4)
  draw(ht, merge_legend = TRUE,
       heatmap_legend_side     = "right",
       annotation_legend_side  = "right")
  invisible(dev.off())

  message(sprintf("[%s] wrote %s  (%d mutations x %d conditions)%s",
                  TITLE, pdf_path, n_rows, length(cond_order),
                  if (nzchar(title_suffix)) paste0("  [", title_suffix, "]") else ""))
}

# -----------------------------------------------------------------------------
# Output 1: full heatmap
# Output 2: heatmap restricted to mutations in genes annotated RAS84 / MAPK-ERK
#           / NCG_cg, or in the RAS gene family (HRAS/KRAS/NRAS/RRAS/RRAS2/MRAS)
# -----------------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

render_heatmap(
  mat, row_labels,
  ann_ras84, ann_mapk_erk, ann_ncg_cg,
  title_suffix = "",
  file_suffix  = "")

keep2 <- ann_ras84 | ann_mapk_erk | ann_ncg_cg | in_ras_family
message(sprintf(
  "[%s] variant 2 (RAS84 | MAPK-ERK | NCG_cg | RAS family): keeping %d / %d rows (RAS84=%d, MAPK-ERK=%d, NCG_cg=%d, RAS-family=%d)",
  TITLE, sum(keep2), length(keep2),
  sum(ann_ras84), sum(ann_mapk_erk), sum(ann_ncg_cg), sum(in_ras_family)))

if (sum(keep2) >= 1L) {
  render_heatmap(
    mat[keep2, , drop = FALSE],
    row_labels[keep2],
    ann_ras84[keep2], ann_mapk_erk[keep2], ann_ncg_cg[keep2],
    title_suffix = "RAS84 / MAPK-ERK / NCG_cg signature or RAS gene family",
    file_suffix  = "_sig_or_ras")
} else {
  message(sprintf(
    "[%s] variant 2 skipped: no rows match RAS84/MAPK-ERK/NCG_cg/RAS-family.",
    TITLE))
}

# -----------------------------------------------------------------------------
# Per-gene-set heatmaps (2026-04-25 extension)
# Filters retained mutations to genes in each scope (case-insensitive match
# against gene_sym, which already resolved from SYMBOL -> gene_name ->
# gene_id_nv). Mouse symbols are upper-cased before matching, so Trp53 / Kras
# / Srsf1..12 hit TP53 / KRAS / SRSF1..12 in the gene set.
# -----------------------------------------------------------------------------

for (gs_name in names(GENE_SET_SCOPES)) {
  gs   <- GENE_SET_SCOPES[[gs_name]]
  syms_upper <- toupper(gs$symbols)

  in_gs <- !is.na(gene_sym) & toupper(gene_sym) %in% syms_upper
  message(sprintf(
    "[%s] gene-set %-12s: %d / %d retained mutations",
    TITLE, gs_name, sum(in_gs), length(in_gs)))

  if (sum(in_gs) == 0L) {
    message(sprintf("[%s] gene-set %s skipped: no retained mutations.",
                    TITLE, gs_name))
    next
  }

  render_heatmap(
    mat[in_gs, , drop = FALSE],
    row_labels[in_gs],
    ann_ras84[in_gs], ann_mapk_erk[in_gs], ann_ncg_cg[in_gs],
    title_suffix = sprintf("Gene set: %s", gs$title),
    file_suffix  = paste0("_", gs$file_tag))
}
