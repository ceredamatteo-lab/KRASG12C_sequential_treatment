#!/usr/bin/env Rscript
# =============================================================================
# MEP_circos.R
# Publication-quality circos plot of the MEP genomic landscape.
# Radial layout (outside → inside):
#   OUTSIDE gene labels : genes from the narrow set hit by CNV (Amp or Loss)
#   ideogram            : mouse (mm39) chromosome ideogram (names only)
#   CNV                 : half-height CNt bars, coloured by category, capped 95th pctile
#   INSIDE  gene labels : genes from the narrow set hit by HIGH-impact SV
#   SV                  : DEL/DUP tick track (half height)
#   BND arcs            : innermost (line width ∝ TUMOR_AF)
#
# Narrow gene set (050-000 baseline, union):
# (ALL members shown only if hit by Amplification | Loss | HIGH-impact SV;
#  labels use the Title-case mouse symbols from the CNV column mouse_gene_name,
#  e.g. Kras, Hras, Nras, Mras, Myc, Trp53):
#   - Myc, Trp53 (explicit callouts)
#   - RAS family (Kras, Hras, Nras, Mras)
#   - every RAS84-signature gene
#   (exact counts printed to stdout at runtime; split into CNV-hit / SV-hit)
#
# Dependencies: R >= 4.x, circlize >= 0.4.16, dplyr, tibble, purrr, stringr.
# Gene coordinates come from the CNV table; no TxDb / org.Mm.eg.db needed.
# Genome: mouse (GRCm39 / mm39). Sibling H23 & Calu scripts use human hg38.
#
# Usage (defaults reproduce the MEP figure when script + .rds are co-located):
#   Rscript MEP_circos.R                      # all defaults, CWD-independent
#   Rscript MEP_circos.R /path/to/mep_cnv_snv_sv.rds 000 MEP /path/to/outdir
#
#   arg 1 : path to the per-cell-model .rds (default: <scriptdir>/mep_cnv_snv_sv.rds)
#   arg 2 : condition code to filter to     (default: "000" — baseline)
#   arg 3 : output prefix + display title   (default: "MEP")
#   arg 4 : output directory                (default: <scriptdir>)
#
# <scriptdir> is resolved from the Rscript --file= argument, so the defaults
# work no matter which directory you invoke Rscript from.
#
# Sibling scripts (same schema, same plotting logic): H23_circos.R, Calu_circos.R.
# NOTE: as of 2026-04-22 mep_cnv_snv_sv.rds in the PanRAS folder is byte-identical
#       to calu_cnv_snv_sv.rds (MD5 77552785fac866870948e4786ea097cc). Replace it
#       with the real MEP data before running, otherwise this script will render
#       a figure matching Calu-1.
# This script does NOT modify the input .rds.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(circlize)
})

# --------------------------- 0. Arguments -----------------------------------
# Resolve the directory this script lives in, so default paths don't depend
# on the caller's working directory.
.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})

# Project layout: PanRAS/{Script,Figures,input}/.
# .rds inputs live in ../input, PDFs are written to ../Figures.
.project_root <- normalizePath(file.path(.script_dir, "..", ".."),
                               mustWork = FALSE)
INPUT_DIR   <- file.path(.project_root, "input")
FIGURES_DIR <- file.path(.project_root, "Figures")

args       <- commandArgs(trailingOnly = TRUE)
rds_path   <- if (length(args) >= 1) args[1] else file.path(INPUT_DIR,
                                                            "mep_cnv_snv_sv.rds")
cond_code  <- if (length(args) >= 2) args[2] else "000"
sample_ttl <- if (length(args) >= 3) args[3] else "MEP"
out_dir    <- if (length(args) >= 4) args[4] else file.path(FIGURES_DIR, "circos")

stopifnot(file.exists(rds_path))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message("Input  : ", rds_path)
message("Cond.  : ", cond_code)
message("Title  : ", sample_ttl)
message("Outdir : ", out_dir)

# --------------------------- 1. Load & validate -----------------------------
dat <- readRDS(rds_path)
if (!all(c("CNV","SNV","SV") %in% names(dat)))
  stop("Expected list elements CNV, SNV, SV — got: ",
       paste(names(dat), collapse = ", "))
cnv_raw <- as_tibble(dat$CNV); sv_raw <- as_tibble(dat$SV, rownames = "._rn")

# --------------------------- 2. Filter to chosen condition -------------------
cnv_cond <- cnv_raw %>% filter(as.character(condition) == cond_code)
if (nrow(cnv_cond) == 0)
  stop("No CNV rows with condition == '", cond_code, "'")

# Derive the sample barcode for this condition from CNV (authoritative)
barcode <- unique(as.character(cnv_cond$barcode))
if (length(barcode) != 1)
  stop("Expected a single CNV barcode for condition ", cond_code,
       ", got: ", paste(barcode, collapse = ", "))
message("Barcode: ", barcode)

# SV table has no sample column — infer via row-name prefix (== barcode).
sv_cond <- sv_raw %>% filter(startsWith(._rn, .env$barcode)) %>% select(-._rn)
if (nrow(sv_cond) == 0)
  warning("No SV rows match barcode prefix '", barcode,
          "'. SV track will be empty.")

# --------------------------- 3. CNV: collapse to segments --------------------
# The CNV table is exploded one row per overlapping gene.
# A segment is uniquely defined by (chromosome, start.pos, end.pos) within a
# barcode; we collapse to one row per segment.
seg <- cnv_cond %>%
  distinct(chromosome, start.pos, end.pos, .keep_all = TRUE) %>%
  transmute(chr      = as.character(chromosome),
            start    = as.integer(start.pos),
            end      = as.integer(end.pos),
            CNt      = as.integer(CNt),
            category = factor(as.character(category),
                              levels = c("Amplification","Loss","LOH","Neutral")))

# Drop unplaceable contigs (keep canonical mouse autosomes + chrX + chrY)
canon <- c(paste0("chr", 1:19), "chrX", "chrY")
seg <- seg %>% filter(chr %in% canon)

# Cap CNt at 99th percentile so one extreme segment doesn't flatten the rest
cnt_cap <- ceiling(quantile(seg$CNt, 0.95, na.rm = TRUE))
if (!is.finite(cnt_cap) || cnt_cap < 4) cnt_cap <- 6
seg <- seg %>% mutate(CNt_plot = pmin(CNt, cnt_cap))

# Colour-blind-safe palette (Okabe–Ito), semantically consistent with direction
cnv_pal <- c(Amplification = "#D55E00",   # vermillion — gain
             Loss          = "#0072B2",   # blue        — loss
             LOH           = "#CC79A7",   # pink        — LOH
             Neutral       = "#BBBBBB")   # neutral grey
seg$col <- cnv_pal[as.character(seg$category)]

# --------------------------- 4. SV: parse breakpoints ------------------------
# Expected schema (inferred from data):
#   SVTYPE == "BND" : VCF breakend. Mate in ALT as  X[chr:pos[ | X]chr:pos] |
#                     [chr:pos[X | ]chr:pos]X . Parsed via regex; one row per SV
#                     (NOT two mate rows), so no mate-record dedup.
#   SVTYPE == "DEL" : ALT == "<DEL>".         No END/SVLEN column present.
#   SVTYPE == "DUP" : ALT == "<DUP:TANDEM>".  No END/SVLEN column present.
# DEL/DUP therefore cannot be drawn as arcs — rendered as ticks on SV track.
parse_mate <- function(alt) {
  m <- str_match(alt, "[\\[\\]]([^:\\[\\]]+):([0-9]+)[\\[\\]]")
  tibble(mate_chr = m[, 2],
         mate_pos = suppressWarnings(as.integer(m[, 3])))
}

sv_u <- sv_cond %>%
  filter(IMPACT == "HIGH") %>%
  distinct(CHROM, POS, REF, ALT, .keep_all = TRUE) %>%
  mutate(POS = as.integer(POS), TUMOR_AF = as.numeric(TUMOR_AF)) %>%
  filter(CHROM %in% canon)

bnd <- sv_u %>% filter(SVTYPE == "BND") %>%
  bind_cols(parse_mate(.$ALT)) %>%
  filter(!is.na(mate_chr), !is.na(mate_pos), mate_chr %in% canon) %>%
  mutate(pair = pmap_chr(list(CHROM, POS, mate_chr, mate_pos),
                        function(a, b, c, d)
                          paste(sort(c(paste(a,b,sep=":"),
                                       paste(c,d,sep=":"))),
                                collapse = "|"))) %>%
  distinct(pair, .keep_all = TRUE)

tick <- sv_u %>% filter(SVTYPE %in% c("DEL","DUP"))

# --------------------------- 4b. Signature-gene pipeline --------------------
# Three signature flags (string-coded: "rest" = FALSE, signature-name = TRUE).
# Gene-level coordinates live in the CNV table's chr / start / end columns
# (distinct from segment coords chromosome / start.pos / end.pos), already mm39.
# Plotted set (this script): all of these are shown ONLY when hit by
# Amplification, Loss, or HIGH-impact SV:
#   (i)   explicit callouts: Myc, Trp53
#   (ii)  RAS family: Kras, Hras, Nras, Mras
#   (iii) RAS84-signature genes
# NOTE: MEP is mouse-genome (mm39, chr1-19+X+Y). Both tables carry mouse
# Title-case symbols natively — CNV.mouse_gene_name and SV.SYMBOL — so we
# key everything on those (no HUMAN→mouse translation needed). The HUMAN
# orthologue stays available in CNV.gene_name / SV.gene_name if useful.
SIG_COLS       <- c("isRAS84", "isMAPK_signature", "isERK_signature")
SIG_TRUE_VAL   <- c(isRAS84          = "RAS84",
                    isMAPK_signature = "MAPK_signature",
                    isERK_signature  = "ERK_signature")
EXPLICIT_GENES <- c("Myc","Trp53")
RAS_FAMILY     <- c("Kras","Hras","Nras","Mras")

# (a) Per-gene signature membership (from the CNV table — its gene annotation
# is the most comprehensive — collapsed to one row per gene). For MEP we key
# on the Title-case mouse symbol (mouse_gene_name) so plot labels match mouse
# nomenclature; the field is renamed to gene_name for downstream compatibility.
gene_sig <- cnv_cond %>%
  filter(!is.na(mouse_gene_name)) %>%
  transmute(gene_name = mouse_gene_name,
            chr_g   = as.character(chr),
            start_g = as.integer(start),
            end_g   = as.integer(end),
            sig_RAS84          = isRAS84          == SIG_TRUE_VAL["isRAS84"],
            sig_MAPK_signature = isMAPK_signature == SIG_TRUE_VAL["isMAPK_signature"],
            sig_ERK_signature  = isERK_signature  == SIG_TRUE_VAL["isERK_signature"]) %>%
  group_by(gene_name, chr_g, start_g, end_g) %>%
  summarise(across(starts_with("sig_"), any), .groups = "drop") %>%
  filter(chr_g %in% canon)

# (b) Hit-type sets (restricted to the 050-000 subset). Both CNV and SV now
# carry mouse Title-case symbols natively (CNV.mouse_gene_name and SV.SYMBOL),
# so no HUMAN→mouse translation is needed.
amp_genes  <- cnv_cond %>% filter(category == "Amplification") %>% pull(mouse_gene_name) %>% unique()
loss_genes <- cnv_cond %>% filter(category == "Loss")          %>% pull(mouse_gene_name) %>% unique()
sv_genes   <- sv_cond  %>% filter(IMPACT   == "HIGH")          %>% pull(SYMBOL)          %>% unique()
amp_genes  <- amp_genes [!is.na(amp_genes )]
loss_genes <- loss_genes[!is.na(loss_genes)]
sv_genes   <- sv_genes  [!is.na(sv_genes  )]

# (c) Plotted-gene union: every callout-list gene must be hit by
#       Amp | Loss | HIGH-SV; nothing is force-plotted.
hit_genes      <- c(amp_genes, loss_genes, sv_genes)
ras84_hit      <- gene_sig %>%
  filter(sig_RAS84) %>%
  filter(gene_name %in% hit_genes) %>%
  pull(gene_name)
ras_family_hit <- intersect(RAS_FAMILY,     hit_genes)
explicit_hit   <- intersect(EXPLICIT_GENES, hit_genes)
plotted_genes  <- unique(c(explicit_hit, ras84_hit, ras_family_hit))

# (d) Assemble the gene-track data frame (one row per plotted gene)
gene_bed <- gene_sig %>%
  filter(gene_name %in% plotted_genes) %>%
  mutate(
    hit_amp  = gene_name %in% amp_genes,
    hit_loss = gene_name %in% loss_genes,
    hit_sv   = gene_name %in% sv_genes
  ) %>%
  # Report any requested gene we couldn't place
  { missing <- setdiff(plotted_genes, .$gene_name)
    if (length(missing) > 0)
      warning("Could not resolve gene coordinates for: ",
              paste(missing, collapse = ", "))
    . } %>%
  arrange(match(chr_g, canon), start_g)

# (e) Split labels by hit-type:
#     outside the ideogram  = CNV hits (Amplification or Loss)
#     inside  the ideogram  = HIGH-impact SV hits
# A gene with both hit-types is drawn on both sides (rare, but handled).
outer_bed <- gene_bed %>% filter(hit_amp | hit_loss)
inner_bed <- gene_bed %>% filter(hit_sv)

# --------------------------- Manual MEP overrides ----------------------------
# (1) Hras and Nras have CNV-table calls in this condition but the user wants
#     them shown only on the SV side (drug-acquired structural events). Strip
#     them from the CNV outside layer.
# (2) Nras has no SV row in the SV table for this condition, but external
#     evidence flags it as a DEL — synthesise a single-bp DEL tick at the
#     gene's coords and add Nras to inner_bed so it gets a label on the
#     SV ring next to the existing Hras BND.
outer_bed <- outer_bed %>% filter(!gene_name %in% c("Hras", "Nras"))

nras_coords <- gene_sig %>% filter(gene_name == "Nras")
if (nrow(nras_coords) >= 1) {
  nras <- nras_coords[1, ]
  manual_del <- data.frame(
    CHROM    = nras$chr_g,
    POS      = nras$start_g,
    ALT      = "<DEL>",
    SVTYPE   = "DEL",
    IMPACT   = "HIGH",
    TUMOR_AF = NA_real_,
    SYMBOL   = "Nras",
    gene_name = "NRAS",
    stringsAsFactors = FALSE
  )
  for (col in setdiff(colnames(tick), colnames(manual_del)))
    manual_del[[col]] <- NA
  tick <- bind_rows(tick, manual_del[, colnames(tick), drop = FALSE])

  nras_inner <- nras
  nras_inner$hit_amp  <- FALSE
  nras_inner$hit_loss <- FALSE
  nras_inner$hit_sv   <- TRUE
  inner_bed <- bind_rows(inner_bed, nras_inner)
  message("[MEP override] Hras/Nras stripped from CNV layer; Nras added as DEL tick.")
}

# --------------------------- 5. Diagnostics ----------------------------------
cat(sprintf("\n%s (%s) diagnostics\n", sample_ttl, cond_code))
cat(sprintf("  CNV: %d gene-rows -> %d unique segments\n",
            nrow(cnv_cond), nrow(seg)))
cat("       categories: ",
    paste(sprintf("%s=%d", levels(seg$category),
                  table(seg$category)[levels(seg$category)]),
          collapse = ", "), "\n")
cat(sprintf("       CNt cap (95th pctile): %d   max CNt: %d\n",
            cnt_cap, max(seg$CNt)))
cat(sprintf("  SV : %d rows -> %d unique variants  (all HIGH: %d)\n",
            nrow(sv_cond), nrow(sv_u),
            sum(sv_cond$IMPACT == "HIGH")))
cat(sprintf("       BND pairs (arcs): %d    DEL ticks: %d    DUP ticks: %d\n",
            nrow(bnd), sum(tick$SVTYPE == "DEL"), sum(tick$SVTYPE == "DUP")))
# Signature-gene diagnostics
cat(sprintf("  Genes: %d plotted (explicit %d + RAS84-hit %d, unioned)\n",
            nrow(gene_bed), length(EXPLICIT_GENES), length(ras84_hit)))
cat(sprintf("       signature membership: RAS84=%d  MAPK_signature=%d  ERK_signature=%d\n",
            sum(gene_bed$sig_RAS84),
            sum(gene_bed$sig_MAPK_signature),
            sum(gene_bed$sig_ERK_signature)))
cat(sprintf("       hit-type: Amp=%d  Loss=%d  HIGH-SV=%d  (multi-hit possible)\n",
            sum(gene_bed$hit_amp),
            sum(gene_bed$hit_loss),
            sum(gene_bed$hit_sv)))
cat(sprintf("       labels: outside (CNV) = %d, inside (SV) = %d  (overlap = %d)\n",
            nrow(outer_bed), nrow(inner_bed),
            length(intersect(outer_bed$gene_name, inner_bed$gene_name))))
cat(sprintf("       outside-labelled (CNV) genes:\n         %s\n",
            paste(strwrap(paste(sort(outer_bed$gene_name), collapse = ", "),
                          width = 78), collapse = "\n         ")))
cat(sprintf("       inside-labelled  (SV ) genes:\n         %s\n",
            paste(strwrap(paste(sort(inner_bed$gene_name), collapse = ", "),
                          width = 78), collapse = "\n         ")))

# --------------------------- 6. mm39 ideogram --------------------------------
# Mouse GRCm39 (mm39) canonical chromosome sizes — used as a reliable offline
# reference. Source: Ensembl/UCSC GRCm39 primary assembly.
mm39_sizes <- tibble(
  chr   = canon,
  start = 1L,
  end   = c(195154279L,181755017L,159745316L,156860686L,151758149L,
            149588044L,144995196L,130127694L,124359700L,130530862L,
            121973369L,120092757L,120883175L,125139656L,104073951L,
             98008968L, 95294699L, 90720763L, 61420004L,
            169476592L, 91455967L)
)

# --------------------------- 7. Plotter --------------------------------------
af_to_lwd <- function(af) pmax(0.175, pmin(1.5, 0.175 + 1.325 * af))

draw_circos <- function(outfile, device = "pdf",
                        width_mm = 60, dpi = 300, pointsize = 8) {
  # PDF-only (PNG emission intentionally removed).
  stopifnot(device == "pdf")
  width_in <- width_mm / 25.4
  grDevices::pdf(outfile, width = width_in, height = width_in,
                 useDingbats = FALSE, bg = "white",
                 pointsize = pointsize)
  on.exit(grDevices::dev.off(), add = TRUE)

  par(mar = c(1, 1, 2, 1), xpd = NA)

  circos.clear()
  circos.par(gap.after      = c(rep(1.2, length(canon) - 1), 10),
             start.degree   = 90,
             cell.padding   = c(0, 0, 0, 0),
             track.margin   = c(0.004, 0.004))

  # Helper: gene midpoint from a bed-ish data frame
  .gene_mid <- function(df) as.integer((df$start_g + df$end_g) %/% 2L)

  # Build the label-track data frames once (used below for outside + inside)
  .mk_label_bed <- function(df) {
    if (nrow(df) == 0) return(df[0, ])
    data.frame(chr   = df$chr_g,
               start = .gene_mid(df),
               end   = .gene_mid(df) + 1L,
               gene  = df$gene_name,
               stringsAsFactors = FALSE)
  }
  # Colour cascade (priority: SV > Amp > Loss > Explicit > default):
  #   SV-hit              → Okabe-Ito green (#009E73)
  #   Amplification-hit   → CNV palette vermillion (#D55E00)
  #   Loss-hit            → CNV palette blue      (#0072B2)
  #   explicit (Myc/Trp53)            → bold via .mk_font (hit-colour kept)
  #   else                → grey25
  # Tolerate NA in logical hit columns.
  .true_v  <- function(x) !is.na(x) & as.logical(x)
  .mk_col <- function(df) {
    out <- rep("grey25", nrow(df))
    out[df$gene_name %in% EXPLICIT_GENES] <- "black"
    out[.true_v(df$hit_loss)]             <- cnv_pal[["Loss"]]
    out[.true_v(df$hit_amp)]              <- cnv_pal[["Amplification"]]
    out[df$gene_name %in% sv_genes]       <- "#009E73"
    out
  }
  .mk_font <- function(df) ifelse(df$gene_name %in% EXPLICIT_GENES, 2L, 1L)

  # Initialize sectors using offline mm39 sizes; no karyotype/cytoband drawn.
  circos.genomicInitialize(as.data.frame(mm39_sizes), plotType = NULL)

  # -------- (A) Outside gene labels (ALL plotted genes: CNV + SV hits) --------
  # Combined outer_bed ∪ inner_bed, deduped by gene_name so genes with
  # multiple hit-types appear once on the outside.
  outer_all <- dplyr::distinct(dplyr::bind_rows(outer_bed, inner_bed),
                               gene_name, .keep_all = TRUE)
  if (nrow(outer_all) > 0) {
    circos.genomicLabels(.mk_label_bed(outer_all),
                         labels.column = 4,
                         side          = "outside",
                         cex           = 0.53,
                         col           = .mk_col(outer_all),
                         font          = .mk_font(outer_all),
                         line_col      = "grey60",
                         line_lwd      = 0.4,
                         padding       = 0.4,
                         connection_height = convert_height(1.25, "mm"),
                         labels_height     = convert_height(5.5, "mm"))
  }

  # -------- (C) Chromosome-name track (no coordinates, no karyotype) --------
  circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
               panel.fun = function(x, y) {
                 chr <- CELL_META$sector.index
                 circos.text(CELL_META$xcenter, 0.5,
                             gsub("^chr", "", chr),
                             facing = "inside", cex = 0.53, niceFacing = TRUE)
               })

  # -------- (D) CNV track (half-height; Okabe–Ito; 95th-pctile cap) --------
  cnv_bed <- as.data.frame(seg[, c("chr","start","end","CNt_plot","col","category")])
  circos.genomicTrack(
    cnv_bed, ylim = c(0, cnt_cap),
    track.height  = 0.10,
    bg.border     = "grey90", bg.lwd = 0.3,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
                         ytop     = value$CNt_plot,
                         ybottom  = 0,
                         col      = value$col,
                         border   = NA)
      circos.lines(CELL_META$cell.xlim, c(2, 2),
                   col = "grey40", lwd = 0.4, lty = 2)
    })
  # y-axis, drawn once on the first (top) sector only
  circos.yaxis(side          = "left",
               at            = unique(c(0, 2, pretty(c(0, cnt_cap), 3), cnt_cap)),
               sector.index  = canon[1],
               labels.cex    = 0.53)

  # -------- (F) SV track (DEL/DUP ticks; halved again) --------
  circos.track(ylim          = c(0, 1),
               track.height  = 0.03,
               bg.border     = "grey90", bg.lwd = 0.3,
               panel.fun = function(x, y) {
                 chr <- CELL_META$sector.index
                 tks <- tick[tick$CHROM == chr, , drop = FALSE]
                 if (nrow(tks) > 0) {
                   col_tick <- ifelse(tks$SVTYPE == "DEL", "#0072B2", "#D55E00")
                   circos.segments(tks$POS, 0.10, tks$POS, 0.90,
                                   col = col_tick, lwd = 0.8)
                 }
               })

  # -------- BND arcs (innermost) --------
  if (nrow(bnd) > 0) {
    for (i in seq_len(nrow(bnd))) {
      circos.link(sector.index1 = bnd$CHROM[i],
                  point1        = bnd$POS[i],
                  sector.index2 = bnd$mate_chr[i],
                  point2        = bnd$mate_pos[i],
                  col           = adjustcolor("#333333", alpha.f = 0.75),
                  lwd           = af_to_lwd(bnd$TUMOR_AF[i]),
                  h.ratio       = 0.65)
    }
  }

  # -------- Title --------
  title(main = sprintf("%s \u2014 Genomic landscape (CNV + high-impact SV)",
                       sample_ttl),
        cex.main = 0.53, font.main = 1, line = 0.2)

  # -------- Legends (all four at the bottom, stacked in two columns) --------
  # CNV category (bottom-left, upper)
  legend("bottomleft", inset = c(0.00, 0.22),
         legend       = names(cnv_pal),
         fill         = cnv_pal, border = NA,
         bty          = "n",
         cex          = 0.53,
         title        = "CNV category",
         title.adj    = 0, title.font = 2,
         y.intersp    = 1.0)

  # BND arc width (bottom-right, upper)
  afs <- c(0.1, 0.5, 1.0)
  legend("bottomright", inset = c(0.00, 0.22),
         legend       = sprintf("AF = %.1f", afs),
         lwd          = af_to_lwd(afs),
         col          = "#333333",
         bty          = "n",
         cex          = 0.53,
         title        = "BND arc width (linear)",
         title.adj    = 0, title.font = 2,
         seg.len      = 2.5,
         y.intersp    = 1.0)

  # SV DEL/DUP ticks (bottom-right, lower)
  legend("bottomright", inset = c(0.00, 0.02),
         legend       = c("DEL tick", "DUP tick"),
         lwd          = 1.2,
         col          = c("#0072B2","#D55E00"),
         bty          = "n",
         cex          = 0.53,
         title        = "SV (no END in input)",
         title.adj    = 0, title.font = 2,
         y.intersp    = 1.0)

  # Notes + gene-label legend (bottom-left, lower)
  legend("bottomleft", inset = c(0.00, 0.02),
         legend       = c("dashed = CNt = 2 (diploid ref)",
                          sprintf("CNt axis capped at %d (95th pctile)", cnt_cap),
                          sprintf("labels: %d CNV-hit + %d SV-hit genes",
                                  nrow(outer_bed), nrow(inner_bed)),
                          "bold = Myc / Trp53",
                          "label colour = CNV palette; SV = green"),
         lty          = c(2, NA, NA, NA, NA),
         lwd          = c(0.6, NA, NA, NA, NA),
         col          = c("grey40", NA, NA, NA, NA),
         text.col     = c("black", "black", "black", "black", "black"),
         bty          = "n",
         cex          = 0.53,
         title        = "Notes",
         title.adj    = 0, title.font = 2,
         y.intersp    = 1.0)
}

# --------------------------- 8. Render ---------------------------------------
out_pdf <- file.path(out_dir, sprintf("%s_circos.pdf", sample_ttl))
draw_circos(out_pdf, "pdf", width_mm = 60)
message(sprintf("Wrote %s", out_pdf))
