#!/usr/bin/env Rscript
# =============================================================================
# run_ORA_H23_CALU.R
# Over-representation analysis on the H23 / CALU DEG set, split by contrast
# and direction (up/down) and run against every gene-set collection in
# go_list__20260214.rds.
#
# DEG filter: |log2FoldChange| > 0.5 AND padj < 0.01 (NA padj dropped).
# Identifier:  gene_id_nv (unversioned ENSG) — matches go$<set>$ensembl_gene.
#
# Inputs   : input/H23_CALU_deg.rds
#            input/go_list__20260214.rds
# Outputs  : ora_H23.rds, ora_CALU.rds at the project root (full per-row tables)
#
# Sources Script/utils_ORA.R which auto-installs+loads the package set
# (msigdbr, clusterProfiler, dplyr, forcats, ...).
# =============================================================================

.script_dir <- local({
  cargs <- commandArgs(trailingOnly = FALSE)
  f <- cargs[grepl("^--file=", cargs)]
  if (length(f)) normalizePath(dirname(sub("^--file=", "", f[1])),
                               mustWork = FALSE) else getwd()
})
.project_root <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
INPUT_DIR <- file.path(.project_root, "input")

source(file.path(.script_dir, "..", "utils", "utils_ORA.R"))

suppressPackageStartupMessages({
  library(DOSE)
  library(stringi)
})

deg <- readRDS(file.path(INPUT_DIR, "H23_CALU_deg.rds"))
go  <- readRDS(file.path(INPUT_DIR, "go_list__20260214.rds"))

message(sprintf("DEG rows: %d   contrasts: %s",
                nrow(deg),
                paste(sort(unique(deg$contrast)), collapse = ", ")))
message(sprintf("Gene-set collections: %s",
                paste(names(go), collapse = ", ")))

sig <- deg %>%
  filter(!is.na(padj),
         padj < 0.01,
         abs(log2FoldChange) > 0.5)

message(sprintf("\nAfter filter (|log2FC|>0.5 & padj<0.01): %d rows", nrow(sig)))

ora.list <- list()
for (con in sort(unique(sig$contrast))) {
  sub <- sig[sig$contrast == con, , drop = FALSE]
  up   <- unique(sub$gene_id_nv[sub$log2FoldChange >  0])
  down <- unique(sub$gene_id_nv[sub$log2FoldChange <  0])
  all  <- unique(sub$gene_id_nv)
  ora.list[[con]] <- list(up = up, down = down, all = all)
  message(sprintf("  %-28s up = %5d   down = %5d   all = %5d",
                  con, length(up), length(down), length(all)))
}

ora <- lapply(ora.list, function(x){
  lapply(x, function(y){
    ora <- ORA(y, gene_set_list = go)
    ora <- add_DOSE_measure_to_ORA(ora)
    ora$Description <- gsub('KEGG_|HALLMARK_|REACTOME_|BIOCARTA_', '', ora$Description)
    ora$Description <- gsub('_', ' ', ora$Description)
    ora$Description <- stringi::stri_trans_totitle(ora$Description)
    ora
  }) %>% bind_rows(.id = 'de.status')
})

ora.h23  <- bind_rows(ora[grepl('H23',  names(ora))], .id = 'contrast')
ora.calu <- bind_rows(ora[grepl('CALU', names(ora))], .id = 'contrast')

out_h23  <- file.path(.project_root, "ora_H23.rds")
out_calu <- file.path(.project_root, "ora_CALU.rds")
saveRDS(ora.h23,  out_h23)
saveRDS(ora.calu, out_calu)

message(sprintf("\nWrote %s (%d rows)",  out_h23,  nrow(ora.h23)))
message(sprintf("Wrote %s (%d rows)",     out_calu, nrow(ora.calu)))

cat("\n=== H23 ORA: top hits (padj<0.05) by contrast/de.status/ontology ===\n")
print(ora.h23  %>% filter(p.adjust < 0.05) %>% count(contrast, de.status, ontology))
cat("\n=== CALU ORA: top hits (padj<0.05) by contrast/de.status/ontology ===\n")
print(ora.calu %>% filter(p.adjust < 0.05) %>% count(contrast, de.status, ontology))
