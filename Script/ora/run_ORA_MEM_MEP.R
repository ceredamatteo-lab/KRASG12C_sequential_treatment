#!/usr/bin/env Rscript
# =============================================================================
# run_ORA_MEM_MEP.R
# Mouse counterpart of run_ORA_H23_CALU.R for the MEF lines MEM and MEP.
#
# Over-representation analysis on the MEM / MEP DEG set, split by contrast
# and direction (up / down / all). Mouse gene-set lists are fetched from
# msigdbr (species = "Mus musculus") at runtime — no precomputed go_list
# is read since input/go_list__20260214.rds is human-only.
#
# DEG filter: |log2FoldChange| > 0.5 AND padj < 0.01 (NA padj dropped).
# Identifier:  unversioned ENSMUSG (deg$gene_id_nv).
#
# Inputs   : input/MEM_MEP_deg.rds
# Outputs  : ora_MEM.rds, ora_MEP.rds at the project root
#
# Sources Script/utils_ORA.R for ORA() and add_DOSE_measure_to_ORA().
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
  library(msigdbr)
})

# ---- mouse gene-set list (built on the fly) ---------------------------------
# Same five ontology keys used downstream by plot_ORA_MEM_MEP.R: kegg,
# hallmarks, go-bp, reactome, wiki. Each value is a 2-column tibble of
# (gs_name cleaned, ENSMUSG) — the format ORA() / clusterProfiler::enricher
# expects via TERM2GENE.
build_mouse_go <- function() {
  pull_set <- function(coll, sub = NULL) {
    args <- list(species = "Mus musculus", collection = coll)
    if (!is.null(sub)) args$subcollection <- sub
    do.call(msigdbr, args) %>%
      transmute(gs_name = gs_name, ensembl_gene = ensembl_gene) %>%
      distinct()
  }
  list(
    kegg      = pull_set("C2", "CP:KEGG_LEGACY") %>%
                  mutate(gs_name = sub("^KEGG_", "", gs_name)),
    hallmarks = pull_set("H") %>%
                  mutate(gs_name = sub("^HALLMARK_", "", gs_name)),
    `go-bp`   = pull_set("C5", "GO:BP") %>%
                  mutate(gs_name = sub("^GOBP_", "", gs_name)),
    reactome  = pull_set("C2", "CP:REACTOME") %>%
                  mutate(gs_name = sub("^REACTOME_", "", gs_name)),
    wiki      = pull_set("C2", "CP:WIKIPATHWAYS") %>%
                  mutate(gs_name = sub("^WP_", "", gs_name))
  )
}

deg <- readRDS(file.path(INPUT_DIR, "MEM_MEP_deg.rds"))
go  <- build_mouse_go()

message(sprintf("DEG rows: %d   contrasts: %s",
                nrow(deg),
                paste(sort(unique(deg$contrast)), collapse = ", ")))
message(sprintf("Gene-set collections (mouse, msigdbr): %s",
                paste(names(go), collapse = ", ")))
message("Per-collection set counts:")
for (nm in names(go)) {
  message(sprintf("  %-10s  %5d sets   %5d (set, gene) pairs",
                  nm,
                  dplyr::n_distinct(go[[nm]]$gs_name),
                  nrow(go[[nm]])))
}

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

ora <- lapply(ora.list, function(x) {
  lapply(x, function(y) {
    ora <- ORA(y, gene_set_list = go)
    ora <- add_DOSE_measure_to_ORA(ora)
    ora$Description <- gsub('KEGG_|HALLMARK_|REACTOME_|BIOCARTA_|WP_|GOBP_',
                            '', ora$Description)
    ora$Description <- gsub('_', ' ', ora$Description)
    ora$Description <- stringi::stri_trans_totitle(ora$Description)
    ora
  }) %>% bind_rows(.id = 'de.status')
})

ora.mem <- bind_rows(ora[grepl('MEM', names(ora))], .id = 'contrast')
ora.mep <- bind_rows(ora[grepl('MEP', names(ora))], .id = 'contrast')

out_mem <- file.path(.project_root, "ora_MEM.rds")
out_mep <- file.path(.project_root, "ora_MEP.rds")
saveRDS(ora.mem, out_mem)
saveRDS(ora.mep, out_mep)

message(sprintf("\nWrote %s (%d rows)", out_mem, nrow(ora.mem)))
message(sprintf("Wrote %s (%d rows)",     out_mep, nrow(ora.mep)))

cat("\n=== MEM ORA: hits at p.adjust<0.05 by contrast/de.status/ontology ===\n")
print(ora.mem %>% filter(p.adjust < 0.05) %>%
        count(contrast, de.status, ontology))
cat("\n=== MEP ORA: hits at p.adjust<0.05 by contrast/de.status/ontology ===\n")
print(ora.mep %>% filter(p.adjust < 0.05) %>%
        count(contrast, de.status, ontology))
