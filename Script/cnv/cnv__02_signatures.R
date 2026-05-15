## ***************************
##
## Script name: CNV__03_signatures.R
## Purpose of script:
## Author: mc
## Date Created: 2026-04-17
## R version: R version 4.4.2 (2024-10-31)
##
## ***************************
##
## Notes:
##   
## ***************************



# setting -------
source("Scripts/config/environment.R")
setwd(DATA_DIR)

source("Scripts/config/config_signatures.R")

library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(biomaRt)
library(ComplexHeatmap)

pal_treatment = c('AMG510' = '#511d5a', 'MRTX' = '#d2b693', 'NAIVE' = '#636d9b'
                  , '010' = '#511d5a', '020' = '#d2b693', '000' = '#636d9b'
                  , 'soto' = '#511d5a', 'ada' = '#d2b693', 'naive' = '#636d9b')


# data ---------

ncg = list(cgc        = read.delim2(paste0(GIT_LOCAL_DIR, 'Pipelines/data/NCG7/NCG_cancerdrivers_annotation_supporting_evidence.tsv')) 
           , healthy  = read.delim2(paste0(GIT_LOCAL_DIR,'Pipelines/data/NCG7/NCG_healthydrivers_annotation_supporting_evidence.tsv')) )

cnv.list = list('h23'   = readRDS('Rdata//CNV_sequenza_H23_PC_gene_level_annotated_depmap.rds')
                , 'cal' = readRDS('Rdata//CNV_sequenza_CALU_PC_gene_level_annotated_depmap.rds')
                , 'mep' = subset(readRDS('Rdata/CNV_sequenza_PC_gene_level.rds'), ! model %in% c('h2301', 'cal01')))

cnv.list$mep$gene_id_nv = sapply(strsplit(cnv.list$mep$gene_id, '\\.'), '[[', 1)
colnames(cnv.list$mep)[which(colnames(cnv.list$mep) %in% c('gene_id', 'gene_id_nv', 'gene_name'))] = paste0('mouse_', colnames(cnv.list$mep)[which(colnames(cnv.list$mep) %in% c('gene_id', 'gene_id_nv', 'gene_name'))])

mouse_genes = unique(cnv.list$mep$mouse_gene_id_nv)
# mouse = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://useast.ensembl.org")
mouse = readRDS('Rdata/mouse__biomart.rds')
human_map = getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), filters = "ensembl_gene_id", values = mouse_genes, mart = mouse)
human_map = subset(human_map, hsapiens_homolog_associated_gene_name != "" )

cnv.list$mep$gene_name = human_map$hsapiens_homolog_associated_gene_name[match(cnv.list$mep$mouse_gene_id_nv, human_map$ensembl_gene_id)]

# Add signatures --------- 

cnv.list[c('h23', 'cal')] = lapply(cnv.list[c('h23', 'cal')], function(x){
  hc.cnv                      = x
  hc.cnv$isKRAS_upregulated   = ifelse(hc.cnv$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_up-regulated' )$gene_name, 'KRAS_up-regulated', 'rest')
  hc.cnv$isKRAS_downregulated = ifelse(hc.cnv$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_down-regulated' )$gene_name, 'KRAS_down-regulated', 'rest')
  hc.cnv$isRAS84              = ifelse(hc.cnv$gene_id_nv %in% RAS84$Uppsala_feature_ids, 'RAS84', 'rest')
  hc.cnv$isMAPK_ERK_signature = ifelse(hc.cnv$gene_id_nv %in% mapk_erk_signature$ensembl_gene, 'MAPK_ERK_signature', 'rest')
  hc.cnv$isMAPK_signature     = ifelse(hc.cnv$gene_id_nv %in% mapk_signature$ensembl_gene, 'MAPK_signature', 'rest')
  hc.cnv$isERK_signature      = ifelse(hc.cnv$gene_id_nv %in% erk_signature$ensembl_gene, 'ERK_signature', 'rest')
  hc.cnv                      = isNCGdriver(hc.cnv, ncg = ncg, gene_identifier = 'symbol', primarySite = 'lung' )
  hc.cnv                      = isGC(hc.cnv, gc_file = paste0(GIT_LOCAL_DIR,"/Pipelines/GeneCategories/Rdata/gene_categories.20251120.rds"), gene_identifier = 'gene_id_nv')
  hc.cnv

})


cnv.list[c('mep')] = lapply(cnv.list[c('mep')], function(x){
  m.cnv= x
  m.cnv$isKRAS_upregulated   = ifelse(m.cnv$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_up-regulated' )$gene_name, 'KRAS_up-regulated', 'rest')
  m.cnv$isKRAS_downregulated = ifelse(m.cnv$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_down-regulated' )$gene_name, 'KRAS_down-regulated', 'rest')
  m.cnv$isRAS84              = ifelse(m.cnv$gene_name %in% RAS84$Lambrechts, 'RAS84', 'rest')
  m.cnv$isMAPK_ERK_signature = ifelse(m.cnv$gene_name %in% mapk_erk_signature$gene_symbol, 'MAPK_ERK_signature', 'rest')
  m.cnv$isMAPK_signature     = ifelse(m.cnv$gene_name %in% mapk_signature$gene_symbol, 'MAPK_signature', 'rest')
  m.cnv$isERK_signature      = ifelse(m.cnv$gene_name %in% erk_signature$gene_symbol, 'ERK_signature', 'rest')
  m.cnv                      = isNCGdriver(m.cnv, ncg = ncg, gene_identifier = 'symbol', primarySite = 'lung' )
  m.cnv                      = isGC(m.cnv, gc_file = paste0(GIT_LOCAL_DIR,"/Pipelines/GeneCategories/Rdata/gene_categories.20251120.rds"), gene_identifier = 'gene_name')
  m.cnv
})

saveRDS(cnv.list$h23, 'Rdata/CNV_sequenza_H23_PC_gene_level_annotated_depmap_signatures.rds')
saveRDS(cnv.list$cal, 'Rdata/CNV_sequenza_CALU_PC_gene_level_annotated_depmap_signatures.rds')
saveRDS(cnv.list$mep, 'Rdata/CNV_sequenza_MEP_PC_gene_level_annotated_depmap_signatures.rds')







