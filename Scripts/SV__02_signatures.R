## ***************************
##
## Script name: SV__01_processing.R
## Purpose of script:
## Author: mc
## Date Created: 2026-04-23
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
library(biomaRt)


# data --------

sv.list = list('h2301'   = readRDS('Rdata/caks-h2301_vs_mgpf-af001-10-050-000.filtered__overlap_manta_tiddit_no_filters.rds')
               , 'cal01' = readRDS('Rdata/caks-cal01_vs_mgpf-ib001-10-050-000.filtered__overlap_manta_tiddit_no_filters.rds')
               , 'mep01' = readRDS('Rdata/caks-mep01_vs_mef_BL6Nx129s1.filtered__overlap_manta_tiddit_no_filters.rds' )
)


sv = lapply(sv.list, function(x){subset(x, BIOTYPE == 'protein_coding' & FILTER == 'PASS' & IMPACT == 'HIGH')})
sapply(sv, dim)
# h2301 cal01 mep01
# [1,]   710   303   624
# [2,]   120   120   120


mouse_genes = unique(sv$mep01$Gene)
# mouse = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://useast.ensembl.org")
# saveRDS(mouse, 'Rdata/mouse__biomart.rds')
mouse = readRDS('Rdata/mouse__biomart.rds')
human_map = getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name", 'hsapiens_homolog_ensembl_gene'), filters = "ensembl_gene_id", values = mouse_genes, mart = mouse)
human_map = subset(human_map, hsapiens_homolog_ensembl_gene != "" )

# For mouse, gene_name and gene_id_nv will correspond to the human orthologs (used to annotate signatures)
sv$mep01$gene_name = human_map$hsapiens_homolog_associated_gene_name[match(sv$mep01$Gene, human_map$ensembl_gene_id)]
sv$mep01$gene_id_nv = human_map$hsapiens_homolog_ensembl_gene[match(sv$mep01$Gene, human_map$ensembl_gene_id)]

sv$h2301$gene_id_nv = sv$h2301$Gene
sv$cal01$gene_id_nv = sv$cal01$Gene
sv$h2301$gene_name = sv$h2301$SYMBOL
sv$cal01$gene_name = sv$cal01$SYMBOL





# Computing AF + annotating ------

sv.anno = lapply(sv, function(x){
  x = x %>%
    tidyr::separate(NORMAL_PR, into = c("normal_PR_ref", "normal_PR_alt"), sep = ",", convert = TRUE) %>%
    tidyr::separate(TUMOR_PR,  into = c("tumor_PR_ref",  "tumor_PR_alt"),  sep = ",", convert = TRUE) %>%
    tidyr::separate(NORMAL_SR, into = c("normal_SR_ref", "normal_SR_alt"), sep = ",", convert = TRUE, fill = "right") %>%
    tidyr::separate(TUMOR_SR,  into = c("tumor_SR_ref",  "tumor_SR_alt"),  sep = ",", convert = TRUE, fill = "right") %>%
    
    mutate(across(c(normal_SR_ref, normal_SR_alt, tumor_SR_ref, tumor_SR_alt), ~ replace_na(., 0))) %>%
    mutate(TUMOR_AF  = (tumor_PR_alt  + tumor_SR_alt)  / (tumor_PR_ref  + tumor_PR_alt  + tumor_SR_ref  + tumor_SR_alt) )
  
  x$isKRAS_upregulated   = ifelse(x$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_up-regulated' )$gene_name, 'KRAS_up-regulated', 'rest')
  x$isKRAS_downregulated = ifelse(x$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_down-regulated' )$gene_name, 'KRAS_down-regulated', 'rest')
  x$isRAS84              = ifelse(x$gene_id_nv %in% RAS84$Uppsala_feature_ids, 'RAS84', 'rest')
  x$isMAPK_ERK_signature = ifelse(x$gene_id_nv %in% mapk_erk_signature$ensembl_gene, 'MAPK_ERK_signature', 'rest')
  x$isMAPK_signature     = ifelse(x$gene_id_nv %in% mapk_signature$ensembl_gene, 'MAPK_signature', 'rest')
  x$isERK_signature      = ifelse(x$gene_id_nv %in% erk_signature$ensembl_gene, 'ERK_signature', 'rest')
  x                      = isNCGdriver(x, ncg = ncg, gene_identifier = 'symbol', primarySite = 'lung' )
  x                      = isGC(x, gc_file = paste0(GIT_LOCAL_DIR,"/Pipelines/GeneCategories/Rdata/gene_categories.20251120.rds"), gene_identifier = 'gene_id_nv')
  x$barcode              = x$sample
  
  return(x)
  
  
})


saveRDS(sv.anno, 'Rdata/SV__selected_annotated.rds')



sv.anno.sel = lapply(sv.anno, function(x){
  x = x[,c('CHROM','POS', 'REF', 'ALT','ID'
           , 'SVTYPE', 'Consequence', 'IMPACT'
           , 'gene_name', 'gene_id_nv'
           , 'VARIANT_CLASS', 'TUMOR_AF'
           , 'barcode'
           , colnames(x)[grepl('SYMBOL$|Gene$', colnames(x))]
           , colnames(x)[grepl('^is', colnames(x))]
           
  )]
  x$Consequence = gsub('missense_variant&', '', x$Consequence)
  
  return(x)
})

saveRDS(sv.anno, 'Rdata/SV__selected_annotated_vaf.rds')










