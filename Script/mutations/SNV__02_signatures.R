## ***************************
##
## Script name: SNV__01_signatures.R
## Purpose of script:
## Author: 
## Date Created: 2026-04-21
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

# data ----
snv.list = list('h2301' = readRDS('Rdata/caks-h2301_vs_mgpf-af001-10-050-000__filtered_overlap_mutect2_strelka.rds')
                , 'cal01' = readRDS('Rdata/caks-cal01_vs_mgpf-ib001-10-050-000__filtered_overlap_mutect2_strelka.rds')
                , 'mep01' = readRDS('Rdata/caks-mep01_vs_mef_BL6Nx129s1__filtered_overlap_mutect2_strelka.rds' ))

sapply(snv.list, dim)
# h2301 cal01   mep01
# [1,] 275829 48168 3091847

snv = lapply(snv.list, function(x){subset(x, BIOTYPE == 'protein_coding' )})
snv = lapply(snv, function(x){x$condition = sapply(strsplit(x$sample, '\\-'), '[[', 5); x })

sapply(snv, dim)
# h2301 cal01   mep01
# [1,] 84320 13131 1156464


# We filter mutations for having IMPACT = high or moderate and NOT missense variant (that in the annovar selection are identified as the synonymous snvs)
snv.impact = lapply(snv, function(x){subset(x, Consequence != 'missense_variant' & IMPACT %in% c('HIGH', 'MODERATE') )})
sapply(snv.impact, dim)
# h2301 cal01 mep01
# [1,]   150    11   493
saveRDS(snv.impact, 'Rdata/SNV__overlapStrelkaMutect_VEPHighModerate.rds')


# mouse genes ----
mouse_genes = unique(snv.impact$mep01$Gene)

mouse = readRDS('Rdata/mouse__biomart.rds')
human_map = getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name", 'hsapiens_homolog_ensembl_gene'), filters = "ensembl_gene_id", values = mouse_genes, mart = mouse)
human_map = subset(human_map, hsapiens_homolog_ensembl_gene != "" )


# For mouse, gene_name and gene_id_nv will correspond to the human orthologs (used to annotate signatures)
snv.impact$mep01$gene_name = human_map$hsapiens_homolog_associated_gene_name[match(snv.impact$mep01$Gene, human_map$ensembl_gene_id)]
snv.impact$mep01$gene_id_nv = human_map$hsapiens_homolog_ensembl_gene[match(snv.impact$mep01$Gene, human_map$ensembl_gene_id)]

snv.impact$h2301$gene_id_nv = snv.impact$h2301$Gene
snv.impact$cal01$gene_id_nv = snv.impact$cal01$Gene
snv.impact$h2301$gene_name = snv.impact$h2301$SYMBOL
snv.impact$cal01$gene_name = snv.impact$cal01$SYMBOL



# Add signatures --------- 

snv.impact = lapply(snv.impact, function(x){
  anno                      = x
  anno$isKRAS_upregulated   = ifelse(anno$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_up-regulated' )$gene_name, 'KRAS_up-regulated', 'rest')
  anno$isKRAS_downregulated = ifelse(anno$gene_name %in% subset(kras_signature, kras_regulation == 'KRAS_down-regulated' )$gene_name, 'KRAS_down-regulated', 'rest')
  anno$isRAS84              = ifelse(anno$gene_id_nv %in% RAS84$Uppsala_feature_ids, 'RAS84', 'rest')
  anno$isMAPK_ERK_signature = ifelse(anno$gene_id_nv %in% mapk_erk_signature$ensembl_gene, 'MAPK_ERK_signature', 'rest')
  anno$isMAPK_signature     = ifelse(anno$gene_id_nv %in% mapk_signature$ensembl_gene, 'MAPK_signature', 'rest')
  anno$isERK_signature      = ifelse(anno$gene_id_nv %in% erk_signature$ensembl_gene, 'ERK_signature', 'rest')
  anno                      = isNCGdriver(anno, ncg = ncg, gene_identifier = 'symbol', primarySite = 'lung' )
  anno                      = isGC(anno, gc_file = paste0(GIT_LOCAL_DIR,"/Pipelines/GeneCategories/Rdata/gene_categories.20251120.rds"), gene_identifier = 'gene_id_nv')
  anno$barcode              = anno$sample
  return(anno)
})

saveRDS(snv.impact, 'Rdata/SNV__selected_annotated.rds')


vaf.filt = lapply(snv.impact, function(x){
  x$ID = paste0(x$ID, '_', x$gene_name)
  x = x[,c('CHROM', 'POS', 'ID', 'TUMOR_AF', 'TUMOR_DP', 'sample', 'condition'
           , 'Consequence', 'IMPACT'
           , 'gene_name', 'gene_id_nv'
           , 'EXON', 'cDNA_position', 'CDS_position', 'Protein_position', 'Codons'
           , 'Existing_variation', 'VARIANT_CLASS'
           , 'barcode'
           , colnames(x)[grepl('SYMBOL$|Gene$', colnames(x))]
           , colnames(x)[grepl('^is', colnames(x))]
           
           #, 'annovar_ExonicFunc.refGen'
  )]
  x$Consequence = gsub('missense_variant&', '', x$Consequence)
  
  return(x)
})


saveRDS(vaf.filt, 'Rdata/SNV__selected_annotated_vaf.rds')


