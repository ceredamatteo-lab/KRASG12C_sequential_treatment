## ***************************
##
## Script name: collecting_model_cnv_snv_sv.R
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


library(plyr)
library(dplyr)


# data ----------
cnv.list = list('h2301'   = readRDS('Rdata//CNV_sequenza_H23_PC_gene_level_annotated_depmap_signatures.rds')
                , 'cal01' = readRDS('Rdata//CNV_sequenza_CALU_PC_gene_level_annotated_depmap_signatures.rds')
                , 'mep01' = readRDS('Rdata//CNV_sequenza_MEP_PC_gene_level_annotated_depmap_signatures.rds'))


sv.anno = readRDS('Rdata/SV__selected_annotated.rds')
vaf.filt = readRDS('Rdata/SNV__selected_annotated_vaf.rds')


# H23 
h23 = list('CNV'   = cnv.list$h2301
           , 'SNV' = vaf.filt$h2301
           , 'SV'  = sv.anno$h2301)
saveRDS(h23, 'Rdata/h23_cnv_snv_sv.rds')



# CALU
calu = list('CNV'   = cnv.list$cal01
            , 'SNV' = vaf.filt$cal01
            , 'SV'  = sv.anno$cal01)
saveRDS(calu, 'Rdata/calu_cnv_snv_sv.rds')


# MEP
mep = list('CNV'   = cnv.list$mep01
           , 'SNV' = vaf.filt$mep01
           , 'SV'  = sv.anno$mep01)
saveRDS(mep, 'Rdata/mep_cnv_snv_sv.rds')


