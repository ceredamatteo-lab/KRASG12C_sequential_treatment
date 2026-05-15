## ***************************
##
## Script name: KRASi_mouse_model_DESeq.R 
## Purpose of script:
## Author: 
## Date Created: 2026-03-11
## R version: R version 4.4.2 (2024-10-31)
##
## ***************************
##
## Notes:
##   
## ***************************



library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(ggrepel)
library(ComplexHeatmap)
library(biomaRt)



# setting -----
source('Script/utils//environment.R')
source("Script/utils/config_signatures.R")
source("Script/utils/utils_RNA.R")




theme_big2 <- function(base_size = 12, base_family = "sans"){
  theme(
    legend.text = element_text( size = 10),
    plot.title = element_text(size = 18),
    axis.title = element_text(size=11),
    legend.title = element_text( size = 15),
    axis.text = element_text(size=10),#,legend.title=element_blank()
    strip.text = element_text(size = 12)
  )
}



# palette 
pal_cancertype = c('LUSC' = '#d23694', 'LUAD' = '#0a757d')
pal_cell_line = c('CALU' = '#d23694', 'H23' = '#0a757d')


pal_treatment = c('AMG510' = '#511d5a', 'MRTX' = '#d2b693', 'NAIVE' = '#636d9b')




# data ------

path_to_FC_counts = "Results/Expression"
pattern = paste(c('mep', 'mem'), collapse = "|")

f     = list.files(path =  path_to_FC_counts, pattern = pattern, recursive =F, full.names = T)
fname = list.files(path =  path_to_FC_counts, pattern = pattern, recursive =F, full.names = F)

counts = get_FC_counts(f, fname)
colnames(counts)[grepl('caks', colnames(counts))] = gsub(".FC.gene.counts", "", colnames(counts)[grepl('caks', colnames(counts))])
counts$gene_id_nv = sapply(strsplit(counts$Geneid, '\\.'), '[[', 1)
rownames(counts) = counts$Geneid
counts = subset(counts, gene_type == 'protein_coding')

saveRDS(counts, 'Rdata/MEM_MEP_counts.rds')

meta = data.frame('barcode' = gsub(".FC.gene.counts", "", colnames(counts)[grepl('caks', colnames(counts))]) )
meta$treatment = ifelse(grepl('00$', meta$barcode), 'NAIVE', ifelse(grepl('10$', meta$barcode), 'AMG510', 'MRTX'))
meta$treatment =  relevel(factor(meta$treatment), ref = 'NAIVE')

meta$condition = ifelse(grepl('mep', meta$barcode), 'mep','mem')
rownames(meta) = meta$barcode


dds = differential_expression_analysis( count_data  = counts[,grepl('mem|mep', colnames(counts))]
                                            ,sample_data = meta
                                            ,model_design = ~ condition + treatment + condition:treatment
                                            ,cores = 1
                                            ,quiet = T) 



# PCA -----

vsd = vst(dds, blind=FALSE)
dim(vsd)
DESeq2::plotPCA(vsd, intgroup = c("treatment")) + geom_text_repel(aes(label = colnames(vsd)), size = 3) + theme_bw()


# DEG ------
FC_th = 0.5
fdr_th = 0.01

resultsNames(dds)

deg.list = list()
deg.list$MEM__AMG510_vs_NAIVE = results(dds, name = 'treatment_AMG510_vs_NAIVE')
deg.list$MEM__MRTX_vs_NAIVE = results(dds, name = 'treatment_MRTX_vs_NAIVE')
deg.list$MEP__AMG510_vs_NAIVE = results(dds, contrast = list(c("treatment_AMG510_vs_NAIVE", "conditionmep.treatmentAMG510")))
deg.list$MEP__MRTX_vs_NAIVE = results(dds, contrast = list(c("treatment_MRTX_vs_NAIVE", "conditionmep.treatmentMRTX")))


deg.list = lapply(deg.list, function(x){
  x$gene_id = sapply(strsplit(rownames(x), '\\_'), '[[', 1)
  x$gene_id_nv = sapply(strsplit(x$gene_id, '\\.'), '[[', 1)
  x$formula = x@elementMetadata[2,2]
  rownames(x) = NULL
  return(as.data.frame(x))
})

deg = bind_rows(deg.list, .id = 'contrast')

deg$ss = ifelse(deg$padj <= fdr_th, TRUE, FALSE )
deg$de = ifelse(abs(deg$log2FoldChange) >= FC_th & deg$ss, TRUE, FALSE )
deg$gene_name = counts$gene_name[match(deg$gene_id, counts$Geneid)]






saveRDS(deg, 'Rdata/MEM_MEP_deg.rds')
write.csv(deg, 'Tables/MEM_MEP_deg.csv', quote = F, row.names = F)


supp = subset(deg, padj < 0.01 & abs(log2FoldChange) > 0.5)
supp = supp[,c(1:9, 14)]
write.csv(supp, 'Tables/Supplementary__MEM_MEP_DEG.csv', quote = F, row.names = F)


deg = readRDS( 'Rdata/MEM_MEP_deg.rds')

