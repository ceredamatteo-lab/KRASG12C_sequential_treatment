## ***************************
##
## Script name: KRASi_models_DESeq.R
## Purpose of script:
## Author: mc
## Date Created: 2026-02-05
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



# setting -----
source('Script/utils//environment.R')
source("Script/utils/config_signatures.R")
source("Script/utils/utils_RNA.R")






theme_big2 = function(base_size = 12, base_family = "sans"){
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

fc.gene.counts_path = paste0(DATA_DIR,'Results/Expression')
treat_cntr_patten = "02.FC"
treat_cntr_label = "NAIVE"

treat_case_patterns = c("12.FC|13.FC","22.FC")
treat_case_labels = c("AMG510","MRTX")

cancer_cntr_patten = "h2301"
cancer_cntr_label = "LUAD"

cancer_case_patterns = "cal01"
cancer_case_labels = "LUSC"


f = list.files(path =  fc.gene.counts_path,full.names = T)
fname = list.files(path =  fc.gene.counts_path, full.names = F)

treatment = rep(treat_cntr_label, length(fname))
for (i in seq_along(treat_case_patterns)) { treatment[grep(treat_case_patterns[i], fname)] = treat_case_labels[i]}

cancertype = rep(cancer_cntr_label, length(fname))
for (i in seq_along(cancer_case_patterns)) {cancertype[grep(cancer_case_patterns[i], fname)] = cancer_case_labels[i]}

meta = data.frame(
  barcode = gsub(".FC.gene.counts","",fname)
  , treatment = relevel(factor(treatment), ref = treat_cntr_label)
  , cancertype = cancertype
)

meta$sample = paste0(substr(meta$barcode, 6,10), '-',  substr(meta$barcode, 19, 21) )
rownames(meta) = meta$barcode
meta$cell_line = ifelse(meta$cancertype == 'LUAD', 'H23', 'CALU')

saveRDS(meta, 'Rdata/H23_CALU_metadata.rds')



# get counts -----
counts = get_FC_counts(f, fname)
rownames(counts) = paste0(counts$Geneid,"_rowId_",1:nrow(counts))
counts = subset(counts, gene_type == "protein_coding")
counts_col_ids = grep(paste(fname,collapse = "|"), colnames(counts))
colnames(counts)[counts_col_ids] = gsub(".FC.gene.counts","",colnames(counts)[counts_col_ids])

# identify sample columns (all except annotation columns)
sample_cols <- !(colnames(counts) %in% c("Geneid", "Length", "gene_name", "gene_type"))

# filter genes: keep genes expressed >5 counts in at least 9 samples
keep <- rowSums(counts[, sample_cols] > 5) >= 9

# subset
counts <- counts[keep, ]
dim(counts)

saveRDS(counts, 'Rdata/H23_CALU_counts.rds')
counts = readRDS('Rdata/H23_CALU_counts.rds')


# DESeq2 ---------
model_design = formula(~ cancertype + treatment + cancertype:treatment)

dds = differential_expression_analysis( count_data  = counts[,c(counts_col_ids)]
                                        ,sample_data = sample_data
                                        ,model_design = model_design
                                        ,cores = 1
                                        ,quiet = T) 

saveRDS(dds, 'Rdata/H23_CALU_dds.rds')



# PCA -----
vsd = vst(dds, blind=FALSE)
dim(vsd)
DESeq2::plotPCA(vsd, intgroup = c("treatment", "cancertype")) + geom_text_repel(aes(label = colnames(vsd)), size = 3) + theme_bw()





# DEG -----

resultsNames(dds)


FC_th = 0.5
fdr_th = 0.01

deg.list = list()
deg.list$H23__AMG510_vs_NAIVE = results(dds, name = 'treatment_AMG510_vs_NAIVE')
deg.list$H23__MRTX_vs_NAIVE = results(dds, name = 'treatment_MRTX_vs_NAIVE')
deg.list$CALU__AMG510_vs_NAIVE = results(dds, contrast = list(c("treatment_AMG510_vs_NAIVE", "cancertypeLUSC.treatmentAMG510")))
deg.list$CALU__MRTX_vs_NAIVE = results(dds, contrast = list(c("treatment_MRTX_vs_NAIVE", "cancertypeLUSC.treatmentMRTX")))

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


deg$isRAS84              = ifelse(deg$gene_name %in% RAS84$Lambrechts, 'RAS84', 'rest')
deg$isMAPK_ERK_signature = ifelse(deg$gene_name %in% mapk_erk_signature$gene_symbol, 'MAPK_ERK_signature', 'rest')
deg$isMAPK_signature     = ifelse(deg$gene_name %in% mapk_signature$gene_symbol, 'MAPK_signature', 'rest')
deg$isERK_signature      = ifelse(deg$gene_name %in% erk_signature$gene_symbol, 'ERK_signature', 'rest')

table(deg$contrast, deg$de, sign(deg$log2FoldChange))

saveRDS(deg, 'Rdata/H23_CALU_deg.rds')
deg = readRDS('Rdata/H23_CALU_deg.rds')


supp = subset(deg, padj <= 0.01 & abs(log2FoldChange) > 0.5)
supp = supp[,c(1:9, 13)]
write.csv(supp, 'Tables/Supplementary__H23_CALU_DEG.csv', quote = F, row.names = F)

table(supp$contrast, sign(supp$log2FoldChange))

