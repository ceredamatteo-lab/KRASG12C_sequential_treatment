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




source('Scripts/environment.R')
source(paste0(GIT_LOCAL_DIR, "Pipelines/Utils/Utils_RNA.R"))
source(paste0(GIT_LOCAL_DIR, "Pipelines/Utils/utils_RNA_DESeq2.R"))

source(paste0(GIT_LOCAL_DIR, "Pipelines/Utils/Utils_ORA.R"))
source(paste0(GIT_LOCAL_DIR, "Pipelines/Utils/Utils_statistical_testing.R"))



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


supp = subset(deg, padj <= 0.01 & abs(log2FoldChange) > 0.5)
supp = supp[,c(1:9, 14)]
write.csv(supp, 'Tables/Supplementary__MEM_MEP_DEG.csv', quote = F, row.names = F)





# ORA --------

go.mouse = readRDS('Rdata/go_list_mouse_kegg_hallmarks.rds')
ora.list = dlply(supp, .(contrast)) %>% lapply(., function(x){
  deg = list('deg' = x$gene_id_nv
             
  )
})


sub.go = list('kegg' = go.mouse$kegg, 'hallmarks' = go.mouse$hallmarks)

ora = lapply(ora.list, function(x){
  lapply(x, function(y){ 
    ora = ORA(y, gene_set_list = sub.go)
    ora = add_DOSE_measure_to_ORA(ora)
    # ora = retrieving_gene_name_ora(ora, gene_info = gene_info_table)
    ora$Description = gsub('KEGG_|HALLMARK_|REACTOME_|BIOCARTA_', '', ora$Description)
    ora$Description = gsub('_', ' ', ora$Description)
    ora$Description = stringi::stri_trans_totitle(ora$Description)
    return(ora) }) %>% bind_rows(., .id = 'de.status')
}) %>% bind_rows(., .id = 'contrast')

ora$condition = ifelse(grepl('AMG510', ora$contrast), 'sotorasib', 'adagrasib')
ora$cell = ifelse(grepl('MEM', ora$contrast), 'MEM', 'MEP')
ora$geneRatio.mod = as.numeric(sapply(strsplit(as.character(ora$GeneRatio), '\\/'), '[[', 1))
ora$cell = factor(ora$cell, levels = c('MEM', 'MEP'))

ora %>%  
  filter(grepl('kegg$|hallmarks$', ontology), p.adjust < 0.1 ) %>% 
  dplyr::group_by(ontology, condition, cell) %>%
  slice_head(n = 10) %>%
  ggplot(., aes(richFactor, Description, fill = condition, size = geneRatio.mod, shape = condition )) +
  geom_point(shape = 21) +
  scale_color_manual(values = c('black', 'transparent')) +
  # scale_fill_viridis_c(option = 'D', guide = guide_colorbar(reverse = T, draw.llim = T), direction = -1) +
  facet_grid(ontology ~ cell, scales = 'free', space = 'free') +
  theme_minimal() +
  xlab("") +
  ylab(NULL) +
  theme(text = element_text(size = 12),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90))



