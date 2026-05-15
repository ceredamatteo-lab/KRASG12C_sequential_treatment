## ***************************
##
## Script name: CNV__01_segment2gene_levels.R
## Purpose of script:
## Author: mc
## Date Created: 2026-04-15
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
library(rtracklayer)
library(biomaRt)
library(tibble)
library(tidyr)

saveBED = function(data, file){write.table(data, file, quote = FALSE, col.names  = FALSE, sep = '\t', row.names = FALSE) }




# Load data -----
selected = c(
  'Results/Sequenza/fromVCF//caks-h2301-50-050-000_vs_mgpf-af001-10-050-000_COV5/sequenza_results_1/caks-h2301-50-050-000_segments_logR.txt'
  , 'Results/Sequenza/fromVCF//caks-h2301-50-050-010_vs_mgpf-af001-10-050-000_COV5/sequenza_results_1/caks-h2301-50-050-010_segments_logR.txt'
  , 'Results/Sequenza/fromVCF//caks-h2301-50-050-020_vs_mgpf-af001-10-050-000_COV5/sequenza_results_1/caks-h2301-50-050-020_segments_logR.txt'
  
  , "Results/Sequenza/fromVCF//caks-cal01-51-050-000_vs_mgpf-ib001-10-050-000_COV5/sequenza_results_1/caks-cal01-51-050-000_segments_logR.txt"
  , "Results/Sequenza/fromVCF//caks-cal01-51-050-010_vs_mgpf-ib001-10-050-000_COV5/sequenza_results_1/caks-cal01-51-050-010_segments_logR.txt"
  , "Results/Sequenza/fromVCF//caks-cal01-51-050-020_vs_mgpf-ib001-10-050-000_COV5/sequenza_results_1/caks-cal01-51-050-020_segments_logR.txt"
  
  , "Results/Sequenza/fromVCF//caks-mep01-53-050-000_vs_mef_COV5/sequenza_results_1/caks-mep01-53-050-000_segments_logR.txt"
  , "Results/Sequenza/fromVCF//caks-mep01-53-050-010_vs_mef_COV5/sequenza_results_1/caks-mep01-53-050-010_segments_logR.txt"
  , "Results/Sequenza/fromVCF//caks-mep01-53-050-020_vs_mef_COV5/sequenza_results_1_altsol2/caks-mep01-53-050-020_segments_logR.txt"
  
  
  
)


cnv = lapply(selected, function(x){
  cnv.table = read.delim2(x)
  cnv.table$name = x
  return(cnv.table)
}) %>% bind_rows()

cnv$comparison = sapply(strsplit(cnv$name, '/'), '[[', 5)
cnv$barcode = sapply(strsplit(cnv$comparison, '_'), '[[', 1)
cnv$model = sapply(strsplit(cnv$barcode, '-'), '[[', 2)
cnv$id =  paste0(cnv$chromosome, '_', cnv$start.pos, '_', cnv$end.pos)
cnv$width = cnv$end.pos - cnv$start.pos

head(cnv)

# ggplot(cnv, aes(width, col = model )) + geom_density() + scale_x_log10() + theme_bw()
# ggplot(cnv, aes(width, col = barcode )) + geom_density() + scale_x_log10() + theme_bw()
# ggplot(cnv, aes(barcode, width, fill = model) ) + geom_boxplot(notch = TRUE) + scale_y_log10() + theme_bw()
# ggplot(cnv, aes(barcode, N.BAF, fill = model) ) + geom_boxplot(notch = TRUE) + scale_y_log10() + theme_bw()
# ggplot(cnv, aes(barcode, as.numeric( depth.ratio), fill = model) ) + geom_boxplot(notch = TRUE) + scale_y_log10() + theme_bw()
# ggplot(cnv, aes(barcode, as.numeric( sd.ratio), fill = model) ) + geom_boxplot(notch = TRUE) + scale_y_log10() + theme_bw()
# ggplot(cnv, aes(barcode, as.numeric( CNt), fill = model) ) + geom_boxplot(notch = TRUE) + scale_y_log10() + theme_bw()
# 
# nsgments = ddply(cnv, .(barcode), summarise , n = length(id), mean.width = mean(width))
# nchr = ddply(cnv, .(barcode, chromosome), summarise , n = length(id), mean.width = mean(width))
# ggplot(nchr , aes(barcode, n, fill = chromosome) ) + geom_bar(stat = 'identity', position = 'fill', col = 'black')


# Retrieving gene-level ---------
cnv.list = dlply(cnv, ~ model )

# Generating bed for gene level intersection
cnv.bed = lapply(cnv.list, function(x){
  x$id = paste0(x$chromosome, '_', x$start.pos, '_', x$end.pos)
  model = unique(x$model)
  bed = x[,c('chromosome','start.pos','end.pos', 'id')]
  bed = unique(bed)
  saveBED(bed, paste0('Tables/CNV_sequenza_selected_', model, '.bed') )
  
})

# Retrieve gene annotation bed files
# gtf = import('Tables/gencode.v46.annotation.gtf.gz')
# gtf = subset(gtf , type == 'gene')
# gtf = as.data.frame(gtf)
# 
# saveBED(gtf[,c('seqnames', 'start', 'end', 'gene_id','gene_name', 'gene_type')], 'Tables/gencode.v46.annotation.genes.bed' )
# 
# 
# gtf = import('Tables/gencode.vM36.annotation_add_hKRAS_G12C.gtf')
# gtf = subset(gtf , type == 'gene')
# gtf = as.data.frame(gtf)
# 
# saveBED(gtf[,c('seqnames', 'start', 'end', 'gene_id','gene_name', 'gene_type')], 'Tables/gencode.vM36.annotation_add_hKRAS_G12C.genes.bed' )

# bedtools intersect -a 'Tables/gencode.v46.annotation.genes.bed' -b 'Tables/CNV_sequenza_selected_cal01.bed' -f 0.8 -wo >  'Tables/CNV_sequenza_selected_cal01_genes.bed'
# bedtools intersect -a 'Tables/gencode.v46.annotation.genes.bed' -b 'Tables/CNV_sequenza_selected_h2301.bed' -f 0.8 -wo >  'Tables/CNV_sequenza_selected_h2301_genes.bed'
# bedtools intersect -a 'Tables/gencode.vM36.annotation_add_hKRAS_G12C.genes.bed' -b 'Tables/CNV_sequenza_selected_mep01.bed' -f 0.8 -wo >  'Tables/CNV_sequenza_selected_mep01_genes.bed'
# 

cnv.genes.list = lapply(cnv.list, function(x){
  model = unique(x$model)
  print(length(unique(x$id)))
  
  segments.genes = read.table(paste0('Tables/CNV_sequenza_selected_', model, '_genes.bed'))
  colnames(segments.genes) = c('chr', 'start', 'end', 'gene_id', 'gene_name', 'gene_type', 'chromosome', 'start.pos', 'end.pos', 'id', 'overlap')
  segments.genes$gene_id_nv = sapply(strsplit(segments.genes$gene_id, '\\.'), '[[', 1)
  print(length(unique(segments.genes$id)))
  
  
  cnv.genes = left_join(x[,c('Bf', 'N.BAF','sd.BAF', 'depth.ratio', 'N.ratio', 'sd.ratio', 'CNt', 'A', 'B', 'LPP', 'logR', 'name', 'comparison', 'barcode', 'model', 'id')], segments.genes)
  cnv.genes$gene_annotated = ifelse(is.na(cnv.genes$gene_id), FALSE, TRUE)
  cnv.genes = cnv.genes %>% mutate(category = case_when(
    CNt > 2 ~ 'Amplification', 
    CNt < 2 ~ "Loss",
    CNt == 2 & A == 1 & B == 1 ~ "Neutral",
    CNt == 2 &  B == 0 ~ "LOH",
    .default = NA))
  
})



# savings -----
cnv.genes.all = bind_rows(cnv.genes.list)
cnv.genes.all$condition = sapply(strsplit(cnv.genes.all$barcode, '\\-'), '[[', 5)

cnv.genes.all.pc = subset(cnv.genes.all, gene_type == 'protein_coding')
cnv.genes.all.pc$condition = sapply(strsplit(cnv.genes.all.pc$barcode, '\\-'), '[[', 5)


saveRDS(cnv.genes.all, 'Rdata/CNV_sequenza_gene_level.rds')
saveRDS(cnv.genes.all.pc, 'Rdata/CNV_sequenza_PC_gene_level.rds')


# 
# mouse = cnv.genes.list$mep01
# mouse = subset(mouse, gene_type == 'protein_coding')
# 
# table(mouse$barcode, mouse$category)
# 
# 
# # Amplification   LOH  Loss Neutral
# # caks-mep01-53-050-000         20652   125    10      15
# # caks-mep01-53-050-010         20616   111     4      11
# # caks-mep01-53-050-020         20818    20    10      22
# 
# 
# cnt.mouse = cnv.genes.all.pc[which(cnv.genes.all.pc$model == 'mep01'),c('barcode', 'CNt', 'category','gene_id','gene_name','chr', 'start', 'end')]
# cncnt.mouset = reshape2::dcast(cnt.mouse, category + chr + start + end  + gene_id + gene_name ~ barcode, value.var = 'CNt')



# Depmap check ------------

# 25Q3 depmap release
depmap = read.csv('Tables/OmicsGlobalSignatures.csv')
depmap_segments     = read.csv("Tables/OmicsCNSegmentsWGS.csv")

## H23 -----------
depmap_h23 = subset(depmap, ModelID == "ACH-000900" )
depmap_segments_h23     = depmap_segments[depmap_segments$ModelID == "ACH-000900", ]


# https://forum.depmap.org/t/how-to-calculate-absolute-copy-number-from-relative-copy-number/1002
depmap_segments_h23$CNt = depmap_segments_h23$SEGMENT_COPY_NUMBER*2.67
depmap_segments_h23$CNt = as.numeric(depmap_segments_h23$CNt)
depmap_segments_h23 = depmap_segments_h23 %>% mutate(category = case_when(
  CNt > 2 ~ 'Amplification', 
  CNt < 2 ~ "Loss",
  CNt == 2  ~ "Neutral/LOH",
  .default = NA
))

depmap_segments_h23$CNt_discrete = round(depmap_segments_h23$CNt)
depmap_segments_h23 = depmap_segments_h23 %>% mutate(category_discrete = case_when(
  CNt_discrete > 2 ~ 'Amplification', 
  CNt_discrete < 2 ~ "Loss",
  CNt_discrete == 2  ~ "Neutral/LOH",
  .default = NA
))


depmap_segments_h23$id = paste0(depmap_segments_h23$CONTIG, '_', depmap_segments_h23$START, '_', depmap_segments_h23$END)
bed = depmap_segments_h23[,c('CONTIG','START','END', 'id')]
bed = unique(bed)
saveBED(bed, paste0('Tables/DepMap_H23_OmicsCNSegmentsWGS', '.bed') )


# bedtools intersect -a 'Tables/gencode.v46.annotation.genes.bed' -b 'Tables/DepMap_H23_OmicsCNSegmentsWGS.bed' -f 0.8 -wo >  'Tables/CNV_DepMap_H23_OmicsCNSegmentsWGS_genes.bed'


segments.genes = read.table( 'Tables/CNV_DepMap_H23_OmicsCNSegmentsWGS_genes.bed')
colnames(segments.genes) = c('chr', 'start', 'end', 'gene_id', 'gene_name', 'gene_type', 'CONTIG', 'START', 'END', 'id', 'overlap')
print(length(unique(segments.genes$id)))


cnv.genes = left_join(depmap_segments_h23[,c('ModelID', 'ModelConditionID','state', 'CNt', 'CNt_discrete', 'category', 'category_discrete', 'id')], segments.genes)
cnv.genes$gene_annotated = ifelse(is.na(cnv.genes$gene_id), FALSE, TRUE)

cnv.genes = subset(cnv.genes, gene_type == 'protein_coding')

saveRDS(cnv.genes, 'Rdata//CNV_DepMap_H23_OmicsCNSegmentsWGS_genes_annotated.rds')






## CALU -----------

depmap_calu = subset(depmap, ModelID == "ACH-000511" )
depmap_segments_calu     = depmap_segments[depmap_segments$ModelID == "ACH-000511", ]


# https://forum.depmap.org/t/how-to-calculate-absolute-copy-number-from-relative-copy-number/1002
depmap_segments_calu$CNt = depmap_segments_calu$SEGMENT_COPY_NUMBER*2.98
depmap_segments_calu$CNt = as.numeric(depmap_segments_calu$CNt)
depmap_segments_calu = depmap_segments_calu %>% mutate(category = case_when(
  CNt > 2 ~ 'Amplification', 
  CNt < 2 ~ "Loss",
  CNt == 2  ~ "Neutral/LOH",
  .default = NA
))

depmap_segments_calu$CNt_discrete = round(depmap_segments_calu$CNt)
depmap_segments_calu = depmap_segments_calu %>% mutate(category_discrete = case_when(
  CNt_discrete > 2 ~ 'Amplification', 
  CNt_discrete < 2 ~ "Loss",
  CNt_discrete == 2  ~ "Neutral/LOH",
  .default = NA
))


depmap_segments_calu$id = paste0(depmap_segments_calu$CONTIG, '_', depmap_segments_calu$START, '_', depmap_segments_calu$END)
bed = depmap_segments_calu[,c('CONTIG','START','END', 'id')]
bed = unique(bed)
saveBED(bed, paste0('Tables/DepMap_CALU_OmicsCNSegmentsWGS', '.bed') )


# bedtools intersect -a 'Tables/gencode.v46.annotation.genes.bed' -b 'Tables/DepMap_CALU_OmicsCNSegmentsWGS.bed' -f 0.8 -wo >  'Tables/CNV_DepMap_CALU_OmicsCNSegmentsWGS_genes.bed'

segments.genes = read.table( 'Tables/CNV_DepMap_CALU_OmicsCNSegmentsWGS_genes.bed')
colnames(segments.genes) = c('chr', 'start', 'end', 'gene_id', 'gene_name', 'gene_type', 'CONTIG', 'START', 'END', 'id', 'overlap')
print(length(unique(segments.genes$id)))


cnv.genes = left_join(depmap_segments_calu[,c('ModelID', 'ModelConditionID','state', 'CNt', 'CNt_discrete', 'category', 'category_discrete', 'id')], segments.genes)
cnv.genes$gene_annotated = ifelse(is.na(cnv.genes$gene_id), FALSE, TRUE)

cnv.genes = subset(cnv.genes, gene_type == 'protein_coding')

saveRDS(cnv.genes, 'Rdata//CNV_DepMap_CALU_OmicsCNSegmentsWGS_genes_annotated.rds')



# Annotating caks with depmap -----

hc.cnv = readRDS('Rdata/CNV_sequenza_PC_gene_level.rds')

h.cnv = subset(hc.cnv, model == 'h2301')
dep.h = readRDS( 'Rdata//CNV_DepMap_H23_OmicsCNSegmentsWGS_genes_annotated.rds')

h.cnv$depmap_CNt               = dep.h$CNt[match(h.cnv$gene_id, dep.h$gene_id)]
h.cnv$depmap_CNt_discrete      = dep.h$CNt_discrete[match(h.cnv$gene_id, dep.h$gene_id)]
h.cnv$depmap_state             = dep.h$state[match(h.cnv$gene_id, dep.h$gene_id)]
h.cnv$depmap_category          = dep.h$category[match(h.cnv$gene_id, dep.h$gene_id)]
h.cnv$depmap_category_discrete = dep.h$category_discrete[match(h.cnv$gene_id, dep.h$gene_id)]
h.cnv[which(h.cnv$condition != '000'),c('depmap_CNt','depmap_CNt_discrete','depmap_state','depmap_category','depmap_category_discrete')] = NA

h.cnv = h.cnv[,c('chromosome','start.pos','end.pos','id', 'Bf','N.BAF','sd.BAF', 'depth.ratio','N.ratio','sd.ratio','CNt','A','B','LPP','logR', 'category'
                 ,'gene_id','gene_name','gene_id_nv', 'gene_type','gene_annotated','chr', 'start', 'end', 'overlap'
                 ,'barcode','model','condition'
                 , 'depmap_CNt', 'depmap_CNt_discrete','depmap_state', 'depmap_category', 'depmap_category_discrete'
                 , 'comparison', 'name')]


saveRDS(h.cnv, 'Rdata//CNV_sequenza_H23_PC_gene_level_annotated_depmap.rds')
write.csv(h.cnv, 'Tables/CNV_sequenza_H23_PC_gene_level_annotated_depmap.csv', row.names = F)


c.cnv = subset(hc.cnv, model == 'cal01')
dep.c = readRDS( 'Rdata//CNV_DepMap_CALU_OmicsCNSegmentsWGS_genes_annotated.rds')

c.cnv$depmap_CNt               = dep.c$CNt[match(c.cnv$gene_id, dep.c$gene_id)]
c.cnv$depmap_CNt_discrete      = dep.c$CNt_discrete[match(c.cnv$gene_id, dep.c$gene_id)]
c.cnv$depmap_state             = dep.c$state[match(c.cnv$gene_id, dep.c$gene_id)]
c.cnv$depmap_category          = dep.c$category[match(c.cnv$gene_id, dep.c$gene_id)]
c.cnv$depmap_category_discrete = dep.c$category_discrete[match(c.cnv$gene_id, dep.c$gene_id)]
c.cnv[which(c.cnv$condition != '000'),c('depmap_CNt','depmap_CNt_discrete','depmap_state','depmap_category','depmap_category_discrete')] = NA

c.cnv = c.cnv[,c('chromosome','start.pos','end.pos','id', 'Bf','N.BAF','sd.BAF', 'depth.ratio','N.ratio','sd.ratio','CNt','A','B','LPP','logR', 'category'
                 ,'gene_id','gene_name','gene_id_nv', 'gene_type','gene_annotated','chr', 'start', 'end', 'overlap'
                 ,'barcode','model','condition'
                 , 'depmap_CNt', 'depmap_CNt_discrete','depmap_state', 'depmap_category', 'depmap_category_discrete'
                 , 'comparison', 'name')]

saveRDS(c.cnv, 'Rdata//CNV_sequenza_CALU_PC_gene_level_annotated_depmap.rds')
write.csv(c.cnv, 'Tables/CNV_sequenza_CALU_PC_gene_level_annotated_depmap.csv', row.names = F)

