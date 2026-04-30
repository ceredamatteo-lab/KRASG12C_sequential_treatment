
## ***************************
##
## Script name: check_mep_ras_snv.R
## Purpose of script:
## Author: mc
## Date Created: 2026-04-26
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



ras.snv = read.table('Tables/mep_RAS_snv_rescue.tsv')

strelka = ras.snv[which(grepl('strelka', ras.snv$V1 )),]
mutect = ras.snv[which(grepl('mutect', ras.snv$V1 )),]
# strelka VCF: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
# mutect VCF: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mep01_caks-mep01-53-050-000	mep01_mef_BL6Nx129s1

colnames(strelka) = c('VCF', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR')
colnames(mutect) = c('VCF', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL')



all = rbind.data.frame(strelka, mutect)
all$tool = ifelse(grepl('strelka', all$VCF), 'strelka', 'mutect')

hKRAS = subset(all, CHROM == 'hKRAS_G12C')
hKRAS$comparison = sapply(strsplit(hKRAS$VCF, '\\/variant_calling'), '[[', 1)
hKRAS$comparison = gsub('/hpcnfs/data/cgb/Ambrogio/WGS260402//Results/sarek/', '', hKRAS$comparison)

write.table(hKRAS, 'Tables/mep_hKRAS_acquired_snv.tsv', col.names = TRUE, row.names = FALSE, sep = '\t')



# ras -------


ras.snv = read.table('Tables/mep_RAS_snv_rescue.tsv')
ras.snv = subset(ras.snv, V2 != 'hKRAS_G12C')

# strelka VCF: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
# mutect VCF: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mep01_caks-mep01-53-050-000	mep01_mef_BL6Nx129s1

strelka = ras.snv[which(grepl('strelka', ras.snv$V1 )),]
colnames(strelka) = c('VCF', 'SAMPLE1', 'SAMPLE2', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR')


strelka.list = list('snv' = strelka[which(grepl('snv', strelka$VCF )),]
                    , 'indels' = strelka[which(grepl('indels', strelka$VCF )),])

strelka.list = lapply(names(strelka.list), function(x, y ){
  strelka = y[[x]]
  strelka = as.data.frame(parse_csq(strelka) )
  strelka$ID = paste0(strelka$CHROM, "_", strelka$POS, "_", strelka$REF, "_", strelka$ALT)
  format_fields      = strsplit(strelka$FORMAT[1], ":")[[1]]
  strelka = strelka %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  strelka = strelka %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  
  if(x == 'snv'){
  strelka = strelka %>%
    mutate(
      refCounts = case_when(REF == "A" ~ TUMOR_AU, REF == "C" ~ TUMOR_CU, REF == "G" ~ TUMOR_GU, REF == "T" ~ TUMOR_TU)
      , altCounts = case_when(ALT == "A" ~ TUMOR_AU, ALT == "C" ~ TUMOR_CU, ALT == "G" ~ TUMOR_GU, ALT == "T" ~ TUMOR_TU)
     
      , ref1 = as.numeric(sapply(strsplit(refCounts, "\\,"), '[[', 1))
      , alt1 = as.numeric(sapply(strsplit(altCounts, "\\,"), '[[', 1))
      
      , TUMOR_AF = alt1 / (alt1 + ref1)
      , TUMOR_DP = alt1 + ref1
      
    )} else if(x == 'indels'){
      strelka = strelka %>%
        mutate(
          ref1   = as.numeric(sapply(strsplit(TUMOR_TAR, "\\,"), '[[', 1))
          , alt1 = as.numeric(sapply(strsplit(TUMOR_TIR, "\\,"), '[[', 1))
          
          , TUMOR_AF = alt1 / (alt1 + ref1)
          , TUMOR_DP = alt1 + ref1
        )
    }
  
  return(strelka)
  
  
}, y = strelka.list)
names(strelka.list) = c('snv', 'indels')

strelka.parsed = bind_rows(strelka.list, .id = 'type')


mutect = ras.snv[which(grepl('mutect', ras.snv$V1 )),]
colnames(mutect) = c('VCF', 'SAMPLE1', 'SAMPLE2', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'FORMAT1', 'FORMAT2')

is_normal = function(x) grepl('mef|050-000|NORMAL', x)

mutect$NORMAL = ifelse(is_normal(mutect$SAMPLE1), mutect$FORMAT1, mutect$FORMAT2)
mutect$TUMOR  = ifelse(is_normal(mutect$SAMPLE1), mutect$FORMAT2, mutect$FORMAT1)

mutect = as.data.frame(parse_csq(mutect) )
mutect$ID = paste0(mutect$CHROM, "_", mutect$POS, "_", mutect$REF, "_", mutect$ALT)

mutect.list = list('snv' = subset(mutect, FORMAT == "GT:AD:AF:DP:F1R2:F2R1:FAD:SB")
                    , 'indels' = subset(mutect, FORMAT == "GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB"))

mutect.list = lapply(mutect.list, function(x){

  mutect = x
  format_fields      = strsplit(mutect$FORMAT[1], ":")[[1]]
  mutect = mutect %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  mutect = mutect %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  return(mutect)
})

mutect.parsed = bind_rows(mutect.list, .id = 'type')






all = rbind.data.frame(strelka.parsed[,c('VCF', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL', 'TUMOR_AF', 'TUMOR_DP')]
                       , mutect.parsed[,c('VCF', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL', 'TUMOR_AF', 'TUMOR_DP')])
all$tool = ifelse(grepl('strelka', all$VCF), 'strelka', 'mutect')

hKRAS = subset(all, CHROM == 'hKRAS_G12C')
hKRAS$comparison = sapply(strsplit(hKRAS$VCF, '\\/variant_calling'), '[[', 1)
hKRAS$comparison = gsub('/hpcnfs/data/cgb/Ambrogio/WGS260402//Results/sarek/', '', hKRAS$comparison)

write.table(hKRAS, 'Tables/mep_hKRAS_acquired_snv.tsv', col.names = TRUE, row.names = FALSE, sep = '\t')











# MODIFIER are genes in a range of 5kb of the reported positions. Thus, I exclude those.


all = rbind.data.frame(strelka.parsed[,c('VCF', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL', 'TUMOR_AF', 'TUMOR_DP'
                                         , 'SYMBOL', 'IMPACT',  "Amino_acids", "Codons", "Existing_variation", 'HGNC_ID', "SIFT", "PolyPhen" )]
                       , mutect.parsed[,c('VCF', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL', 'TUMOR_AF', 'TUMOR_DP'
                                          , 'SYMBOL', 'IMPACT',  "Amino_acids", "Codons", "Existing_variation", 'HGNC_ID', "SIFT", "PolyPhen" )]
)
all$tool = ifelse(grepl('strelka', all$VCF), 'strelka', 'mutect')

all.ras = subset(all, SYMBOL %in% c('Hras', 'Nras', "Mras", "Kras") & IMPACT != 'MODIFIER')
all.ras$comparison = sapply(strsplit(all.ras$VCF, '\\/annotation'), '[[', 1)
all.ras$comparison = gsub('/hpcnfs/data/cgb/Ambrogio/WGS260402//Results/sarek/', '', all.ras$comparison)


write.table(all.ras, 'Tables/mep_ras_acquired_snv.tsv', col.names = TRUE, row.names = FALSE, sep = '\t')


all.ras$downsampling = ifelse(grepl('downsampling', all.ras$comparison), FALSE, TRUE)
all.ras$control = ifelse(grepl('vs_mef', all.ras$comparison), 'mef_BL6Nx129s1', 'caks-mep01-53-050-000')
all.ras$sample = sapply(strsplit(all.ras$comparison, '\\_v|\\_downsampledBL6'), '[[', 1)
all.ras$sample = gsub('vc__no_downsampling/', '', all.ras$sample)


pdf('Figures/ras_acquired_mutations.pdf', width = unit(7, 'cm'), height = unit(7, 'cm'))
ggplot(all.ras, aes(sample, paste0(SYMBOL, '_', ID), col  = as.numeric(TUMOR_AF), shape = FILTER )) + geom_point() + 
  scale_colour_viridis_c() + facet_grid(control~downsampling) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
dev.off()




