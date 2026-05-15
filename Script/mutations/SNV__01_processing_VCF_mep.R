## ***************************
##
## Script name: mep01_processing_vcf_files_SNV_no_downsampling.R
## Purpose of script: filtering, annotations (VEP + annovar), pathogeneicity, oncogenic patheways
## Author: mc
## Date Created: 2026-04-18
##
## ***************************
##
## Notes: based on point_00_mutation_catalogue.R, done on @PERALTA
## conda env /hpcnfs/data/cgb/conda_envs/R4
##   
## ***************************

# setting -------
DATA_DIR = "/hpcnfs/data/cgb/Ambrogio/WGS260402/"  
source('/hpcnfs/data/cgb/Ambrogio/WGS260402/Scripts/config_filtering_point__PERALTA.R')


# libraries --------
library(openxlsx)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(readr)



# data -----
model = 'caks-mep01'
samples = c('caks-mep01-53-050-000', 'caks-mep01-53-050-010','caks-mep01-53-050-020')
control = 'mef_BL6Nx129s1'
is_mouse = TRUE

message('[*] model: ' , model)
message('[*] samples: ' , paste0(samples, collapse = ' '))
message('[*] control: ' , control)
message('[*] is_mouse: ' , is_mouse)


mutect.path = paste0(DATA_DIR, "Results/vep/vc__no_downsampling/annotation/mutect2/")
file_suffix_mutect2_VEP = ".mutect2.filtered_VEP.ann.vcf.gz"

strelka.path = paste0(DATA_DIR, "Results/sarek/vc__no_downsampling/")
file_suffix_strelka_snvs_VEP = ".strelka.somatic_snvs_VEP.ann.vcf.gz"
file_suffix_strelka_indels_VEP = ".strelka.somatic_indels_VEP.ann.vcf.gz"


if (! is_mouse){
  onc_pws.file      = paste0(DATA_DIR, "Tables/oncogenic_sig_patwhays.tsv")
  onc_pws           = read_tsv(onc_pws.file)
  colnames(onc_pws) = c("oncogenic_pathway","SYMBOL", "OG_TSG")
  
  kegg_gs_list.file = paste0(DATA_DIR, 'Tables/kegg_gs_list.csv')
  kegg_gs_list      = read.csv(kegg_gs_list.file)
  
  hallmarks_gs_list.file = paste0(DATA_DIR, 'Tables/hallmarks_gs_list.csv')
  hallmarks_gs_list      = read.csv(hallmarks_gs_list.file)
                                  
                                  
}else{
  onc_pws.file      = paste0(DATA_DIR, "Tables/oncogenic_sig_patwhays_mouse.csv") # in sv_00....R there is script to obtain this file 
  onc_pws           = read.csv(onc_pws.file) # in sv_00....R there is script to obtain this file 
  onc_pws           = onc_pws[, c("oncogenic_pathway", "mouse_gene_symbol", "OG_TSG")]
  colnames(onc_pws) = c("oncogenic_pathway","SYMBOL", "OG_TSG")
  onc_pws           = onc_pws[!is.na(onc_pws$SYMBOL),]
  
  kegg_gs_list.file = paste0(DATA_DIR, 'Tables/kegg_gs_list_mm.csv')
  kegg_gs_list      = read.csv(kegg_gs_list.file)
  
  hallmarks_gs_list.file = paste0(DATA_DIR, 'Tables/hallmarks_gs_list_mm.csv')
  hallmarks_gs_list      = read.csv(hallmarks_gs_list.file)
}






# Mutect2 ------
message('[*] 1. Parsing Mutect2 ... ' )

mutect2_list  =  list()

for (sample in samples){
  
  message('[*] Processing ', sample, '...')
  file_path_mutect2        = paste0(mutect.path, sample, '/', sample, "_vs_", control, file_suffix_mutect2_VEP )
  
  message('[*] Input: ', file_path_mutect2)
  vcf_mutect2              = fread(file_path_mutect2, skip = "#CHROM", header = TRUE, sep = "\t")
  colnames(vcf_mutect2)[1] = 'CHROM'
  
  vcf_mutect2$ID     = paste0(vcf_mutect2$CHROM, "_", vcf_mutect2$POS, "_", vcf_mutect2$REF, "_", vcf_mutect2$ALT)
  vcf_mutect2$sample = sub("_.*", "", sample)
  
  # change colnames
  colnames(vcf_mutect2)[grepl(control, colnames(vcf_mutect2))] = "NORMAL"
  colnames(vcf_mutect2)[grepl('caks', colnames(vcf_mutect2))]  = "TUMOR"
  
  vcf_mutect2_snv    = subset(vcf_mutect2, FORMAT == "GT:AD:AF:DP:F1R2:F2R1:FAD:SB")
  vcf_mutect2_indels = subset(vcf_mutect2, FORMAT == "GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB")
  
  
  message('[*] SNVs ...')
  # columns in FORMAT x snv
  format_fields = strsplit(vcf_mutect2_snv$FORMAT[1], ":")[[1]]
  # separate TUMOR
  vcf_mutect2_snv.2 = vcf_mutect2_snv %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  # separate NORMAL
  vcf_mutect2_snv.2 = vcf_mutect2_snv.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  
  message('[*] InDels ...')
  # columns in FORMAT x indels
  format_fields = strsplit(vcf_mutect2_indels$FORMAT[1], ":")[[1]]
  # separate TUMOR
  vcf_mutect2_indels.2 = vcf_mutect2_indels %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  # separate NORMAL
  vcf_mutect2_indels.2 = vcf_mutect2_indels.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  
  message('[*] Binding SNVs + InDels ...')
  vcf_mutect2_all = bind_rows(vcf_mutect2_snv.2, vcf_mutect2_indels.2)
  message("[*] Mutect unique ID: ", length(unique(vcf_mutect2_all$ID) ))

  
  message('[*] Parsing CSQ field ...')
  vcf_mutect2_all_vep = as.data.frame(parse_csq(vcf_mutect2_all) )
  message("[*] Mutect after parsing vep annotaion unique ID: ", length(unique(vcf_mutect2_all_vep$ID)))
  
  mutect2_list[[sample]] <- vcf_mutect2_all_vep
  
}


# Strelka -------
message('[*] 2. Parsing Strelka ... ' )

strelka_list  =  list()

for (sample in samples){
  message('[*] Processing ', sample, '...')
  message('[*] SNVs ... ')
  
  file_path_strelka_snvs = paste0(strelka.path, sample, "_vs_", control, '/annotation/strelka/', sample, "_vs_", control, '/', sample, "_vs_", control, file_suffix_strelka_snvs_VEP )
  
  message('[*] Input: ', file_path_strelka_snvs)
  vcf_strelka_snvs              = fread(file_path_strelka_snvs, skip = "#CHROM", header = TRUE, sep = "\t")
  colnames(vcf_strelka_snvs)[1] = 'CHROM'
  
  vcf_strelka_snvs$ID = paste0(vcf_strelka_snvs$CHROM, "_", vcf_strelka_snvs$POS, "_", vcf_strelka_snvs$REF, "_", vcf_strelka_snvs$ALT)

  # columns in FORMAT
  format_fields      = strsplit(vcf_strelka_snvs$FORMAT[1], ":")[[1]]
  # separate TUMOR
  vcf_strelka_snvs.2 = vcf_strelka_snvs %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  ## separate NORMAL
  vcf_strelka_snvs.2 = vcf_strelka_snvs.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  message('[*] Compute AF ... ')

  # compute AF https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md
  vcf_strelka_snvs.2 = vcf_strelka_snvs.2 %>%
    mutate(
      refCounts = case_when(
        REF == "A" ~ TUMOR_AU
        , REF == "C" ~ TUMOR_CU
        , REF == "G" ~ TUMOR_GU
        , REF == "T" ~ TUMOR_TU
      )
      
      , altCounts = case_when(
        ALT == "A" ~ TUMOR_AU
        , ALT == "C" ~ TUMOR_CU
        , ALT == "G" ~ TUMOR_GU
        , ALT == "T" ~ TUMOR_TU
      )
      
      , ref1 = as.numeric(sapply(strsplit(refCounts, "\\,"), '[[', 1))
      , alt1 = as.numeric(sapply(strsplit(altCounts, "\\,"), '[[', 1))
      
      , TUMOR_AF_manual = alt1 / (alt1 + ref1)
      , TUMOR_DP_manual = alt1 + ref1
    )
  
  message('[*] InDels ... ')
  file_path_strelka_indels = paste0(strelka.path, sample, "_vs_", control, '/annotation/strelka/', sample, "_vs_", control, '/', sample, "_vs_", control, file_suffix_strelka_indels_VEP )
  
  message('[*] Input: ', file_path_strelka_indels)
  vcf_strelka_indels              = fread(file_path_strelka_indels, skip = "#CHROM", header = TRUE, sep = "\t")
  colnames(vcf_strelka_indels)[1] = 'CHROM'
  
  vcf_strelka_indels$ID = paste0(vcf_strelka_indels$CHROM, "_", vcf_strelka_indels$POS, "_", vcf_strelka_indels$REF, "_", vcf_strelka_indels$ALT)

  # columns in FORMAT
  format_fields <- strsplit(vcf_strelka_indels$FORMAT[1], ":")[[1]]
  ## separate TUMOR
  vcf_strelka_indels.2 = vcf_strelka_indels %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  ## separate NORMAL
  vcf_strelka_indels.2 = vcf_strelka_indels.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  message('[*] Compute AF ... ')
  # compute AF 
  vcf_strelka_indels.2 = vcf_strelka_indels.2 %>%
    mutate(
      ref1   = as.numeric(sapply(strsplit(TUMOR_TAR, "\\,"), '[[', 1))
      , alt1 = as.numeric(sapply(strsplit(TUMOR_TIR, "\\,"), '[[', 1))
      
      , TUMOR_AF_manual = alt1 / (alt1 + ref1)
      , TUMOR_DP_manual = alt1 + ref1
    )
  
  
  vcf_strelka        =  bind_rows(vcf_strelka_snvs.2, vcf_strelka_indels.2)
  vcf_strelka$sample = sub("_.*","",sample)
  
  message("Strelka unique ID: ", length(unique(vcf_strelka$ID)))
  
  vcf_strelka_vep = as.data.frame(parse_csq(vcf_strelka)) 
  message("Strelka after parsing vep annotaion unique ID: ", length(unique(vcf_strelka_vep$ID)) )
  
  strelka_list[[sample]] <- vcf_strelka_vep
}  


mutect2_df = do.call(rbind, mutect2_list)
strelka_df = do.call(rbind, strelka_list)


saveRDS(mutect2_df, paste0(DATA_DIR, "Rdata/", model, '_vs_', control, '_no_downsampling__', 'mutect2_no_filters.rds'))
saveRDS(strelka_df, paste0(DATA_DIR, "Rdata/", model, '_vs_', control, '_no_downsampling__', 'strelka_no_filters.rds'))

lapply(mutect2_list, function(x){length(unique(x$ID))})
lapply(strelka_list, function(x){length(unique(x$ID))})


# PASS filtering -------
message('[*] 3. PASS filtering ... ')

mutect2_df_pass = subset(mutect2_df, FILTER == "PASS")
strelka_df_pass = subset(strelka_df, FILTER == "PASS")

lapply(unique(mutect2_df_pass$sample), function(x){length(unique(mutect2_df_pass$ID[mutect2_df_pass$sample==x]))})
lapply(unique(strelka_df_pass$sample), function(x){length(unique(strelka_df_pass$ID[strelka_df_pass$sample==x]))})


# AF -------
message('[*] 4. AF ≥ 0.05 ... ')

mutect2_df_af   = subset(mutect2_df_pass, as.numeric(TUMOR_AF) >= 0.05)
strelka_df_af   = subset(strelka_df_pass, as.numeric(TUMOR_AF_manual) >= 0.05)

lapply(unique(mutect2_df_af$sample), function(x){length(unique(mutect2_df_af$ID[mutect2_df_af$sample==x]))})
lapply(unique(strelka_df_af$sample), function(x){length(unique(strelka_df_af$ID[strelka_df_af$sample==x]))})


# DP -------
message('[*] 5. DP ≥ 10 ... ')
mutect2_df_dp   = subset(mutect2_df_af, as.numeric(TUMOR_DP) >= 10)
strelka_df_dp   = subset(strelka_df_af, as.numeric(TUMOR_DP_manual) >= 10)

lapply(unique(mutect2_df_dp$sample), function(x){length(unique(mutect2_df_dp$ID[mutect2_df_dp$sample==x]))})
lapply(unique(strelka_df_dp$sample), function(x){length(unique(strelka_df_dp$ID[strelka_df_dp$sample==x]))})



# Overlap -----
message('[*] 6. Mutect2+strelka overlapping mutations ... ')
overlap_mutect2_strelka_list  =  list()

for (sample in unique(mutect2_df_dp$sample)){ 
  message('[*] Sample: ', sample, '... ')
  
  print(sample)

  mutect_sample  = mutect2_df_dp[mutect2_df_dp$sample==sample,]
  strelka_sample = strelka_df_dp[strelka_df_dp$sample==sample,]
  
  common = intersect(mutect_sample$ID, strelka_sample$ID)
  vcf_overlap_vep = subset(mutect_sample, ID %in% common  )

  overlap_mutect2_strelka_list[[sample]] = vcf_overlap_vep
  message("Overlap unique ID: ", length(unique(vcf_overlap_vep$ID))) 
  
}


lapply(overlap_mutect2_strelka_list, function(x){length(unique(x$ID))})

overlap_df = do.call(rbind, overlap_mutect2_strelka_list)
saveRDS(overlap_df, paste0(DATA_DIR, "Rdata/", model,  '_vs_', control, '_no_downsampling__', 'filtered_overlap_mutect2_strelka.rds'))




# Autosomal chromsomes ----
message('[*] 7. Autosomal chromosomes selection ... ')
table(overlap_df$CHROM)
selected = subset(overlap_df, !CHROM %in% c("chrX", "chrY", "chrM") & grepl("^chr[0-9]+$", CHROM))
table(selected$CHROM)
lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample == x]))})

# PC genes ----
message('[*] 8. PC genes ... ')
selected = subset(selected, BIOTYPE == "protein_coding")
lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample==x]))})
saveRDS(selected, paste0(DATA_DIR, "Rdata/", model, '_vs_', control, '_no_downsampling__', 'filtered_overlap_mutect2_strelka_PCgenes_autochrom.rds'))

# impact HIGH or MODERATE ----
message('[*] 9. impact HIGH or MODERATE ... ')
selected = subset(selected, IMPACT %in% c("HIGH", "MODERATE"))

lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample==x]))})
saveRDS(selected, paste0(DATA_DIR, "Rdata/", model, '_vs_', control, '_no_downsampling__', 'filtered_overlap_mutect2_strelka_PCgenes_autochrom_impact.rds'))



if (! is_mouse){
# Damaging level (human cell lines only) ------
message('[*] Damaging level ... ')
selected = assign_damaging_level(df = selected, SIFT_col = "SIFT", PolyPhen_col = "PolyPhen", ClinVar_col = "CLIN_SIG") # stringent or permissive label

}else{message('[*] Damaging level not available for mouse ')}



# Pathways annotation -----------------
message('[*] 10. Pathways annotation ... ')

message('[*] Oncogenic pathways: ', onc_pws.file)
message('[*] KEGG pathways: ', kegg_gs_list.file)
message('[*] Hallmarks pathways: ', hallmarks_gs_list.file)


selected = left_join(selected, onc_pws, by = "SYMBOL") 
selected = add_kegg_info(anno_df = selected, kegg_gs_list = kegg_gs_list)
selected = add_hallmarks_info(anno_df = selected, hallmarks_gs_list = hallmarks_gs_list)

selected[!is.na(selected$hallmark),]

lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample==x]))})

saveRDS(selected, paste0(DATA_DIR, "Rdata/", model,  '_vs_', control, '_no_downsampling__', 'filtered_overlap_mutect2_strelka_PCgenes_autochrom_impact_pathways.rds'))
