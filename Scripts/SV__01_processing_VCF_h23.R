## ***************************
##
## Script name: caks_h2301_processing_vcf_files_SV.R
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
source('/hpcnfs/data/cgb/Ambrogio/WGS260402/Scripts/config_filtering_sv__PERALTA.R')


# libraries --------
library(openxlsx)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(readr)



# data -----
model = 'caks-h2301'
samples = c('caks-h2301-50-050-000_downsampledNormal.filtered', 'caks-h2301-50-050-010_downsampledNormal.filtered', 'caks-h2301-50-050-020_downsampledNormal.filtered')
control = 'mgpf-af001-10-050-000.filtered'
is_mouse = FALSE


message('[*] model: ' , model)
message('[*] samples: ' , paste0(samples, collapse = ' '))
message('[*] control: ' , control)
message('[*] is_mouse: ' , is_mouse)

svdb.path = paste0(DATA_DIR, "Results/SVDB/")
file_suffix_svdb = "_merged_caller.vcf.gz"

manta.path = paste0(DATA_DIR,"Results/vep/annotation/manta/")


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



# SVDB ------------
message('[*] 1. Parsing SVDB ... ' )

svdb_list = list()

for (sample in samples){

  message('[*] Processing ', sample, '...')
  file_path_svdb = paste0(svdb.path, sample, "_vs_", control, file_suffix_svdb)
  
  message('[*] Input: ', file_path_svdb)
  vcf_svdb              = fread(file_path_svdb, skip = "#CHROM", header = TRUE, sep = "\t")
  colnames(vcf_svdb)[1] = 'CHROM'

  colnames(vcf_svdb)[colnames(vcf_svdb)=="ID"] = "ID_tool"
  vcf_svdb$ID     = paste0(vcf_svdb$CHROM, "_", vcf_svdb$POS, "_", vcf_svdb$REF, "_", vcf_svdb$ALT)
  vcf_svdb$sample = sub("_.*", "", sample)

  # change colnames
  colnames(vcf_svdb)[grepl(sub(".filtered*","",control), colnames(vcf_svdb))] = "NORMAL"
  colnames(vcf_svdb)[grepl(sub("_.*","",sample), colnames(vcf_svdb))] = "TUMOR"
  
  
  message('[*] Extracting overlapping SVs: SUPP_VEC=11 ... ' )
  vcf_svdb_info_col = vcf_svdb$INFO
  supp_vec_info = sapply(strsplit(vcf_svdb_info_col, '\\;'), function(x){tail(x, n = 1)})
  print(table(supp_vec_info))
  overlap.index = which(supp_vec_info == 'SUPP_VEC=11')
  
  vcf_svdb_overlap = vcf_svdb[overlap.index,]
  
  message("Unique formats in svdb ",sample," df:", length(unique(vcf_svdb_overlap$FORMAT)), ". In the analysis 2 is considered.")
  
  vcf_svdb_pr_only = subset(vcf_svdb_overlap, FORMAT=="PR")
  vcf_svdb_pr_sr = subset(vcf_svdb_overlap, FORMAT=="PR:SR")
  
  # columns in FORMAT x vcf_svdb_pr_only
  format_fields = strsplit(vcf_svdb_pr_only$FORMAT[1], ":")[[1]]
  ## separate TUMOR
  vcf_svdb_pr_only.2 = vcf_svdb_pr_only %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  ## separate NORMAL
  vcf_svdb_pr_only.2 = vcf_svdb_pr_only.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  # columns in FORMAT x vcf_svdb_pr_sr
  format_fields = strsplit(vcf_svdb_pr_sr$FORMAT[1], ":")[[1]]
  ## separate TUMOR
  vcf_svdb_pr_sr.2 = vcf_svdb_pr_sr %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
  ## separate NORMAL
  vcf_svdb_pr_sr.2 = vcf_svdb_pr_sr.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
  
  vcf_svdb_all = bind_rows(vcf_svdb_pr_only.2, vcf_svdb_pr_sr.2)
  message("SVDB unique ID: ", length(unique(vcf_svdb_all$ID)))
  colnames(vcf_svdb_all)[colnames(vcf_svdb_all) == "INFO"] = "INFO_merged"
  
  
  message('[*] Retrieve information from manta ... ' )
  # use INFO manta, which also contains VEP annotation
  
  file_path_manta = paste0(manta.path, sample, "/somaticSV_VEP.ann.vcf.gz")
  message('[*] manta input: ', file_path_manta )
  
  vcf_manta = fread(file_path_manta, skip = "#CHROM", header = TRUE, sep = "\t")
  
  colnames(vcf_manta)[1]                           = 'CHROM'
  colnames(vcf_manta)[colnames(vcf_manta) == "ID"] = "ID_tool"
  vcf_manta$ID                                     = paste0(vcf_manta$CHROM, "_", vcf_manta$POS, "_", vcf_manta$REF, "_", vcf_manta$ALT)
  
  
  message('[*] Joining manta information to instersection set ... ' )
  vcf_svdb_all_info_from_manta = left_join(vcf_svdb_all, vcf_manta[,c("ID", "INFO")], by = "ID")
  
  
  message('[*] Parsing manta info ... ' )
  vcf_svdb_all_vep = as.data.frame(parse_manta_info(vcf_svdb_all_info_from_manta))
  message("SVDB after parsing vep annotaion unique ID: ", length(unique(vcf_svdb_all_vep$ID)))
  
  
  message('[*] Add mate information ... ' )
  
  vcf_svdb_all_vep = vcf_svdb_all_vep %>%
    mutate(
      mateCHROM = ifelse(str_detect(ALT, "\\[|\\]"),
                         str_extract(ALT, "chr[^:\\[\\]]+"),
                         NA),
      matePOS   = ifelse(str_detect(ALT, "\\[|\\]"),
                         as.numeric(str_extract(ALT, "(?<=:)[0-9]+")),
                         NA)
    )
  
  svdb_list[[sample]] = vcf_svdb_all_vep
}


svdb_df = do.call(rbind, svdb_list)
saveRDS(svdb_df, paste0(DATA_DIR, "Rdata/", model,  '_vs_', control, '__', 'overlap_manta_tiddit_no_filters.rds'))
lapply(svdb_list, function(x){length(unique(x$ID))})



# PASS FILTERING --------------------------
message('[*] 2. PASS filtering ... ')
selected = subset(svdb_df, FILTER == "PASS")

lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample==x]))})


## autosomal chr ----
message('[*] 3. Autosomal chromosomes ... ')
table(selected$CHROM)
selected = subset(selected, !CHROM %in% c("chrX", "chrY", "chrM") & grepl("^chr[0-9]+$", CHROM))
table(selected$CHROM)

lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample == x]))})



## protein coding genes-------
message('[*] 4. PC genes ... ')
selected = subset(selected, BIOTYPE == "protein_coding")
lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample == x]))})


### check SYMBOL
# if (is_mouse){
#   genome <- read.table("Tables/gencode.vM36.annotation_add_hKRAS_G12C.genes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# }else{
#   genome <- read.table("Tables/gencode.v46.annotation.genes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# }
# 
# colnames(genome) <- c("chr", "start", "end", "gene_id", "gene_name", "gene_type")
# genome_pc <- genome[genome$gene_type=="protein_coding",]
# 
# prova_check <- annotate_sv_genes_advanced(sv_df=svdb_df_pc_mated, genes_df=genome_pc, upstream = 5000, downstream = 5000)
# prova_check_clean <- prova_check %>%
#   filter(!grepl("^ENSG", gene_name))
# 
# check_gene <- prova_check_clean %>% 
#   group_by(sv_id) %>% 
#   mutate(
#     vep_in_check = gene_name %in% svdb_df_pc_mated$SYMBOL[
#       svdb_df_pc_mated$ID == unique(sv_id)
#     ]
#   ) %>%
#   ungroup()

saveRDS(selected, paste0(DATA_DIR, "Rdata/", model,  '_vs_', control, '__', 'overlap_manta_tiddit_PCgenes_autochrom.rds'))



## impact HIGH or MODERATE ----
message('[*] 5. impact HIGH  ... ')
selected = subset(selected, IMPACT %in% c("HIGH"))

lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample == x]))})

saveRDS(selected, paste0(DATA_DIR, "Rdata/", model,  '_vs_', control, '__', 'overlap_manta_tiddit_PCgenes_autochrom_impact.rds'))


## Pathways ---------------------------
message('[*] 10. Pathways annotation ... ')


message('[*] Oncogenic pathways: ', onc_pws.file)
message('[*] KEGG pathways: ', kegg_gs_list.file)
message('[*] Hallmarks pathways: ', hallmarks_gs_list.file)


selected = left_join(selected, onc_pws, by = "SYMBOL") 
selected = add_kegg_info(anno_df = selected, kegg_gs_list = kegg_gs_list)
selected = add_hallmarks_info(anno_df = selected, hallmarks_gs_list = hallmarks_gs_list)

selected[!is.na(selected$hallmark),]

lapply(unique(selected$sample), function(x){length(unique(selected$ID[selected$sample == x]))})

saveRDS(selected, paste0(DATA_DIR, "Rdata/", model,  '_vs_', control, '__', 'overlap_manta_tiddit_PCgenes_autochrom_impact_pathways.rds'))






# garage -----

###
# tool<-c('manta','tiddit')
# 
# file_suffix_manta_VEP<-".manta.somatic_sv_VEP.ann.vcf.gz"
# file_suffix_tiddit_VEP<-".tiddit_sv_merge_VEP.ann.vcf.gz"
# 
# manta_list  =  list()
# tiddit_list  =  list()
# PARSING RAW --------------
# data_path = "Results/sarek/"
# ## Manta -------------------
# for (sample in samples){
#   print(sample)
#   file_path_manta = paste0(data_path, "annotation/", tool[1], "/",sample, "_vs_", control ,"/",sample, "_vs_", control, file_suffix_manta_VEP)
#   # manta
#   vcf_manta <- fread(file_path_manta, skip = "#CHROM", header = TRUE, sep = "\t")
#   colnames(vcf_manta)[1]<-'CHROM'
#   # ID
#   colnames(vcf_manta)[colnames(vcf_manta)=="ID"] <- "ID_tool"
#   vcf_manta = vcf_manta %>% mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))
#   #sample
#   vcf_manta$sample <- sub("_.*","",sample)
#   # change colnames
#   colnames(vcf_manta)[grepl(sub(".filtered*","",control), colnames(vcf_manta))] <- "NORMAL"
#   colnames(vcf_manta)[grepl(sub("_.*","",sample), colnames(vcf_manta))] <- "TUMOR"
#   print(paste0("Unique formats in manta ",sample," df:", length(unique(vcf_manta$FORMAT)), ". In the analysis 2 is considered."))
#   vcf_manta_pr_only <- subset(vcf_manta, FORMAT=="PR")
#   vcf_manta_pr_sr <- subset(vcf_manta, FORMAT=="PR:SR")
#   # columns in FORMAT x vcf_manta_pr_only
#   format_fields <- strsplit(vcf_manta_pr_only$FORMAT[1], ":")[[1]]
#   ## separate TUMOR
#   vcf_manta_pr_only.2 <- vcf_manta_pr_only %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
#   ## separate NORMAL
#   vcf_manta_pr_only.3 <- vcf_manta_pr_only.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
#   # columns in FORMAT x vcf_manta_pr_sr
#   format_fields <- strsplit(vcf_manta_pr_sr$FORMAT[1], ":")[[1]]
#   ## separate TUMOR
#   vcf_manta_pr_sr.2 <- vcf_manta_pr_sr %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
#   ## separate NORMAL
#   vcf_manta_pr_sr.3 <- vcf_manta_pr_sr.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
#   vcf_manta_all <- bind_rows(vcf_manta_pr_only.3, vcf_manta_pr_sr.3)
#   print(paste0("Manta unique ID: ", length(unique(vcf_manta_all$ID))))
#   vcf_manta_all_vep = vcf_manta_all %>% parse_manta_info() %>% as.data.frame()
#   print(paste0("Manta after parsing vep annotaion unique ID: ", length(unique(vcf_manta_all_vep$ID))))
#   manta_list[[sample]] <- vcf_manta_all_vep
# }
# manta_df <- do.call(rbind, manta_list)
# saveRDS(manta_df, paste0(DATA_DIR, "Rdata/", paste0(sub("_.*","",samples), collapse = "+"),  '_vs_', control, '__', 'manta_no_filters.rds'))
# lapply(manta_list, function(x){length(unique(x$ID))})
# 
# ## tiddit -----------------
# for (sample in samples){
#   print(sample)
#   file_path_tiddit = paste0(data_path, "annotation/", tool[2], "/",sample, "_vs_", control, "/", sample, "_vs_", control, file_suffix_tiddit_VEP)
#   # tiddit
#   vcf_tiddit <- fread(file_path_tiddit, skip = "#CHROM", header = TRUE, sep = "\t")
#   colnames(vcf_tiddit)[1]<-'CHROM'
#   # ID 
#   colnames(vcf_tiddit)[colnames(vcf_tiddit)=="ID"] <- "ID_tool"
#   vcf_tiddit = vcf_tiddit %>% mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))
#   #sample
#   vcf_tiddit$sample <- sub("_.*","",sample)
#   # change colnames
#   colnames(vcf_tiddit)[grepl(sub(".filtered*","",control), colnames(vcf_tiddit))] <- "NORMAL"
#   colnames(vcf_tiddit)[grepl(sub("_.*","",sample), colnames(vcf_tiddit))] <- "TUMOR"
#   print(paste0("Unique formats in tiddit ",sample," df:", length(unique(vcf_tiddit$FORMAT)), ". In the analysis 1 is considered."))
#   # columns in FORMAT x vcf_tiddit_pr_only
#   format_fields <- strsplit(vcf_tiddit$FORMAT[1], ":")[[1]]
#   ## separate TUMOR
#   vcf_tiddit.2 <- vcf_tiddit %>% separate(TUMOR, into = paste0("TUMOR_", format_fields), sep = ":", remove = FALSE)
#   ## separate NORMAL
#   vcf_tiddit.3 <- vcf_tiddit.2 %>% separate(NORMAL, into = paste0("NORMAL_", format_fields), sep = ":", remove = FALSE)
#   print(paste0("Tiddit unique ID: ", length(unique(vcf_tiddit.3$ID))))
#   vcf_tiddit_all_vep = vcf_tiddit.3 %>% parse_tiddit_info() %>% as.data.frame()
#   print(paste0("Tiddit after parsing vep annotaion unique ID: ", length(unique(vcf_tiddit_all_vep$ID))))
#   tiddit_list[[sample]] <- vcf_tiddit_all_vep
# }
# tiddit_df <- do.call(rbind, tiddit_list)
# saveRDS(tiddit_df, paste0(DATA_DIR, "Rdata/", paste0(sub("_.*","",samples), collapse = "+"),  '_vs_', control, '__', 'tiddit_no_filters.rds'))
# lapply(tiddit_list, function(x){length(unique(x$ID))})

