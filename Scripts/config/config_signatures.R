## ***************************
##
## Script name: config_signatures.R
## Purpose of script:
## Author: mc
## Date Created: 2026-04-17
## R version: R version 4.4.2 (2024-10-31)
##
## ***************************
##
## Notes:
##   
## ***************************


source("Scripts/config/environment.R")
setwd(DATA_DIR)


# functions ------
# NCG annotation
isNCGdriver = function(x, ncg, gene_identifier, organ = NULL, primarySite = NULL){
  # ncg is a list containing 'cgc' and 'healthy' tables from NCG like:
  ### ncg = list(cgc = read.delim2(paste0(GIT_LOCAL_DIR, 'Pipelines/data/NCG7/NCG_cancerdrivers_annotation_supporting_evidence.tsv'))
  ###           , healthy  = read.delim2(paste0(GIT_LOCAL_DIR,'Pipelines/data/NCG7/NCG_healthydrivers_annotation_supporting_evidence.tsv')) )
  
  # gene_identifier can be: entrez | symbol
  
  # organ can be: 
  ### "Hematologic and lymphatic" "Multiple" "Urologic"  "Gynecologic"   "Soft tissue" "Developmental gastrointestinal"
  ### "Core gastrointestinal" "Thoracic"  "Head and neck" "Skin" 
  #### note:  organ_system of cancer genes and healthy drivers names match
  
  
  # primarySite can be: 
  ###  "blood"  "multiple"  "bladder"   "breast"  "brain"   "pancreas"  "colorectal"                  
  ###  "esophagus"  "small_intestine"   "hepatobiliary"   "pan-gynecological and breast" "stomach"  "lung"    "head_and_neck"               
  ###  "prostate" "soft_tissue" "uterus"    "peripheral_nervous_system"    "thymus"     "skin"  "kidney"                      
  ###  "adrenal_gland"  "pan-gastric"   "thyroid"   "bone"   "parathyroid_gland"    "ovary"   "penis"                       
  ###  "cervix"   "testis"  "uvea"  "retina"  "pleura"    "vascular_system"     
  #### note:  primary site of cancer genes and organ site of healthy drivers names match in exception of "colorectal" and "colon" 
  
  if(gene_identifier == 'entrez'){
    
    x$isNCG_cg = ifelse(x$entrez %in% unique(ncg$cgc$entrez), 'NCG_cg','rest')
    x$isNCG_hd = ifelse(x$entrez %in% unique(ncg$healthy$entrez), 'NCG_hd','rest')
    x$isNCG_ccd = ifelse(x$entrez %in% unique(subset(ncg$cgc, type == 'Canonical Cancer Driver' )$entrez), 'NCG_ccd','rest')
    
    if(!is.null(primarySite)){x$ncg_primary_site = ifelse(x$entrez %in% unique(subset(ncg$cgc, primary_site == primarySite )$entrez), primarySite ,'rest')}
  }
  if(gene_identifier == 'symbol'){
    x$isNCG_cg = ifelse(x$gene_name %in% unique(ncg$cgc$symbol), 'NCG_cg','rest')
    x$isNCG_hd = ifelse(x$gene_name %in% unique(ncg$healthy$symbol), 'NCG_hd','rest')
    x$isNCG_ccd = ifelse(x$gene_name %in% unique(subset(ncg$cgc, type == 'Canonical Cancer Driver' )$symbol), 'NCG_ccd','rest')
    
    if(!is.null(primarySite)){
      x$isNCG_cg_primary_site = ifelse(x$gene_name %in% unique(subset(ncg$cgc, primary_site == primarySite )$symbol), primarySite ,'rest')
      
      if (primarySite == 'colorectal'){
        organSite = 'colon'
        x$isNCG_hd_organ_site = ifelse(x$gene_name %in% unique(subset(ncg$healthy, organ_site == primarySite )$symbol), primarySite ,'rest')
      } else { x$isNCG_hd_organ_site = ifelse(x$gene_name %in% unique(subset(ncg$healthy, organ_site == primarySite )$symbol), primarySite ,'rest') }
    }  
    
    if(!is.null(organ)){
      x$isNCG_cg_organ_system = ifelse(x$gene_name %in% unique(subset(ncg$cgc, organ_system == organ )$symbol), oragn ,'rest')
      if (primarySite == 'colorectal'){
        organSite = 'colon'
        x$isNCG_hd_organ_system = ifelse(x$gene_name %in% unique(subset(ncg$healthy, organ_system == organ )$symbol), organ ,'rest')
      } else { x$ncg_hd_organ_system = ifelse(x$gene_name %in% unique(subset(ncg$healthy, organ_system == organ )$symbol), organ ,'rest') }
    }  
    
    
    
  }
  return(x)
}


isGC = function(x, gc_file = paste0(GIT_LOCAL_DIR,"/Pipelines/GeneCategories/Rdata/gene_categories.20240913.rds"), gene_identifier){
  # gc_file is gene categories file in Pipelines 
  # gene_identifier can be gene_name | entrez | gene_id
  
  gc = readRDS(gc_file)
  gc$gene_id_nv = gc$gene_id
  
  categories = c("RBP", "EM", "TF", "ETM", "PT") # add more in case
  subtypes   = c("splicing", "pioneer") # add more in case 
  
  sel_categories = categories[categories %in% unique(gc$category)]
  sel_subtypes = subtypes[subtypes %in% unique(gc$subtype)]
  
  skipped_categories = setdiff(categories, sel_categories)
  skipped_subtypes   = setdiff(subtypes,   sel_subtypes)
  
  if(length(skipped_categories) > 0) message("[*] Categories not included: ", paste(skipped_categories, collapse=", "))
  if(length(skipped_subtypes) > 0) message("[*] Subtypes not included: ", paste(skipped_subtypes, collapse=", "))
  
  colnames(x)[grep("gene_id_nv|gene_id_no_version$", colnames(x))] = "gene_id_nv"
  
  for (variable in sel_categories) {x[[paste0("is", variable)]] = ifelse(x[[gene_identifier]] %in% unique(gc[[gene_identifier]][gc$category == variable]), variable, "rest")}
  for (variable in sel_subtypes) {x[[paste0("is", variable)]] = ifelse(x[[gene_identifier]] %in% unique(gc[[gene_identifier]][gc$subtype == variable]), variable, "rest")}
  
  colnames(x)[grep("issplicing", colnames(x))] = "isSplicing"
  colnames(x)[grep("ispioneer", colnames(x))] = "isPTF"
  
  return(x)
}




# signatures --------
message("[*] kras_signature = KRAS regulated genes in PDAC (Klomb Science 2024) ")
# KRAS regulated genes in PDAC (Klomb Science 2024)  
kras_signature = readRDS(paste0(GIT_LOCAL_DIR,"Pipelines/GeneCategories/Rdata/KRAS_signature.20240806.rds"))
kras_signature$gene_id_nv = sapply(strsplit(kras_signature$gene_id, '\\.'), '[[', 1)



message("[*] RAS84 = RAS84 oncogenic genes (East NatCom 2022) ")
# RAS84 oncogenic genes (East NatCom 2022)  
RAS84 = read.csv(paste0(GIT_LOCAL_DIR,"Pipelines/GeneCategories/KRAS_signatures/RAS_84_PhilipEast_et_al_NatCom_2022.csv"))



message("[*] mapk_erk_signature = MAPK and ERK1/2 GO signature (East NatCom 2022) GO human")
# MAPK/ERK SIGNATURE (East NatCom 2022)  
mapk_erk_signature = readRDS(paste0(GIT_LOCAL_DIR,"GIRP/Rdata/mapk_erk_signature_East_NatCom_2022.rds"))
# mapk_signature=subset(msigdbr(species = "human", category = 'C5', subcategory = "BP") ,gs_exact_source=="GO:0000165")
# saveRDS(mapk_signature, paste0(DATA_DIR, "Rdata/mapk_signature.rds"))
mapk_signature = readRDS(paste0(GIT_LOCAL_DIR, "GIRP/Rdata/mapk_signature.rds"))
# erk_signature=subset(msigdbr(species='human',category = 'C5', subcategory = 'BP'), gs_exact_source=='GO:0070371')
# saveRDS(erk_signature, paste0(DATA_DIR, "Rdata/erk_signature.rds"))
erk_signature = readRDS(paste0(GIT_LOCAL_DIR, "GIRP/Rdata/erk_signature.rds"))




message("[*] NCG7")
ncg = list(cgc        = read.delim2(paste0(GIT_LOCAL_DIR, 'Pipelines/data/NCG7/NCG_cancerdrivers_annotation_supporting_evidence.tsv')) 
           , healthy  = read.delim2(paste0(GIT_LOCAL_DIR,'Pipelines/data/NCG7/NCG_healthydrivers_annotation_supporting_evidence.tsv')) 
)




