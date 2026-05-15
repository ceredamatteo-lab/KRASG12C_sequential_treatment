# Requirments --------

pkgs = c(
  'msigdbr'
  #,'GOxploreR' #version > 1.2.7
  ,'clusterProfiler'
  #,'pathfindR'
  ,'readr'
  ,'dplyr'
  , 'ggplot2'
  , 'forcats'
)


lib = installed.packages()
installed=pkgs %in% rownames(lib)
not_installed = which(!installed)
if(length(not_installed)>0){
  for(i in not_installed) install.packages(pkgs[i], dependencies=T)
}

lib = installed.packages()
installed=pkgs %in% rownames(lib)
not_installed = which(!installed)
if(length(not_installed)>0){
  for(i in not_installed) BiocManager::install(pkgs[i])
}

for(i in pkgs) suppressPackageStartupMessages(library(i, character.only =T))

ORA = function(gene_identifier, gene_set_list, uni=NULL, quiet=F, unique_gene = T){
  enr = list()
  for (i in names(gene_set_list)) {
    if(!quiet) print(i)
    if(!is.null(uni)){
      if(unique_gene){
        y = enricher(unique(gene_identifier), TERM2GENE = gene_set_list[[i]]
                     , pAdjustMethod = "BH"
                     , pvalueCutoff = 1
                     , qvalueCutoff = 1, universe = uni )
      } else if(unique_gene == F) {
        y = enricher(gene_identifier, TERM2GENE = gene_set_list[[i]]
                     , pAdjustMethod = "BH"
                     , pvalueCutoff = 1
                     , qvalueCutoff = 1, universe = uni )
      }
    } else{
      if(unique_gene){
        y = enricher(unique(gene_identifier), TERM2GENE = gene_set_list[[i]]
                     , pAdjustMethod = "BH"
                     , pvalueCutoff = 1
                     , qvalueCutoff = 1 )
      } else if(unique_gene == F) {
        y = enricher(gene_identifier, TERM2GENE = gene_set_list[[i]]
                     , pAdjustMethod = "BH"
                     , pvalueCutoff = 1
                     , qvalueCutoff = 1 )
      }
    }
    if(!is.null(y)){
      y = y@result
      y$Description=fct_reorder(y$Description, as.numeric(sub("/\\d+", "", y$GeneRatio)), .desc = F)
      y$rank = rank(y$p.adjust)
    }
    enr[[i]] = y
  }
  
  enrx = mapply(function(x,y) {x$ontology=y; return(x)}, enr, names(enr),SIMPLIFY=F)
  enrx = do.call(rbind,enrx)
  
  return(enrx)
}

add_DOSE_measure_to_ORA = function(y, pattern='chr'){
  require(DOSE)
  y = mutate(  y
               , Count = as.numeric(sapply(strsplit(GeneRatio,'/'),`[`,1))
               , geneRatio = parse_ratio(GeneRatio)
               , richFactor = as.numeric(sapply(strsplit(GeneRatio,'/'),`[`,1)) / as.numeric(sub("/\\d+", "", BgRatio))
               , FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio)
  )
  y$ES = y$FoldEnrichment /y$p.adjust
  # y$ID = tolower(gsub(pattern,"", y$ID))
  y=y[order(y$rank,decreasing = F),]
  y$ID=factor(y$ID,levels=unique(y$ID))
  y
}
