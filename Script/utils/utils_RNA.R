# Tools for RNA-seq data analysis
# you MUST specify here all the package required in the script

pkgs = c( 'tidyverse'
          ,'data.table'
          ,'reshape2'
          ,'clusterProfiler'
          ,'RColorBrewer'
          , 'DESeq2'
          , 'msigdbr'
          , 'dplyr'
          , 'BiocParallel'
          , 'ggpp'
          , "ggplot2"
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



# Heatmap ====
get_heatmap3 <- function(m
                         , annotDF  = NULL
                         , annotCol = NULL
                         , fig_out  = NULL
                         , fig_h = unit(4,'cm')
                         , fig_w = unit(4,'cm')
                         , retHm    = F
                         , rowAnnot = NULL
                         , show_rownames = T
                         , myPalette = c("blue","black","red")
                         , ...){
  
  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)
  set.seed(30580)
  
  base_mean <- rowMeans(m)
  m_scaled <- t(apply(m, 1, scale))
  colnames(m_scaled) <- colnames(m)
  
  bPalette <- colorRampPalette(brewer.pal(9, "Reds"))
  ramp <- colorRamp2(c(-2, 0, 2), myPalette)
  
  hm = NULL
  if (!is.null(annotDF) & is.null(rowAnnot)) {
    if (!is.null(annotCol)) {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     col = annotCol,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     , annotation_height = unit(1, "mm")
                                     , annotation_name_side = "left"
                                     , annotation_name_gp = gpar(fontsize=8)
                                     , height = unit(10, "mm")
      )
    } else {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     , annotation_height = unit(1, "mm")
                                     , annotation_name_side = "left"
                                     , annotation_name_gp = gpar(fontsize=8)
                                     , height = unit(10, "mm")
      )
    }
    
    hm <- Heatmap(m_scaled,
                  col = ramp,
                  show_row_dend = T,
                  show_row_names = show_rownames,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize=8),
                  column_names_gp = gpar(fontsize=8),
                  column_title_gp = gpar(fontsize=10, fontface="bold"),
                  heatmap_legend_param = list(title = "row Z-score",
                                              title_gp = gpar(fontsize=8),
                                              title_position = "topcenter",
                                              # legend_width  = unit(4, "cm"),
                                              # legend_height = unit(0.5, "mm"),
                                              values_gp     = gpar(fontsize=8),
                                              legend_direction = "vertical")
                  , top_annotation = ha_column
                  , ...)
    
    
  }else if (!is.null(annotDF) & !is.null(rowAnnot)) {
    if (!is.null(annotCol)) {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     col = annotCol,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     , annotation_height = unit(1, "mm")
                                     , annotation_name_side = "left"
                                     , annotation_name_gp = gpar(fontsize=8)
                                     , height = unit(10, "mm")
      )
    } else {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     , annotation_height = unit(1, "mm")
                                     , annotation_name_side = "left"
                                     , annotation_name_gp = gpar(fontsize=8)
                                     , height = unit(10, "mm")
      )
    }
    
    hm <- Heatmap(m_scaled,
                  col = ramp,
                  show_row_dend = T,
                  show_row_names = show_rownames,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize=8),
                  column_names_gp = gpar(fontsize=8),
                  column_title_gp = gpar(fontsize=10, fontface="bold"),
                  heatmap_legend_param = list(title = "row Z-score",
                                              title_gp = gpar(fontsize=8),
                                              title_position = "topcenter",
                                              # legend_width  = unit(4, "cm"),
                                              # legend_height = unit(0.5, "mm"),
                                              values_gp     = gpar(fontsize=8),
                                              legend_direction = "vertical")
                  , top_annotation = ha_column
                  , right_annotation = rowAnnot
                  , ...)
    
    
  }  else if (is.null(annotDF) & !is.null(rowAnnot)){
    # ha_column <- new("HeatmapAnnotation")
    hm <- Heatmap(m_scaled,
                  col = ramp,
                  show_row_dend = T,
                  show_row_names = show_rownames,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize=8),
                  column_names_gp = gpar(fontsize=8),
                  column_title_gp = gpar(fontsize=10, fontface="bold"),
                  heatmap_legend_param = list(title = "row Z-score",
                                              title_gp = gpar(fontsize=8),
                                              title_position = "topcenter",
                                              # legend_width  = unit(4, "cm"),
                                              # legend_height = unit(0.5, "mm"),
                                              values_gp     = gpar(fontsize=8),
                                              legend_direction = "vertical")
                  , right_annotation = rowAnnot
                  , ...)
    
  }else {
    # ha_column <- new("HeatmapAnnotation")
    hm <- Heatmap(m_scaled,
                  col = ramp,
                  show_row_dend = T,
                  show_row_names = show_rownames,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize=8),
                  column_names_gp = gpar(fontsize=8),
                  column_title_gp = gpar(fontsize=10, fontface="bold"),
                  heatmap_legend_param = list(title = "row Z-score",
                                              title_gp = gpar(fontsize=8),
                                              title_position = "topcenter",
                                              # legend_width  = unit(4, "cm"),
                                              # legend_height = unit(0.5, "mm"),
                                              values_gp     = gpar(fontsize=8),
                                              legend_direction = "vertical")
                  , ...)
    
  }
  
  
  bmscale <- summary(base_mean)
  bmramp <- colorRamp2(c(bmscale[1],bmscale[3],bmscale[5]), bPalette(3))
  bmh <- Heatmap(base_mean
                 , name = "Mean Expression"
                 , column_names_gp = gpar(fontsize=8)
                 , show_row_names = FALSE
                 , width = unit(3, "mm")
                 , col = bmramp
                 , heatmap_legend_param = list(title = "Base Mean",title_gp = gpar(fontsize=8)))
  
  hmOut <- hm + bmh
  
  if(!is.null(fig_out)){
    pdf(file = fig_out, useDingbats = F, h=fig_h, w=fig_w, paper = "a4")
    draw(hmOut, heatmap_legend_side = "right")
    dev.off()
  } else{
    draw(hmOut, heatmap_legend_side = "right")
  }
  
  if(retHm) return(hmOut)
}

get_clusters <- function(m, hm){
  # Retrieve clusters from K-means
  clusters <- lapply(row_order(hm),
                     function(x){
                       rownames(m[x,])
                     }
  )
  return(clusters)
}



get_FC_counts = function(filenames, colnames){
  cnts = lapply(filenames, read.delim2, skip=1)
  info = cnts[[1]][,c("Geneid", "Length","gene_name", "gene_type")]
  cnts = lapply(cnts, function(x) x[,ncol(x)])
  names(cnts) = colnames
  cnts = do.call(cbind.data.frame, cnts)
  cbind.data.frame(info, cnts)
}



differential_expression_analysis <- function(count_data, sample_data, model_design,
                                             rowsums_counts_th = 1,  cores=1
                                             ,quiet = T, ...) {
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData   = sample_data,
                                design    = model_design
  )
  
  cat("[*] Model Matrix:\n" )
  print(model.matrix(model_design, sample_data))
  
  dds <- dds[ rowSums(counts(dds)) > rowsums_counts_th, ]
  
  cat("[*] Running DESeq with rowsums raw count >= ", rowsums_counts_th)
  
  if(cores>1 && !(.Platform$OS.type=="windows")){
    dds <- DESeq(dds, parallel = T, BPPARAM = BiocParallel::MulticoreParam(workers = cores), quiet = quiet, ...)
  }else if(cores>1 && .Platform$OS.type=="windows") {
    dds <- DESeq(dds, parallel = T, BPPARAM = BiocParallel::SnowParam(workers = cores), quiet = quiet, ...)
  } else{
    dds <- DESeq(dds, quiet=quiet, ...)
  }
  
  cat("\n[*] resultsNames: ", paste(resultsNames(dds), collapse=" "))
  dds
}
