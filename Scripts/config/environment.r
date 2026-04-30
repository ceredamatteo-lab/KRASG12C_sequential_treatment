## ***************************
##
## Script name: environment.R
## Purpose of script: config file to set environment variables used in tbe oproject 
## Author: mc
## Date Created: 
##
## ***************************
##
## Notes:
##   
## ***************************




if( Sys.info()['nodename']=='cgb01'){
  
  
}else if( Sys.info()['nodename']=='psychobook-2.local'){
  
  
}else if( Sys.info()['user'] == 'cgb01'){
  DATA_DIR='/Users/cgb01/repo/KRAS_WGS/'
  CGB_DIR='/Users/cgb01/'
  GIT_LOCAL_DIR='/Users/cgb01/repo/'
  CGB_SHARED='~/Dropbox (HuGeF)/'
  BEDTOOLS='/Users/cgb01/miniconda/bin/bedtools'
  
}else if( Sys.info()['user'] =='mariachiara.grieco' ){
  DATA_DIR="/hpcnfs/data/cgb/Ambrogio/WGS/" 
  CGB_DIR='/hpcnfs/data/cgb/'
  GIT_LOCAL_DIR='/hpcnfs/home/mariachiara.grieco/'
  
  
}else if( Sys.info()['user'] =='serena.peirone' ){
  DATA_DIR="/hpcnfs/data/cgb/Ambrogio/WGS/" 
  CGB_DIR='/hpcnfs/data/cgb/'
  GIT_LOCAL_DIR='/hpcnfs/home/serena.peirone/'
  
  
}else if(Sys.info()['user'] =="jacop"){
  DATA_DIR="/Users/jacop/Laura/IFOM_analysis/genomics_laura/WGS/"
  WORK_DIR="/Users/jacop/Laura/IFOM_analysis/genomics_laura/WGS/"
  CGB_DIR='/Users/jacop/Laura/CGBlab/'
  GIT_LOCAL_DIR='/Users/jacop/Laura/CGBlab/repo/'
  
}else{
}

# HELPER functions --------

read_csv_paste = function(path, file, ...){
  require(tidyverse)
  read_csv(file=paste0(path, file), ...)
}

pdf_paste = function(path, file, ...){
  pdf(file=paste0(path, file), ...)
}