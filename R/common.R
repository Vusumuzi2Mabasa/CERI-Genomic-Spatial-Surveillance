suppressPackageStartupMessages({
  library(ggplot2);library(dplyr);library(tidyr);library(readxl)
  library(ape);library(sf);library(patchwork);library(jsonlite)
})
options(stringsAsFactors=FALSE)
palette_ceri <- c('#126782','#ee9b00','#9b2226','#4c956c','#6d597a','#5c677d')
theme_set(theme_minimal(base_size=12)+theme(
  plot.title=element_text(face='bold',colour='#172b4d',size=15),
  plot.subtitle=element_text(colour='#52616b',size=10),
  panel.grid.minor=element_blank(),legend.position='bottom',
  axis.title=element_text(face='bold'),plot.background=element_rect(fill='white',colour=NA)))
save_figure <- function(plot,name,width=12,height=8){
  ggsave(file.path('results/figures',paste0(name,'.png')),plot,width=width,height=height,dpi=240,bg='white')
  ggsave(file.path('results/figures',paste0(name,'.pdf')),plot,width=width,height=height,bg='white')
}
write_table <- function(x,name)write.csv(as.data.frame(x),file.path('results/tables',paste0(name,'.csv')),row.names=FALSE,na='')
check <- function(condition,description){
  qc_records[[length(qc_records)+1L]] <<- data.frame(check=description,passed=isTRUE(condition))
  if(!isTRUE(condition))stop('QC failed: ',description,call.=FALSE)
}
wilson <- function(x,n,z=qnorm(.975)){
  ifelse(n==0,NA_real_,0) # zero denominators are handled explicitly below
  if(n==0)return(c(estimate=NA,lower=NA,upper=NA))
  p<-x/n;den<-1+z^2/n;c(estimate=p,lower=(p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/den,
    upper=(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/den)
}
input_path<-function(name)file.path('data/example',name)
label_missing<-function(x){x<-as.character(x);x[is.na(x)|!nzchar(trimws(x))]<-'Not recorded';x}
read_tree_metadata<-function(path){
 m<-read.csv(path,fileEncoding='UTF-8-BOM',check.names=FALSE)
 names(m)[1]<-'name';m
}
