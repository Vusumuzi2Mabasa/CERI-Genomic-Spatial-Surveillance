args<-commandArgs(trailingOnly=FALSE);fa<-sub('^--file=','',args[grepl('^--file=',args)])
if(length(fa)&&basename(fa)=='build_docs.R')setwd(normalizePath(file.path(dirname(fa),'..'),winslash='/'))
root<-normalizePath('.',winslash='/')
library(jsonlite)
m<-fromJSON('results/tables/run_summary.json')
read_result<-function(n)read.csv(file.path('results/tables',paste0(n,'.csv')),check.names=FALSE)
md_table<-function(d){
 for(n in names(d))if(is.numeric(d[[n]]))d[[n]]<-format(round(d[[n]],2),trim=TRUE,scientific=FALSE)
 lines<-c(paste0('| ',paste(names(d),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(d)),collapse=' | '),' |'))
 for(i in seq_len(nrow(d)))lines<-c(lines,paste0('| ',paste(d[i,],collapse=' | '),' |'))
 paste(lines,collapse='\n')
}
fmt<-function(x)format(x,trim=TRUE,scientific=FALSE)
values<-list(
 VARIANT_RECORDS=m$variants$records,VARIANT_SAMPLES=m$variants$samples,VARIANT_RETAINED=m$variants$retained,VARIANT_RETENTION=sprintf('%.1f',m$variants$retention_percent),
 TRAINING_TIPS=m$phylogeny$training_tips,TRAINING_METADATA=m$phylogeny$training_metadata,TRAINING_EXTRAS=m$phylogeny$training_extra_metadata,TRAINING_HOSTS=m$phylogeny$training_hosts,
 SHIGELLA_TIPS=m$phylogeny$shigella_tips,SHIGELLA_COUNTRIES=m$phylogeny$shigella_countries,MACROLIDE_RECORDED=m$phylogeny$macrolide_gene_recorded,
 LOCATION_RECORDS=m$spatial$location_records,LOCATION_COUNTRIES=m$spatial$location_countries,MAPPED_RECORDS=m$spatial$valid_location_records,UNIQUE_LAB_LOCATIONS=m$spatial$unique_mapped_lab_locations,
 FLOW_RECORDS=m$spatial$flow_records,CROSS_BORDER_RECORDS=m$spatial$cross_border_records,HISTORICAL_MATCHED=m$spatial$historical_map_matched,HISTORICAL_ROWS=m$spatial$historical_map_rows,
 OUTBREAK_CASES=m$outbreak$cases,ONSET_OBSERVED=m$outbreak$onset_observed,ONSET_MISSING=m$outbreak$onset_missing,DEATHS=m$outbreak$deaths,RECOVERED=m$outbreak$recovered,OUTCOME_MISSING=m$outbreak$outcome_missing,
 DEATH_PERCENT=sprintf('%.1f',m$outbreak$known_outcome_death_percent),DEATH_LOWER=sprintf('%.1f',m$outbreak$known_outcome_lower),DEATH_UPPER=sprintf('%.1f',m$outbreak$known_outcome_upper),
 PEAK_WEEK=m$outbreak$peak_week,PEAK_CASES=m$outbreak$peak_week_cases,ADMISSION_MEDIAN=m$outbreak$median_admission_delay,NEGATIVE_DELAYS=m$outbreak$negative_admission_delays,
 RECORDED_LINKS=m$outbreak$recorded_links,MATCHED_LINKS=m$outbreak$matched_endpoint_links,UNMATCHED_LINKS=m$outbreak$unmatched_endpoint_links,DISPLAYED_LINKS=m$outbreak$displayed_links,QC_COUNT=nrow(read_result('qc_checks')),
 VARIANT_TABLE=md_table(read_result('variant_summary')[,c('sample_id','records','median_depth','median_QUAL','retained','retained_percent')]),
 TREE_MATCH_TABLE=md_table(read_result('phylogeny_matching')),SHIGELLA_TABLE=md_table(read_result('shigella_continent_summary')),
 DESTINATION_TABLE=md_table(read_result('sequencing_destination_summary')),OUTCOME_TABLE=md_table(read_result('outbreak_outcomes')))
text<-paste(readLines('references/README.template',warn=FALSE,encoding='UTF-8'),collapse='\n')
for(n in names(values))text<-gsub(paste0('{{',n,'}}'),fmt(values[[n]]),text,fixed=TRUE)
stopifnot(!grepl('{{',text,fixed=TRUE))
writeLines(enc2utf8(text),'README.md',useBytes=TRUE)
if(!requireNamespace('rmarkdown',quietly=TRUE)||!requireNamespace('flexdashboard',quietly=TRUE))stop('Install reporting dependencies first.')
if(!rmarkdown::pandoc_available()){
 candidates<-c(Sys.getenv('RSTUDIO_PANDOC'),Sys.glob('C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/*'),Sys.glob('C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'))
 for(p in candidates)if(file.exists(file.path(p,'pandoc.exe'))||file.exists(file.path(p,'pandoc'))){Sys.setenv(RSTUDIO_PANDOC=p);break}
}
if(!rmarkdown::pandoc_available())stop('Pandoc is required: install RStudio or Quarto, or set RSTUDIO_PANDOC to its directory.')
rmarkdown::render('report.Rmd',output_file='surveillance_report.html',output_dir=file.path(root,'results/reports'),knit_root_dir=root,quiet=TRUE,envir=new.env(parent=globalenv()))
rmarkdown::render('dashboard/outbreak_dashboard.Rmd',output_file='outbreak_dashboard.html',output_dir=file.path(root,'results/reports'),knit_root_dir=root,params=list(project_root=root),quiet=TRUE,envir=new.env(parent=globalenv()))
cat('README, research report and dashboard rendered.\n')
