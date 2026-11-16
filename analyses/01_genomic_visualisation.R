v<-read.csv(input_path('combined_tidy_vcf.csv'))
check(all(c('sample_id','POS','DP','QUAL','REF','ALT','CHROM') %in% names(v)),'Variant schema available')
check(all(is.finite(v$DP)&v$DP>=0)&all(is.finite(v$POS)&v$POS>0),'Nonnegative depth and positive genomic positions')
check(!anyDuplicated(v[c('sample_id','CHROM','POS','REF','ALT')]),'No duplicate sample-specific allele records')
v<-v %>% mutate(retained=DP>=10 & is.finite(QUAL) & QUAL>=30)
variant_summary<-v %>% group_by(sample_id) %>% summarise(records=n(),unique_positions=n_distinct(POS),
 median_depth=median(DP),depth_IQR=IQR(DP),median_QUAL=median(QUAL,na.rm=TRUE),
 retained=sum(retained),retained_percent=100*retained/records,.groups='drop')
check(all(variant_summary$retained_percent>=0 & variant_summary$retained_percent<=100),'Variant retention percentages between zero and 100')
write_table(variant_summary,'variant_summary')
p1<-ggplot(v,aes(POS,DP,colour=sample_id))+geom_point(alpha=.6,size=1.1)+
 scale_colour_manual(values=palette_ceri)+labs(title='A  Variant-call depth across genomic positions',x='Position (bp)',y='Call depth',colour='Sample')
p2<-ggplot(v,aes(sample_id,DP,fill=sample_id))+geom_violin(alpha=.65,trim=FALSE)+geom_boxplot(width=.15,outlier.shape=NA)+
 scale_fill_manual(values=palette_ceri)+guides(fill='none')+labs(title='B  Depth distributions',x=NULL,y='Call depth')
p3<-ggplot(v,aes(sample_id,QUAL+1,fill=sample_id))+geom_boxplot(alpha=.8)+scale_y_log10()+
 scale_fill_manual(values=palette_ceri)+guides(fill='none')+labs(title='C  Variant quality',x=NULL,y='QUAL + 1 (log scale)')
p4<-ggplot(variant_summary,aes(sample_id,retained_percent,fill=sample_id))+geom_col(width=.6)+
 geom_text(aes(label=sprintf('%.1f%%',retained_percent)),vjust=-.5)+ylim(0,110)+
 scale_fill_manual(values=palette_ceri)+guides(fill='none')+labs(title='D  Illustrative filtering sensitivity',subtitle='Depth >=10 and QUAL >=30; not a validated calling rule',x=NULL,y='Records retained (%)')
save_figure((p1+p2)/(p3+p4),'genomic_variant_profiles',13,9)
metrics$variants<-list(records=nrow(v),samples=n_distinct(v$sample_id),retained=sum(v$retained),retention_percent=100*mean(v$retained))
