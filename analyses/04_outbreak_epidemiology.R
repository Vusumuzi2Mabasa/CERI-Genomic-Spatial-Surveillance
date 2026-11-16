linelist<-readRDS(input_path('linelist_cleaned.rds'))
check(all(c('case_id','date_onset','date_hospitalisation','outcome','age_years','infector') %in% names(linelist)),'Fictional outbreak linelist schema available')
check(!anyDuplicated(linelist$case_id),'Fictional case identifiers unique')
check(all(is.na(linelist$age_years)|linelist$age_years>=0),'Recorded ages nonnegative')
missing_table<-data.frame(variable=names(linelist),missing=vapply(linelist,function(x)sum(is.na(x)),integer(1))) %>% mutate(missing_percent=100*missing/nrow(linelist))
write_table(missing_table,'outbreak_missingness')
dated<-linelist %>% filter(!is.na(date_onset)) %>% mutate(week=as.Date(cut(as.Date(date_onset),'week',start.on.monday=TRUE)))
weekly<-dated %>% count(week,name='cases') %>% complete(week=seq(min(week),max(week),by='week'),fill=list(cases=0))
check(sum(weekly$cases)==sum(!is.na(linelist$date_onset)),'Weekly case totals reconcile with available onset dates')
write_table(weekly,'weekly_onset_cases')
outcomes<-linelist %>% mutate(outcome=label_missing(outcome)) %>% count(outcome,name='cases') %>% mutate(percent=100*cases/nrow(linelist));write_table(outcomes,'outbreak_outcomes')
hospital_summary<-linelist %>% group_by(hospital) %>% summarise(cases=n(),deaths=sum(outcome=='Death',na.rm=TRUE),recovered=sum(outcome=='Recover',na.rm=TRUE),missing_outcome=sum(is.na(outcome)),.groups='drop')
ci<-t(vapply(seq_len(nrow(hospital_summary)),function(i)wilson(hospital_summary$deaths[i],hospital_summary$deaths[i]+hospital_summary$recovered[i]),numeric(3)))
hospital_summary<-cbind(hospital_summary,as.data.frame(ci));write_table(hospital_summary,'hospital_outcome_summary')
delay<-as.numeric(as.Date(linelist$date_hospitalisation)-as.Date(linelist$date_onset))
valid_delay<-is.finite(delay)&delay>=0
write_table(data.frame(total_cases=nrow(linelist),missing_onset=sum(is.na(linelist$date_onset)),negative_admission_delays=sum(is.finite(delay)&delay<0),valid_admission_delays=sum(valid_delay),median_valid_delay=median(delay[valid_delay]),IQR_valid_delay=IQR(delay[valid_delay])),'outbreak_date_qc')
p1<-ggplot(weekly,aes(week,cases))+geom_col(fill='#126782',width=6)+
 scale_x_date(date_labels='%b %Y',date_breaks='2 months')+
 labs(title='A  Weekly symptom-onset curve',subtitle='Fictional Ebola outbreak; missing onset dates excluded',x=NULL,y='Cases')
p2<-ggplot(outcomes,aes(reorder(outcome,cases),cases,fill=outcome))+geom_col(show.legend=FALSE)+coord_flip()+
 scale_fill_manual(values=c(Death='#9b2226','Not recorded'='#919ba3',Recover='#4c956c'))+
 labs(title='B  Recorded outcomes',x=NULL,y='Cases')
p3<-ggplot(linelist %>% filter(!is.na(age_years)),aes(age_years))+geom_histogram(binwidth=5,boundary=0,fill='#ee9b00',colour='white')+
 labs(title='C  Age distribution',x='Age (years)',y='Cases')
p4<-ggplot(data.frame(delay=delay[valid_delay]),aes(delay))+geom_histogram(binwidth=1,boundary=-.5,fill='#4c956c',colour='white')+
 labs(title='D  Onset-to-admission interval',subtitle='Nonnegative, observed intervals only',x='Days',y='Cases')
save_figure((p1+p2)/(p3+p4),'fictional_outbreak_overview',13,9)
links_raw<-linelist %>% filter(!is.na(infector)) %>% transmute(from=as.character(infector),to=as.character(case_id))
endpoint_match<-links_raw$from %in% linelist$case_id & links_raw$to %in% linelist$case_id
links<-links_raw[endpoint_match,]
write_table(data.frame(recorded_links=nrow(links_raw),complete_endpoint_links=nrow(links),unmatched_endpoint_links=sum(!endpoint_match)),'outbreak_link_qc')
check(nrow(links)+sum(!endpoint_match)==nrow(links_raw),'Transmission endpoint audit reconciles all recorded links')
check(all(links$from %in% linelist$case_id)&all(links$to %in% linelist$case_id),'Displayed transmission links have matched endpoints')
check(!any(links$from==links$to),'No self-links in supplied transmission contacts')
nodes<-linelist %>% filter(!is.na(date_onset),date_onset>=as.Date('2014-06-01'),date_onset<=as.Date('2014-06-30')) %>% arrange(date_onset,case_id) %>% slice_head(n=80) %>% mutate(y=row_number(),case_id=as.character(case_id),outcome=label_missing(outcome))
edges<-links %>% filter(from %in% nodes$case_id,to %in% nodes$case_id) %>%
 left_join(nodes %>% select(from=case_id,x=date_onset,y),by='from') %>%
 left_join(nodes %>% select(to=case_id,xend=date_onset,yend=y),by='to')
p5<-ggplot()+geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),colour='#8a9ba8',alpha=.65,arrow=grid::arrow(length=grid::unit(.08,'inches')))+
 geom_point(data=nodes,aes(date_onset,y,colour=outcome),size=2.2)+
 scale_colour_manual(values=c(Death='#9b2226','Not recorded'='#919ba3',Recover='#4c956c'))+
 labs(title='Recorded links in a fictional outbreak subset',subtitle='First 80 dated cases in June 2014; only links internal to that subset are displayed',x='Symptom-onset date',y='Case display order',colour='Outcome')+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.major.y=element_blank())
save_figure(p5,'fictional_transmission_links',12,8)
n_death<-sum(linelist$outcome=='Death',na.rm=TRUE);n_recover<-sum(linelist$outcome=='Recover',na.rm=TRUE);overall_ci<-wilson(n_death,n_death+n_recover)
metrics$outbreak<-list(cases=nrow(linelist),onset_observed=nrow(dated),onset_missing=sum(is.na(linelist$date_onset)),deaths=n_death,recovered=n_recover,outcome_missing=sum(is.na(linelist$outcome)),known_outcome_death_percent=100*overall_ci[1],known_outcome_lower=100*overall_ci[2],known_outcome_upper=100*overall_ci[3],peak_week=as.character(weekly$week[which.max(weekly$cases)]),peak_week_cases=max(weekly$cases),median_admission_delay=median(delay[valid_delay]),negative_admission_delays=sum(is.finite(delay)&delay<0),recorded_links=nrow(links_raw),matched_endpoint_links=nrow(links),unmatched_endpoint_links=sum(!endpoint_match),displayed_link_nodes=nrow(nodes),displayed_links=nrow(edges))
