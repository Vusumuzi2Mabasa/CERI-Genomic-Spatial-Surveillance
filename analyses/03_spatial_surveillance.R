world<-rnaturalearth::ne_countries(scale=110,returnclass='sf')
# These maps use longitude/latitude for display only, with no spherical areas or
# distances. Validate the planar plotting rings with GEOS, separately from s2.
invalid_geometry_count<-sum(!st_is_valid(st_set_crs(world,NA)))
if(invalid_geometry_count>0)world<-st_set_crs(st_make_valid(st_set_crs(world,NA)),4326)
check(all(st_is_valid(st_set_crs(world,NA))),'Natural Earth plotting rings valid under GEOS')
local_labs<-read_excel(input_path('local_sequencing_data.xlsx'))
flows<-read_excel(input_path('Africa_sequencing_data.xlsx'))
numeric_coordinate<-function(x)suppressWarnings(as.numeric(gsub(',','.',trimws(as.character(x)),fixed=TRUE)))
local_labs<-local_labs %>% mutate(lon=numeric_coordinate(long_submitting_lab),lat=numeric_coordinate(lat_submitting_lab),
 valid_coordinate=is.finite(lon)&is.finite(lat)&abs(lon)<=180&abs(lat)<=90&!(lon==0&lat==0))
check(nrow(local_labs)>0&nrow(flows)>0,'Sequencing-location workbooks contain records')
map_labs<-local_labs %>% filter(valid_coordinate) %>% distinct(country,submitting_lab,lon,lat)
africa<-world %>% filter(continent=='Africa')
lab_summary<-local_labs %>% group_by(country) %>% summarise(records=n(),mapped_records=sum(valid_coordinate),
 unique_submitting_labs=n_distinct(submitting_lab),.groups='drop') %>% arrange(desc(records))
write_table(lab_summary,'sequencing_location_country_summary')
p1<-ggplot()+geom_sf(data=africa,fill='#edf1f4',colour='white',linewidth=.25)+
 geom_point(data=map_labs,aes(lon,lat),colour='#126782',alpha=.75,size=2.2)+
 coord_sf(xlim=c(-26,60),ylim=c(-40,38),expand=FALSE)+theme_void()+
 labs(title='A  Submitting-laboratory locations',subtitle='Supplied CERI training workbook; overlapping locations may share a marker')
top<-lab_summary %>% slice_head(n=12)
p2<-ggplot(top,aes(reorder(country,records),records))+geom_col(fill='#126782',width=.7)+coord_flip()+
 labs(title='B  Records represented by country',subtitle='Counts of workbook rows, not sequencing throughput',x=NULL,y='Workbook records')
save_figure(p1+p2+plot_layout(widths=c(1.3,1)),'african_sequencing_locations',13,8)
canonical<-function(x){y<-tolower(stringi::stri_trans_general(enc2utf8(trimws(as.character(x))),'Latin-ASCII'));y<-gsub('[^a-z0-9]','',y);y[y %in% c('southafrica','republicofsouthafrica')]<-'southafrica';y[y %in% c('unitedkingdom','uk','greatbritain')]<-'unitedkingdom';y[y %in% c('usa','unitedstates','unitedstatesofamerica')]<-'unitedstates';y}
flows<-flows %>% mutate(origin_key=canonical(country),destination_key=canonical(SequencingCountry),
 origin_lon=numeric_coordinate(long_originating_country),origin_lat=numeric_coordinate(lat_originating_country),
 dest_lon=numeric_coordinate(long_sequencing_lab),dest_lat=numeric_coordinate(lat_sequencing_lab),
 valid_route=is.finite(origin_lon)&is.finite(origin_lat)&is.finite(dest_lon)&is.finite(dest_lat)&abs(origin_lon)<=180&abs(dest_lon)<=180&abs(origin_lat)<=90&abs(dest_lat)<=90,
 cross_border=ifelse(is.na(country)|is.na(SequencingCountry),NA,origin_key!=destination_key))
routes<-flows %>% filter(valid_route) %>% count(country,SequencingCountry,origin_lon,origin_lat,dest_lon,dest_lat,name='records')
write_table(routes,'sequencing_routes')
route_counts<-flows %>% group_by(SequencingContinent) %>% summarise(records=n(),cross_border=sum(cross_border,na.rm=TRUE),.groups='drop');write_table(route_counts,'sequencing_destination_summary')
p3<-ggplot()+geom_sf(data=world,fill='#edf1f4',colour='white',linewidth=.15)+
 geom_curve(data=routes,aes(x=origin_lon,y=origin_lat,xend=dest_lon,yend=dest_lat,linewidth=records),
 curvature=.15,colour='#126782',alpha=.35,arrow=grid::arrow(length=grid::unit(.07,'inches'),type='closed'))+
 scale_linewidth_continuous(range=c(.3,1.8),name='Workbook records')+
 coord_sf(xlim=c(-180,180),ylim=c(-60,85),expand=FALSE)+theme_void()+
 labs(title='Geographic links between origin countries and sequencing laboratories',
 subtitle='Supplied CERI training records; lines are submission links, not pathogen transmission routes')
save_figure(p3,'sequencing_submission_routes',13,7)
covid<-read.csv(input_path('covid_incidence.csv'))
check(all(c('Name','Cumulative_incidence') %in% names(covid)),'Historical mapping-table schema available')
check(!anyDuplicated(covid$Name),'Historical mapping table has unique country labels')
aliases<-c('swaziland'='eswatini','ivorycoast'='cotedivoire','republicofcongo'='republicofthecongo','congo'='republicofthecongo','democraticrepublicofthecongo'='democraticrepublicofthecongo','guineabissau'='guineabissau')
map_key<-function(x){k<-canonical(x);j<-k %in% names(aliases);k[j]<-aliases[k[j]];k}
lookup<-bind_rows(lapply(c('admin','name','name_long','sovereignt'),function(nm)data.frame(key=map_key(world[[nm]]),map_id=world$adm0_a3))) %>% distinct(key,.keep_all=TRUE)
covid$key<-map_key(covid$Name);covid<-left_join(covid,lookup,by='key')
write_table(covid %>% select(Name,Cumulative_incidence,map_id),'historical_map_country_matching')
map_values<-covid %>% filter(!is.na(map_id)) %>% group_by(map_id) %>% summarise(Cumulative_incidence=first(Cumulative_incidence),.groups='drop')
wm<-left_join(world,map_values,by=c('adm0_a3'='map_id'))
p4<-ggplot(wm)+geom_sf(aes(fill=Cumulative_incidence),colour='white',linewidth=.15)+
 scale_fill_viridis_c(option='C',trans='sqrt',na.value='#e5e7eb',name='Illustrative value')+theme_void()+
 labs(title='Historical country-level mapping example',subtitle='Handbook training values; not current or independently verified COVID-19 surveillance estimates')
save_figure(p4,'historical_country_mapping',13,7)
metrics$spatial<-list(repaired_map_geometries=invalid_geometry_count,location_records=nrow(local_labs),valid_location_records=sum(local_labs$valid_coordinate),unique_mapped_lab_locations=nrow(map_labs),location_countries=n_distinct(local_labs$country),flow_records=nrow(flows),cross_border_records=sum(flows$cross_border,na.rm=TRUE),valid_route_records=sum(flows$valid_route),historical_map_rows=nrow(covid),historical_map_matched=sum(!is.na(covid$map_id)))
