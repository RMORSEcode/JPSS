library(ncdf4)
library(dplyr)
library(raster)
library(lubridate)
library(maptools)
library(marmap)
library(rgeos)
library(maps)
library(mapdata)
library(plotrix)

xdt=today()
xdt=gsub('-','',xdt)

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

### read in extracted chl from SeaBASS (fluormetric and hplc)
ischl=read.csv('/home/ryan/Downloads/SeaBASS_extracted_chl (1).txt', stringsAsFactors = F)
ischl$date=as.POSIXct(ischl$datetime, tz="UTC") 


gbk=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_GBKPoly.shp")
gom=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_GOMPoly.shp")
mab=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_MABPoly.shp")
scs=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_SCSPoly.shp")
#extract just lat/lons for lines
gbk.lonlat =as.data.frame(lapply(slot(gbk, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
gom.lonlat =as.data.frame(lapply(slot(gom, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
mab.lonlat =as.data.frame(lapply(slot(mab, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
scs.lonlat =as.data.frame(lapply(slot(scs, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
# create matrix to use in in.out function from package 'mgcv'
gom.mat=as.matrix(gom.lonlat)
gbk.mat=as.matrix(gbk.lonlat)
mab.mat=as.matrix(mab.lonlat)
scs.mat=as.matrix(scs.lonlat)
m4=as.matrix(ischl[,c('longitude','latitude')]) #lon,lat from ZPD
ischl$epu=NA
ischl$epu[which(mgcv::in.out(gbk.mat, m4))]='GB'
ischl$epu[which(mgcv::in.out(gom.mat, m4))]='GOM'
ischl$epu[which(mgcv::in.out(scs.mat, m4))]='SS'
ischl$epu[which(mgcv::in.out(mab.mat, m4))]='MAB'
# test=ischl[is.na(ischl$epu),] #unassigned
GB=ischl %>% filter(epu=='GB', year(date)>1997)
table(year(GB$date))
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70"); map.axes(las=1)
points(GB$longitude, GB$latitude)
test=GB %>% filter(month(GB$date)==6)
table(year(test$date))

# nc1=nc_open('/home/ryan/Downloads/GLORYS_Bottom_Temp_2018-10-11.nc')
nc1=nc_open('/home/ryan/Downloads/global-reanalysis-phy-001-030-monthly_1608734707901.nc')
BTlist=list.files('/home/ryan/Git/NEhabitat/rasters/Fall/BT2/', pattern='RAST_NESREG_')
lon.bt=ncvar_get(nc1, 'longitude')
lat.occi=ncvar_get(nc1, 'latitude')
chl.occi=ncvar_get(nc1, 'chlor_a')
# nc.get.variable.list(nc)
attributes(nc1$var)$names
attributes(nc1$dim)$names
ncatt_get(nc1, attributes(nc1$var)$names[1])
t=ncatt_get(nc1, attributes(nc1$dim)$names[1])


### read in in situ chlorophyll-satellite matchups from NASA SEABASS
sbmodisa=read.csv('/home/ryan/Downloads/1612220222584465_chlor_a.csv', header = F, stringsAsFactors = F, skip = 29)
t=readLines('/home/ryan/Downloads/1612220222584465_chlor_a.csv', n=28)
t2=strsplit(t[27], split=',')
colnames(sbmodisa)=unlist(t2)

sbseawifs=read.csv('/home/ryan/Downloads/161222030910903_chlor_a.csv', header = F, stringsAsFactors = F, skip = 29)
t=readLines('/home/ryan/Downloads/161222030910903_chlor_a.csv', n=28)
t2=strsplit(t[27], split=',')
colnames(sbseawifs)=unlist(t2)

sbmeris=read.csv('/home/ryan/Downloads/1612220388699824_chlor_a.csv', header = F, stringsAsFactors = F, skip = 29)
t=readLines('/home/ryan/Downloads/1612220388699824_chlor_a.csv', n=28)
t2=strsplit(t[27], split=',')
colnames(sbmeris)=unlist(t2)

sbviirssnpp=read.csv('/home/ryan/Downloads/1612220443903143_chlor_a.csv', header = F, stringsAsFactors = F, skip = 29)
t=readLines('/home/ryan/Downloads/1612220443903143_chlor_a.csv', n=28)
t2=strsplit(t[27], split=',')
colnames(sbviirssnpp)=unlist(t2)

sbmodist=read.csv('/home/ryan/Downloads/1612220548231565_chlor_a.csv', header = F, stringsAsFactors = F, skip = 29)
t=readLines('/home/ryan/Downloads/1612220548231565_chlor_a.csv', n=28)
t2=strsplit(t[27], split=',')
colnames(sbmodist)=unlist(t2)

tr1=sbmodisa %>% dplyr::select( latitude, longitude, date_time, insitu_chlor_a, aqua_chlor_a)
tr2=sbmodist %>% dplyr::select( latitude, longitude, date_time, insitu_chlor_a, terra_chlor_a)
tr3=sbmeris %>% dplyr::select( latitude, longitude, date_time, insitu_chlor_a, meris_chlor_a)
tr4=sbseawifs %>% dplyr::select( latitude, longitude, date_time, insitu_chlor_a, seawifs_chlor_a)
tr5=sbviirssnpp %>% dplyr::select( latitude, longitude, date_time, insitu_chlor_a, viirs_chlor_a)

tr=tr1 %>% full_join(tr2, by=c('latitude', 'longitude', 'date_time'), name=NULL)
tr=tr%>% full_join(tr3, by=c('latitude', 'longitude', 'date_time'), name=NULL)
tr=tr%>% full_join(tr4, by=c('latitude', 'longitude', 'date_time'), name=NULL)
tr=tr%>% full_join(tr5, by=c('latitude', 'longitude', 'date_time'), name=NULL)
### merge all in-situ chlorophylls together, drop merged repeats
insitu_chl=tr %>% dplyr::select(insitu_chlor_a, insitu_chlor_a.x, insitu_chlor_a.y, insitu_chlor_a.x.x, insitu_chlor_a.y.y) %>%
  mutate(ischl=rowMeans(., na.rm=T))
tr$insitu_chl_final=insitu_chl$ischl
tr=tr %>% dplyr::select(-insitu_chlor_a, -insitu_chlor_a.x, -insitu_chlor_a.y, -insitu_chlor_a.x.x, -insitu_chlor_a.y.y)
tr$date=as.POSIXct(tr$date_time, tz="UTC") #origin = "1970-01-01")

### read in  satellite and ship matchups from Kim (daily)
satshpv4=read.csv('/home/ryan/Downloads/SATSHIP_L2-JPSS_NEC_SEABASS-CHLOR_A-OCI-OCCCI.CSV', sep = ',', 
                  header = F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Downloads/SATSHIP_L2-JPSS_NEC_SEABASS-CHLOR_A-OCI-OCCCI.CSV', n=1)
t2=strsplit(t, split=',')
colnames(satshpv4)=unlist(t2)
satshpv4$SHIP_DATE=format(satshpv4$SHIP_DATE, scientific = F)
satshpv4$SHIP_DATE=ymd_hms(satshpv4$SHIP_DATE)
satshpv5=read.csv('/home/ryan/Downloads/SATSHIP_L2-JPSS_NEC_SEABASS-CHLOR_A-CCI-OCCCI_V5.0.CSV', 
                  sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Downloads/SATSHIP_L2-JPSS_NEC_SEABASS-CHLOR_A-CCI-OCCCI_V5.0.CSV', n=1)
t2=strsplit(t, split=',')
colnames(satshpv5)=unlist(t2)
satshpv5$SHIP_DATE=format(satshpv5$SHIP_DATE, scientific = F)
satshpv5$SHIP_DATE=ymd_hms(satshpv5$SHIP_DATE)
## separate out 3x3 pixels around center chl, take geometric mean (or median)
### For OCCCI v5 data extractions
satv5chl=satshpv5 %>% tidyr::separate(SAT_CHLOR_A_CCI, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satv5chl$med=apply(satv5chl[,1:9], 1, median, na.rm=T)
satshpv5$SAT_MED_CHLOR=satv5chl$med
tv5=satshpv5[complete.cases(satshpv5$SAT_MED_CHLOR),]
tv5$date=as.POSIXct(tv5$SHIP_DATE, tz="UTC") #origin = "1970-01-01")
### now do for v4 OCCCI data
satv4chl=satshpv4 %>% tidyr::separate(SAT_CHLOR_A_OCI, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satv4chl$med=apply(satv4chl[,1:9], 1, median, na.rm=T)
satshpv4$SAT_MED_CHLOR=satv4chl$med
tv4=satshpv4[complete.cases(satshpv4$SAT_MED_CHLOR),]
tv4$date=as.POSIXct(tv4$SHIP_DATE, tz="UTC") #origin = "1970-01-01")

#__________________________________________________________________________________________
# 20210427 compare newest values of total chl, then just JPSS project Chl
### newest extracted values Apr 2021 from Audreys database of fluorometric + CTD + HPLC
### OCCCI v4.2 first ###
satshpv4=read.csv('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_SEABASS_EXTRACTED_CHL-CHLOR_A-CCI-OCCCI-SHIP_CHL.csv', 
                  sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_SEABASS_EXTRACTED_CHL-CHLOR_A-CCI-OCCCI-SHIP_CHL.csv', n=1)
t2=strsplit(t, split=',')
colnames(satshpv4)=unlist(t2)
satshpv4$SHIP_DATE=format(satshpv4$SHIPDATE, scientific = F)
satshpv4$SHIP_DATE=ymd_hms(satshpv4$SHIP_DATE)
satshpv4$SAT_DATE=format(satshpv4$SATDATE, scientific = F)
satshpv4$SAT_DATE=ymd_hms(satshpv4$SAT_DATE)
## replace Inf with NA for string of 9 values for satellite 3x3 pixels, take median value
satv4chl=satshpv4 %>% tidyr::separate(SATDATA, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satv4chl[sapply(satv4chl, is.infinite)] <- NA
satv4chl$med=apply(satv4chl[,1:9], 1, median, na.rm=T)
### add back into original DF
satshpv4$SATDATA_med=satv4chl$med
### keep values with valid center pixel and at least 4 of 9 satellite cells from 3x3 area
KH4=satshpv4 #%>% filter(!(is.na(SATDATA_0)), N>3)
## replace Inf with NA for SATDATA_0 and SATDATA
KH4$SATDATA_0[KH4$SATDATA_0==Inf]=NA # change to NA from Inf
KH4$SATDATA_all=KH4$SATDATA_0 # copy center values to 'all' column
## replace missing center values with median of 3x3 pixels
KH4$SATDATA_all[is.na(KH4$SATDATA_all)]=KH4$SATDATA_med[is.na(KH4$SATDATA_all)]
# test=KH4[1:100,]
### drop bad satellite data no center or data in 3x3 window
KH4f=KH4[complete.cases(KH4$SATDATA_all),]
### deal with multiple depths and/or replicates from ship data
# tunq=KH4f %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()==1) %>% mutate(shipnum=n())
# tdup=KH4f %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()>1) %>% mutate(shipnum=n())
# KH4mn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# KH4f2=bind_rows(tunq, KH4mn) %>% ungroup() #reassign name
### select just needed cols, look for replicates, take mean
KH4final=KH4f %>% 
  dplyr::select(SHIP_LAT, SHIP_LON, SHIP_DATE, SAT_DATE, SAT_LAT_0, SAT_LON_0, N, SHIPDATA, SATDATA_0, SATDATA_med) %>% 
  group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>%
  mutate_each(funs(mean), -(1:4)) %>%
  distinct
# test2=KH4final[1:100,]

### OCCCI v5.0 ###
satshpv5=read.csv('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_SEABASS_EXTRACTED_CHL-CHLOR_A-CCI-OCCCI_V5.0-SHIP_CHL.csv', 
                  sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_SEABASS_EXTRACTED_CHL-CHLOR_A-CCI-OCCCI_V5.0-SHIP_CHL.csv', n=1)
t2=strsplit(t, split=',')
colnames(satshpv5)=unlist(t2)
satshpv5$SHIP_DATE=format(satshpv5$SHIPDATE, scientific = F)
satshpv5$SHIP_DATE=ymd_hms(satshpv5$SHIP_DATE)
satshpv5$SAT_DATE=format(satshpv5$SATDATE, scientific = F)
satshpv5$SAT_DATE=ymd_hms(satshpv5$SAT_DATE)
## replace Inf with NA for string of 9 values for satellite 3x3 pixels, take median value
satv5chl=satshpv5 %>% tidyr::separate(SATDATA, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satv5chl[sapply(satv5chl, is.infinite)] <- NA
satv5chl$med=apply(satv5chl[,1:9], 1, median, na.rm=T)
### add back into original DF
satshpv5$SATDATA_med=satv5chl$med
### keep values with valid center pixel and at least 4 of 9 satellite cells from 3x3 area
KH5=satshpv5 #%>% filter(!(is.na(SATDATA_0)), N>3)
## replace Inf with NA for SATDATA_0 and SATDATA
KH5$SATDATA_0[KH5$SATDATA_0==Inf]=NA # change to NA from Inf
KH5$SATDATA_all=KH5$SATDATA_0 # copy center values to 'all' column
## replace missing center values with median of 3x3 pixels
KH5$SATDATA_all[is.na(KH5$SATDATA_all)]=KH5$SATDATA_med[is.na(KH5$SATDATA_all)]
# test=KH5[1:100,]
### drop bad satellite data no center or data in 3x3 window
KH5f=KH5[complete.cases(KH5$SATDATA_all),]
### deal with multiple depths and/or replicates from ship data
# tunq=KH5f %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()==1) %>% mutate(shipnum=n())
# tdup=KH5f %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()>1) %>% mutate(shipnum=n())
# KH5mn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# KH5f2=bind_rows(tunq, KH5mn) %>% ungroup() #reassign name
### select just needed cols, look for replicates, take mean
KH5final=KH5f %>% 
  dplyr::select(SHIP_LAT, SHIP_LON, SHIP_DATE, SAT_DATE, SAT_LAT_0, SAT_LON_0, N, SHIPDATA, SATDATA_0, SATDATA_med) %>% 
  group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>%
  mutate_each(funs(mean), -(1:4)) %>%
  distinct
# test2=KH5final[1:100,]

### Read in Hermes extracted data (GSM and AV 20210614) ###
## Hermes GSM
satshpgsm=read.csv('/home/ryan/Downloads/SATSHIP_MATCHIP-JPSS_SEABASS_EXTRACTED_CHL-CHLOR_A-GSM-HERMES-SHIP_CHL.csv',
                   sep = ',', header=F, stringsAsFactors = F, skip=1)
satshpgsm=read.csv('/home/ryan/Downloads/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-GSM-HERMES-SHIP_CHL.csv',
                   sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Downloads/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-GSM-HERMES-SHIP_CHL.csv', n=1)
t2=strsplit(t, split=',')
colnames(satshpgsm)=unlist(t2)
satshpgsm$SHIP_DATE=format(satshpgsm$SHIPDATE, scientific = F)
satshpgsm$SHIP_DATE=ymd_hms(satshpgsm$SHIP_DATE)
satshpgsm$SAT_DATE=format(satshpgsm$SATDATE, scientific = F)
satshpgsm$SAT_DATE=ymd_hms(satshpgsm$SAT_DATE)
## replace Inf with NA for string of 9 values for satellite 3x3 pixels, take median value
satgsm=satshpgsm %>% tidyr::separate(SATDATA, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satgsm[sapply(satgsm, is.infinite)] <- NA
satgsm$med=apply(satgsm[,1:9], 1, median, na.rm=T)
### add back into original DF
satshpgsm$SATDATA_med=satgsm$med
### keep values with valid center pixel and at least 4 of 9 satellite cells from 3x3 area
KHgsm=satshpgsm #%>% filter(!(is.na(SATDATA_0)), N>3)
## replace Inf with NA for SATDATA_0 and SATDATA
KHgsm$SATDATA_0[KHgsm$SATDATA_0==Inf]=NA # change to NA from Inf
KHgsm$SATDATA_all=KHgsm$SATDATA_0 # copy center values to 'all' column
## replace missing center values with median of 3x3 pixels
KHgsm$SATDATA_all[is.na(KHgsm$SATDATA_all)]=KHgsm$SATDATA_med[is.na(KHgsm$SATDATA_all)]
### drop bad satellite data no center or data in 3x3 window
KHgsmf=KHgsm[complete.cases(KHgsm$SATDATA_all),]
### select just needed cols, look for replicates, take mean
KHgsmfinal=KHgsmf %>% 
  dplyr::select(SHIP_LAT, SHIP_LON, SHIP_DATE, SAT_DATE, SAT_LAT_0, SAT_LON_0, N, SHIPDATA, SATDATA_0, SATDATA_med) %>% 
  group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>%
  mutate_each(funs(mean), -(1:4)) %>%
  distinct

## Hermes AV
satshpav=read.csv('/home/ryan/Downloads/SATSHIP_MATCHIP-JPSS_SEABASS_EXTRACTED_CHL-CHLOR_A-AV-HERMES-SHIP_CHL.csv',
                   sep = ',', header=F, stringsAsFactors = F, skip=1)
satshpav=read.csv('/home/ryan/Downloads/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-AV-HERMES-SHIP_CHL.csv',
                   sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Downloads/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-AV-HERMES-SHIP_CHL.csv', n=1)
t2=strsplit(t, split=',')
colnames(satshpav)=unlist(t2)
satshpav$SHIP_DATE=format(satshpav$SHIPDATE, scientific = F)
satshpav$SHIP_DATE=ymd_hms(satshpav$SHIP_DATE)
satshpav$SAT_DATE=format(satshpav$SATDATE, scientific = F)
satshpav$SAT_DATE=ymd_hms(satshpav$SAT_DATE)
## replace Inf with NA for string of 9 values for satellite 3x3 pixels, take median value
satav=satshpav %>% tidyr::separate(SATDATA, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satav[sapply(satav, is.infinite)] <- NA
satav$med=apply(satav[,1:9], 1, median, na.rm=T)
### add back into original DF
satshpav$SATDATA_med=satav$med
### keep values with valid center pixel and at least 4 of 9 satellite cells from 3x3 area
KHav=satshpav #%>% filter(!(is.na(SATDATA_0)), N>3)
## replace Inf with NA for SATDATA_0 and SATDATA
KHav$SATDATA_0[KHav$SATDATA_0==Inf]=NA # change to NA from Inf
KHav$SATDATA_all=KHav$SATDATA_0 # copy center values to 'all' column
## replace missing center values with median of 3x3 pixels
KHav$SATDATA_all[is.na(KHav$SATDATA_all)]=KHav$SATDATA_med[is.na(KHav$SATDATA_all)]
### drop bad satellite data no center or data in 3x3 window
KHavf=KHav[complete.cases(KHav$SATDATA_all),]
### select just needed cols, look for replicates, take mean
KHavfinal=KHavf %>% 
  dplyr::select(SHIP_LAT, SHIP_LON, SHIP_DATE, SAT_DATE, SAT_LAT_0, SAT_LON_0, N, SHIPDATA, SATDATA_0, SATDATA_med) %>% 
  group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>%
  mutate_each(funs(mean), -(1:4)) %>%
  distinct


## satellite center pixel value or median of 3x3 pixels when no center available
# taylor.diagram(KHgsmfinal$SHIPDATA, KHgsmfinal$SATDATA_all, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Central or 3x3 Median')
# taylor.diagram(KHavfinal$SHIPDATA, KHavfinal$SATDATA_all, col='blue', add=T)
# legend(1.5,1.5,legend=c("gsm","av"),pch=19,col=c("red","blue"))

v5=ungroup(KH5final)
v4=ungroup(KH4final)
vg=ungroup(KHgsmfinal)
va=ungroup(KHavfinal)

### merge all data sets and subset to common dates between 4 versions
df2=v5 %>% left_join(dplyr::select(v4, SHIP_DATE, SATDATA_0, SATDATA_med), by="SHIP_DATE")
df3=df2 %>% left_join(dplyr::select(vg, SHIP_DATE, SATDATA_0, SATDATA_med), by="SHIP_DATE")
df4=df3 %>% left_join(dplyr::select(va, SHIP_DATE, SATDATA_0, SATDATA_med), by="SHIP_DATE")
df5=df4[complete.cases(df4),]

### REMOVE replicate values, take mean of replicates
tunq=df5 %>% group_by(SHIP_LAT, SHIP_LON, SAT_DATE) %>% filter(n()==1) %>% mutate(num=n())
tdup=df5 %>% group_by(SHIP_LAT, SHIP_LON, SAT_DATE) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
df5x=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name

pdf(file=paste(wd,'OCCCI_satship_validation_compare_OCCCI_Hermes_2', xdt, '_', '.pdf', sep=''), height=4, width=6)
#plot sat center values
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_0.x, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Central value')
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_0.y, col='blue', add=T)
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_0.x.x, col='green', add=T)
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_0.y.y, col='black', add=T)
legend(4,5,legend=c("OCCCI v 5.0", "OCCCI v 4.2","Hermes GSM","Hermes AV"),pch=19,col=c("red","blue", "green", "black"),bty='n')
#plot median value
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_med.x, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='3x3 Median value')
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_med.y, col='blue', add=T)
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_med.x.x, col='green', add=T)
taylor.diagram(df5x$SHIPDATA, df5x$SATDATA_med.y.y, col='black', add=T)
legend(4,5,legend=c("OCCCI v 5.0", "OCCCI v 4.2","Hermes GSM","Hermes AV"),pch=19,col=c("red","blue", "green", "black"),bty='n')

txtmat=matrix(NA,nrow=4, ncol=2)
colnames(txtmat)=c('Model', 'RMSE')
txtmat[1,1]='OCCCI v5'
txtmat[2,1]='OCCCI v4.2'
txtmat[3,1]='Hermes GSM'
txtmat[4,1]='Hermes AV'
txtmat[1,2]=round(RMSE(df5x$SATDATA_0.x, df5x$SHIPDATA), 3) #v5
txtmat[2,2]=round(RMSE(df5x$SATDATA_0.y, df5x$SHIPDATA), 3)  #v4
txtmat[3,2]=round(RMSE(df5x$SATDATA_0.x.x, df5x$SHIPDATA), 3)  #gsm
txtmat[4,2]=round(RMSE(df5x$SATDATA_0.y.y, df5x$SHIPDATA), 3)  #av
gplots::textplot(txtmat)

barplot(table(year(df5x$SHIP_DATE)), main="N years")
barplot(table(month(df5x$SHIP_DATE)), main="N months")
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(df5x$SHIP_LON, df5x$SHIP_LAT, pch=16)
text(-76, 44, pos=4, paste(dim(df5x)[1],' samples', sep=''))
dev.off()



### Sample validation ###
### plot sample locations for sat-ship matches
## v4
wd='/home/ryan/Git/JPSS_2/JPSS/'
pdf(file=paste(wd,'OCCCI_satship_validation_compare_v42_v50_', xdt, '_', '.pdf', sep=''), height=4, width=6)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(KH4final$SHIP_LON, KH4final$SHIP_LAT, pch=16)
text(-76, 44, pos=4, 'v4.2 samples')
## v5
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(KH5final$SHIP_LON, KH5final$SHIP_LAT, pch=16)
text(-76, 44, pos=4, 'v5.0 samples')
### barplot of year and month coverage
barplot(table(year(KH4final$SHIP_DATE)), main="v4 year matchups")
barplot(table(year(KH5final$SHIP_DATE)), main="v5 year matchups")
barplot(table(month(KH4final$SHIP_DATE)), main="v4 month matchups")
barplot(table(month(KH5final$SHIP_DATE)), main="v5 month matchups")
### plot log in situ vs satellite values
plot(log(KH5final$SATDATA_all) ~log(KH5final$SHIPDATA), type='p', main="OCCCI v5 median vs in situ chl", xlim=c(-4,4), ylim=c(-4,4)); abline(a=0, b=1)
a=RMSE(o=log(KH5final$SHIPDATA) , m=log(KH5final$SATDATA_all))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(KH5final)[1],sep=''))
plot(log(KH4final$SATDATA_all) ~log(KH4final$SHIPDATA), type='p', main="OCCCI v4.2 median vs in situ chl", xlim=c(-4,4), ylim=c(-4,4)); abline(a=0, b=1)
a=RMSE(o=log(KH4final$SHIPDATA) , m=log(KH4final$SATDATA_all))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(KH4final)[1],sep=''))
### Taylor Diagrams ###
## satellite center pixel value or median of 3x3 pixels when no center available
taylor.diagram(KH5final$SHIPDATA, KH5final$SATDATA_all, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Central or 3x3 Median')
taylor.diagram(KH4final$SHIPDATA, KH4final$SATDATA_all, col='blue', add=T)
legend(6,7,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
## vs center pixel only
taylor.diagram(KH5final$SHIPDATA, KH5final$SATDATA_0, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Center Only')
taylor.diagram(KH4final$SHIPDATA, KH4final$SATDATA_0, col='blue', add=T)
legend(6,7,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
## log scaled center or median of 3x3
taylor.diagram(log(KH5final$SHIPDATA), log(KH5final$SATDATA_all), col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Log Central or 3x3 Median')
taylor.diagram(log(KH4final$SHIPDATA), log(KH4final$SATDATA_all), col='blue', add=T)
legend(1.5,1.5,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
## log scaled center only
taylor.diagram(log(KH5final$SHIPDATA), log(KH5final$SATDATA_0), col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Log Center Only')
taylor.diagram(log(KH4final$SHIPDATA), log(KH4final$SATDATA_0), col='blue', add=T)
legend(1.5,1.5,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
dev.off()

#### NOW JUST COMPARE JPSS DATA ####
### newest extracted values Apr 2021 from Audreys database of fluorometric + CTD + HPLC
satshpv4eco=read.csv('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-CCI-OCCCI-SHIP_CHL.csv', 
                     sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-CCI-OCCCI-SHIP_CHL.csv', n=1)
t2=strsplit(t, split=',')
colnames(satshpv4eco)=unlist(t2)
satshpv4eco$SHIP_DATE=format(satshpv4eco$SHIPDATE, scientific = F)
satshpv4eco$SHIP_DATE=ymd_hms(satshpv4eco$SHIP_DATE)
satshpv4eco$SAT_DATE=format(satshpv4eco$SATDATE, scientific = F)
satshpv4eco$SAT_DATE=ymd_hms(satshpv4eco$SAT_DATE)
## replace Inf with NA for string of 9 values for satellite 3x3 pixels, take median value
satv4ecochl=satshpv4eco %>% tidyr::separate(SATDATA, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satv4ecochl[sapply(satv4ecochl, is.infinite)] <- NA
satv4ecochl$med=apply(satv4ecochl[,1:9], 1, median, na.rm=T)
### add back into original DF
satshpv4eco$SATDATA_med=satv4ecochl$med
### keep values with valid center pixel and at least 4eco of 9 satellite cells from 3x3 area
KH4eco=satshpv4eco #%>% filter(!(is.na(SATDATA_0)), N>3)
## replace Inf with NA for SATDATA_0 and SATDATA
KH4eco$SATDATA_0[KH4eco$SATDATA_0==Inf]=NA # change to NA from Inf
KH4eco$SATDATA_all=KH4eco$SATDATA_0 # copy center values to 'all' column
## replace missing center values with median of 3x3 pixels
KH4eco$SATDATA_all[is.na(KH4eco$SATDATA_all)]=KH4eco$SATDATA_med[is.na(KH4eco$SATDATA_all)]
test=KH4eco[1:100,]
### drop bad satellite data no center or data in 3x3 window
KH4ecof=KH4eco[complete.cases(KH4eco$SATDATA_all),]
### deal with multiple depths and/or replicates from ship data
# tunq=KH4ecof %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()==1) %>% mutate(shipnum=n())
# tdup=KH4ecof %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()>1) %>% mutate(shipnum=n())
# KH4ecomn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# KH4ecof2=bind_rows(tunq, KH4ecomn) %>% ungroup() #reassign name
### select just needed cols, look for replicates, take mean
KH4ecofinal=KH4ecof %>% 
  dplyr::select(SHIP_LAT, SHIP_LON, SHIP_DATE, SAT_DATE, SAT_LAT_0, SAT_LON_0, N, SHIPDATA, SATDATA_0, SATDATA_all) %>% 
  group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>%
  mutate_each(funs(mean), -(1:4)) %>%
  distinct
# test2=KH4ecofinal[1:100,]
### newest extracted values Apr 2021 from Audreys database of fluorometric + CTD + HPLC
satshpv5eco=read.csv('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-CCI-OCCCI_V5.0-SHIP_CHL.csv', 
                     sep = ',', header=F, stringsAsFactors = F, skip=1)
t=readLines('/home/ryan/Git/JPSS_2/JPSS/2021 Extracted CHL/SATSHIP_MATCHIP-JPSS_ECOMON_EXTRACTED-CHLOR_A-CCI-OCCCI_V5.0-SHIP_CHL.csv', n=1)
t2=strsplit(t, split=',')
colnames(satshpv5eco)=unlist(t2)
satshpv5eco$SHIP_DATE=format(satshpv5eco$SHIPDATE, scientific = F)
satshpv5eco$SHIP_DATE=ymd_hms(satshpv5eco$SHIP_DATE)
satshpv5eco$SAT_DATE=format(satshpv5eco$SATDATE, scientific = F)
satshpv5eco$SAT_DATE=ymd_hms(satshpv5eco$SAT_DATE)
## replace Inf with NA for string of 9 values for satellite 3x3 pixels, take median value
satv5ecochl=satshpv5eco %>% tidyr::separate(SATDATA, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% dplyr::select(v1:v9)
satv5ecochl[sapply(satv5ecochl, is.infinite)] <- NA
satv5ecochl$med=apply(satv5ecochl[,1:9], 1, median, na.rm=T)
### add back into original DF
satshpv5eco$SATDATA_med=satv5ecochl$med
### keep values with valid center pixel and at least 5eco of 9 satellite cells from 3x3 area
KH5eco=satshpv5eco #%>% filter(!(is.na(SATDATA_0)), N>3)
## replace Inf with NA for SATDATA_0 and SATDATA
KH5eco$SATDATA_0[KH5eco$SATDATA_0==Inf]=NA # change to NA from Inf
KH5eco$SATDATA_all=KH5eco$SATDATA_0 # copy center values to 'all' column
## replace missing center values with median of 3x3 pixels
KH5eco$SATDATA_all[is.na(KH5eco$SATDATA_all)]=KH5eco$SATDATA_med[is.na(KH5eco$SATDATA_all)]
# test=KH5eco[1:100,]
### drop bad satellite data no center or data in 3x3 window
KH5ecof=KH5eco[complete.cases(KH5eco$SATDATA_all),]
### deal with multiple depths and/or replicates from ship data
# tunq=KH5ecof %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()==1) %>% mutate(shipnum=n())
# tdup=KH5ecof %>% group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>% filter(n()>1) %>% mutate(shipnum=n())
# KH5ecomn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# KH5ecof2=bind_rows(tunq, KH5ecomn) %>% ungroup() #reassign name
### select just needed cols, look for replicates, take mean
KH5ecofinal=KH5ecof %>% 
  dplyr::select(SHIP_LAT, SHIP_LON, SHIP_DATE, SAT_DATE, SAT_LAT_0, SAT_LON_0, N, SHIPDATA, SATDATA_0, SATDATA_all) %>% 
  group_by(SHIP_LAT, SHIP_LON, SHIP_DATE) %>%
  mutate_each(funs(mean), -(1:4)) %>%
  distinct
# test2=KH5ecofinal[1:100,]


wd='/home/ryan/Git/JPSS_2/JPSS/'
pdf(file=paste(wd,'OCCCI_satship_validation_compare_v42_v50_JPSS_only', xdt, '_', '.pdf', sep=''), height=4, width=6)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(KH4ecofinal$SHIP_LON, KH4ecofinal$SHIP_LAT, pch=16)
text(-76, 44, pos=4, 'v4.2 samples')
## v5
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(KH5ecofinal$SHIP_LON, KH5ecofinal$SHIP_LAT, pch=16)
text(-76, 44, pos=4, 'v5.0 samples')
### barplot of year and month coverage
barplot(table(year(KH4ecofinal$SHIP_DATE)), main="v4 year matchups")
barplot(table(year(KH5ecofinal$SHIP_DATE)), main="v5 year matchups")
barplot(table(month(KH4ecofinal$SHIP_DATE)), main="v4 month matchups")
barplot(table(month(KH5ecofinal$SHIP_DATE)), main="v5 month matchups")
### plot log in situ vs satellite values
plot(log(KH5ecofinal$SATDATA_all) ~log(KH5ecofinal$SHIPDATA), type='p', main="OCCCI v5 median vs in situ chl", xlim=c(-4,4), ylim=c(-4,4)); abline(a=0, b=1)
a=RMSE(o=log(KH5ecofinal$SHIPDATA) , m=log(KH5ecofinal$SATDATA_all))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(KH5ecofinal)[1],sep=''))
plot(log(KH4ecofinal$SATDATA_all) ~log(KH4ecofinal$SHIPDATA), type='p', main="OCCCI v4.2 median vs in situ chl", xlim=c(-4,4), ylim=c(-4,4)); abline(a=0, b=1)
a=RMSE(o=log(KH4ecofinal$SHIPDATA) , m=log(KH4ecofinal$SATDATA_all))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(KH4ecofinal)[1],sep=''))
### Taylor Diagrams ###
## satellite center pixel value or median of 3x3 pixels when no center available
taylor.diagram(KH5ecofinal$SHIPDATA, KH5ecofinal$SATDATA_all, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Central or 3x3 Median')
taylor.diagram(KH4ecofinal$SHIPDATA, KH4ecofinal$SATDATA_all, col='blue', add=T)
legend(1.5,1.5,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
## vs center pixel only
taylor.diagram(KH5ecofinal$SHIPDATA, KH5ecofinal$SATDATA_0, col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Center Only')
taylor.diagram(KH4ecofinal$SHIPDATA, KH4ecofinal$SATDATA_0, col='blue', add=T)
legend(1.5,1.5,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
## log scaled center or median of 3x3
taylor.diagram(log(KH5ecofinal$SHIPDATA), log(KH5ecofinal$SATDATA_all), col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Log Central or 3x3 Median')
taylor.diagram(log(KH4ecofinal$SHIPDATA), log(KH4ecofinal$SATDATA_all), col='blue', add=T)
legend(1,1,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
## log scaled center only
taylor.diagram(log(KH5ecofinal$SHIPDATA), log(KH5ecofinal$SATDATA_0), col='red', pos.cor=T, sd.arcs=T, show.gamma = T, main='Log Center Only')
taylor.diagram(log(KH4ecofinal$SHIPDATA), log(KH4ecofinal$SATDATA_0), col='blue', add=T)
legend(1,1,legend=c("v5.0","v4.2"),pch=19,col=c("red","blue"))
dev.off()
#________________________________________________________________________________________________

### read sat-ship matchup file from Kyle turner
ktsatmat=readxl::read_excel('/home/ryan/Downloads/nes_insitu_satellite_chl_upper10m.xlsx', na="NaN", 
                            col_types=c('numeric', 'numeric', 'guess', 'numeric', 'text', 'numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))

### Turner merge v5 sat-matchups from Kim
ktest5=left_join(tv5,ktsatmat, by=c("date"="datetime", "SHIP_LAT"="lat", "SHIP_LON"="lon"))
v5full=left_join(ktest5,tr, by=c("date", "SHIP_LAT"="latitude", "SHIP_LON"="longitude"))
v5full$insituchl=apply(v5full[,c("insitu_chl_final","in_situ_chl [mg m-3]")], 1, mean, na.rm=T)
v5full=v5full[complete.cases(v5full$insituchl),]
## now use median sat from 3x3 matrix where center pixel value is missing
v5full$SAT_CENTER_CHLOR_A_CCI[v5full$SAT_CENTER_CHLOR_A_CCI==Inf]=NA # change to NA from Inf
v5full$SAT_CENTER_CHLOR_A_CCI[is.na(v5full$SAT_CENTER_CHLOR_A_CCI)]=v5full$SAT_MED_CHLOR[is.na(v5full$SAT_CENTER_CHLOR_A_CCI)]
table(year(v5full$date))
plot(v5full$SAT_MED_CHLOR ~ v5full$insituchl, type='p', main="OCCCI v5.0 median vs in situ chl", ylim=c(0,10), xlim=c(0,20)); abline(a=0, b=1)
a=RMSE(o=v5full$insituchl , m=v5full$SAT_MED_CHLOR)
text(18,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(18,2, paste('N: ', dim(v5full)[1],sep=''))

### Turner merge v4 sat-matchups from Kim
ktest4=left_join(tv4,ktsatmat, by=c("date"="datetime", "SHIP_LAT"="lat", "SHIP_LON"="lon"))
v4full=left_join(ktest4,tr, by=c("date", "SHIP_LAT"="latitude", "SHIP_LON"="longitude"))
v4full$insituchl=apply(v4full[,c("insitu_chl_final","in_situ_chl [mg m-3]")], 1, mean, na.rm=T)
v4full=v4full[complete.cases(v4full$insituchl),]
## now use median sat from 3x3 matrix where center pixel value is missing
v4full$SAT_CENTER_CHLOR_A_OCI[v4full$SAT_CENTER_CHLOR_A_OCI==Inf]=NA # change to NA from Inf
v4full$SAT_CENTER_CHLOR_A_OCI[is.na(v4full$SAT_CENTER_CHLOR_A_OCI)]=v4full$SAT_MED_CHLOR[is.na(v4full$SAT_CENTER_CHLOR_A_OCI)]
table(year(v4full$date))
plot(v4full$SAT_MED_CHLOR ~ v4full$insituchl, type='p', main="OCCCI v4.2 median vs in situ chl", ylim=c(0,10), xlim=c(0,20)); abline(a=0, b=1)
a=RMSE(o=v4full$insituchl , m=v4full$SAT_MED_CHLOR)
text(18,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(18,2, paste('N: ', dim(v4full)[1],sep=''))

save(tr, file='/home/ryan/Git/JPSS_2/JPSS/NASAcalibration.RData')
save(ktsatmat, file='/home/ryan/Git/JPSS_2/JPSS/TurnerSatelliteMatches.Rdata')

plot(log(v4full$SAT_MED_CHLOR) ~log(v4full$insituchl), type='p', main="OCCCI v4.2 median vs in situ chl", xlim=c(-4,4), ylim=c(-4,4)); abline(a=0, b=1)
a=RMSE(o=log(v4full$insituchl) , m=log(v4full$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(v4full)[1],sep=''))

plot(log(v5full$SAT_MED_CHLOR) ~log(v5full$insituchl), type='p', main="OCCCI v5.0 median vs in situ chl", xlim=c(-4,4), ylim=c(-4,4)); abline(a=0, b=1)
a=RMSE(o=log(v5full$insituchl) , m=log(v5full$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(v5full)[1],sep=''))
par(mar = c(0,0,0,0))
par(oma = c(0,0,0,0))
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON, v4full$SHIP_LAT, pch=16)

### check on high chl value locations
# [which(v4full$SAT_MED_CHLOR>3)]
## satellite vals > 2
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON[which(v4full$SAT_MED_CHLOR>2)], v4full$SHIP_LAT[which(v4full$SAT_MED_CHLOR>2)], pch=16)
text(-68,38, paste('Sat vals > 2ug/L ', '\n', 'N: ', length(v4full$SHIP_LON[which(v4full$SAT_MED_CHLOR>2)]),sep=''))
## sat vals >=1
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON[which(v4full$SAT_MED_CHLOR>=1)], v4full$SHIP_LAT[which(v4full$SAT_MED_CHLOR>=1)], pch=16)
text(-68,38, paste('Sat vals >= 1ug/L ', '\n', 'N: ', length(v4full$SHIP_LON[which(v4full$SAT_MED_CHLOR>=1)]),sep=''))
## sat center vals >=1
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON[which(v4full$SAT_CENTER_CHLOR_A_OCI >=1)], v4full$SHIP_LAT[which(v4full$SAT_CENTER_CHLOR_A_OCI>=1)], pch=16)
text(-68,38, paste('v4 Sat center vals >= 1ug/L ', '\n', 'N: ', length(v4full$SHIP_LON[which(v4full$SAT_CENTER_CHLOR_A_OCI>=1)]),sep=''))
## sat center vals >=1
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v5full$SHIP_LON[which(v5full$SAT_CENTER_CHLOR_A_CCI >=1)], v5full$SHIP_LAT[which(v5full$SAT_CENTER_CHLOR_A_CCI>=1)], pch=16)
text(-68,38, paste('v5 Sat center vals >= 1ug/L ', '\n', 'N: ', length(v5full$SHIP_LON[which(v5full$SAT_CENTER_CHLOR_A_CCI>=1)]),sep=''))


## in situ vals >2
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON[which(v4full$insituchl >2)], v4full$SHIP_LAT[which(v4full$insituchl>2)], pch=16, col='red')
text(-68,38, paste('in situ > 2ug/L ', '\n', 'N: ', length(v4full$SHIP_LON[which(v4full$insituchl>2)]),sep=''))
## in situ vals >=1
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON[which(v4full$insituchl>=1)], v4full$SHIP_LAT[which(v4full$insituchl>=1)], pch=16, col='red')
text(-68,38, paste('in situ >=1 ug/L ', '\n', 'N: ', length(v4full$SHIP_LON[which(v4full$insituchl>=1)]),sep=''))


## all points v4
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(v4full$SHIP_LON, v4full$SHIP_LAT, pch=16)
text(-68,38, paste('N: ', dim(v4full)[1],sep=''))

## median of 3x3 pixels
taylor.diagram(v4full$insituchl, v4full$SAT_MED_CHLOR, col='red', pos.cor=T, sd.arcs=T, show.gamma = T)
taylor.diagram(v5full$insituchl, v5full$SAT_MED_CHLOR, col='blue', add=T)
legend(9,10,legend=c("v4.2","v5.0"),pch=19,col=c("red","blue"))
## vs center pixel (when avail, otherwise median of 3x3 if missing)
taylor.diagram(v4full$insituchl, v4full$SAT_CENTER_CHLOR_A_OCI, col='red', pos.cor=T, sd.arcs=T, show.gamma = T)
taylor.diagram(v5full$insituchl, v5full$SAT_CENTER_CHLOR_A_CCI, col='blue', add=T)
legend(9,10,legend=c("v4.2","v5.0"),pch=19,col=c("red","blue"))
## log scaled median 3x3 value
taylor.diagram(log(v4full$insituchl), log(v4full$SAT_MED_CHLOR), col='red', pos.cor=T, sd.arcs=T, show.gamma = T)
taylor.diagram(log(v5full$insituchl), log(v5full$SAT_MED_CHLOR), col='blue', add=T)
legend(1.8,2,legend=c("v4.2","v5.0"),pch=19,col=c("red","blue"))
## log scaled center pixel value
taylor.diagram(log(v4full$insituchl), log(v4full$SAT_CENTER_CHLOR_A_OCI), col='red', pos.cor=T, sd.arcs=T, show.gamma = T)
taylor.diagram(log(v5full$insituchl), log(v5full$SAT_CENTER_CHLOR_A_CCI), col='blue', add=T)
legend(1.8,2,legend=c("v4.2","v5.0"),pch=19,col=c("red","blue"))
