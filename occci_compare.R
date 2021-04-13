library(ncdf4)
library(dplyr)
library(raster)
library(lubridate)
library(maptools)
library(marmap)
library(rgeos)
library(maps)
library(mapdata)

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

library(plotrix)
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
