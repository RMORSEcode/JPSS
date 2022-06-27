library (dplyr)
library(lubridate)
library(raster)
library(ncdf4)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
# d <- loadRData("~/blah/ricardo.RData")

## geometric mean for extractions
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## log normalized mean for extractions
log_mean = function(x, na.rm=TRUE){
  x2=log10(x)
  x3=mean(x2, na.rm=na.rm)
  x4=10^x3
  return(x4)
}

### Function to extract data using a shapefile
extract_calc_mean=function(x, shp){
  v2=list()
  for(i in 1:dim(x)[3]){
    v=extract(x[[i]], shp)
    v1=lapply(v, function(xx) mean(xx, na.rm=T))
    v2[i]=list(v1)
  }
  m=matrix(unlist(v2), ncol=dim(x)[3], nrow=length(shp@polygons)) # box 0-29 =rows, years =cols
  # colnames(m)=seq(1998, 2016, by=1)
  # rownamse(m)=
  return(m)
}


### Function to extract data using a shapefile, mean or median, na.option is T or F for na.rm call
extract_calc=function(x, shp, fun, na.option){
  if (!exists(fun)){
    print('No function supplied')
    break
  }
  v2=list()
  if (fun=='mean'){
    for(i in 1:dim(x)[3]){
      v=extract(x[[i]], shp)
      v1=lapply(v, function(xx) mean(xx, na.rm=na.option))
      v2[i]=list(v1)
    }
  }
  else if (fun=='median'){
    for(i in 1:dim(x)[3]){
      v=extract(x[[i]], shp)
      v1=lapply(v, function(xx) median(xx, na.rm=na.option))
      v2[i]=list(v1)
    }
  }
  m=matrix(unlist(v2), ncol=dim(x)[3], nrow=length(shp@polygons)) # box 0-29 =rows, years =cols
  return(m)
}
### Function to extract data using a shapefile with geometric mean for Chl
extract_calc_geo=function(x, shp){
  v2=list()
  for(i in 1:dim(x)[3]){
    v=extract(x[[i]], shp)
    v1=lapply(v, function(xx) gm_mean(xx))
    v2[i]=list(v1)
  }
  m=matrix(unlist(v2), ncol=dim(x)[3], nrow=length(shp@polygons)) # box 0-29 =rows, years =cols
  # colnames(m)=seq(1998, 2016, by=1)
  # rownamse(m)=
  return(m)
}

### Function to extract data using a shapefile after log normalize data and take mean for Chl
extract_calc_logmean=function(x, shp){
  v2=list()
  for(i in 1:dim(x)[3]){
    v=extract(x[[i]], shp)
    v1=lapply(v, function(xx) log_mean(xx))
    v2[i]=list(v1)
  }
  m=matrix(unlist(v2), ncol=dim(x)[3], nrow=length(shp@polygons)) # box 0-29 =rows, years =cols
  # colnames(m)=seq(1998, 2016, by=1)
  # rownamse(m)=
  return(m)
}

### Function to load ncdf files as stacked raster
nc2raster=function(x, varn){
  s=stack()
  yy=length(x)
  for (i in 1:length(x)){
    r <- raster(x[i],  varname = varn)
    s=stack(s, r)
    print(paste(i,' of ', yy, ' files', sep=''))
  }
  return(s)
}



### get shapefile for subset
setwd("C:/Users/ryan.morse/Desktop/NES_5area")
nes.five=rgdal::readOGR('nes_gbk_gome_gomw_mabn_mabsPoly.shp')
nes.five=rgdal::readOGR('/home/ryan/Desktop/shapefiles/NES_5_area/nes_gbk_gome_gomw_mabn_mabsPoly.shp')
### Read in OCCI data
# nc1=nc_open('C:/Users/ryan.morse/Documents/GitHub/JPSS/CCI_ALL-v4.2-8DAY.nc') #udpated with data fix 2019
nc1=nc_open('/home/ryan/Git/JPSS/CCI_ALL-v4.2-8DAY.nc') # v4.2 udpated with data fix 2019
nc1=nc_open('/home/ryan/Git/JPSS/CCI_ALL-v5.0-8DAY (3).nc') # v5.0 20201220 downloaded
### get shapefile for subset
setwd("C:/Users/ryan.morse/Desktop/NES_5area")
nes.five=rgdal::readOGR('nes_gbk_gome_gomw_mabn_mabsPoly.shp')

### Read in OCCI data
nc1=nc_open('C:/Users/ryan.morse/Documents/GitHub/JPSS/CCI_ALL-v4.2-8DAY.nc') #udpated with data fix 2019
lon.occi=ncvar_get(nc1, 'lon')
lat.occi=ncvar_get(nc1, 'lat')
chl.occi=ncvar_get(nc1, 'chlor_a')
lg.chl.occi=log10(chl.occi)
chl.occi.bias=ncvar_get(nc1, 'chlor_a_log10_bias') # only if using new file in downloads (2GB)
chl.occi.rmsd=ncvar_get(nc1, 'chlor_a_log10_rmsd') # only if using new file in downloads (2GB)
time.occi=ncvar_get(nc1, 'time') #days since Jan 1, 1970
dim(chl.occi)
colnames(chl.occi)=lat.occi
rownames(chl.occi)=lon.occi
nc_close(nc1)

# test=month.day.year(time.occi, c(1,1,1970)) # these are not 8-days apart.... something odd -see manual
test=as.Date(time.occi, origin="1970-01-01")
occi.date=data.frame(test)
# occi.date$F1=paste(occi.date[,3], occi.date[,1], occi.date[,2], sep='-')
# occi.date$DOY=as.numeric(strftime(occi.date$F1, '%j'))
# ddiff=diff(occi.date$DOY)
# occi.date$diff=c(0, ddiff) # see OC-CCI manual in JPSS/calibration folder for list of missing dates, explains why some are not 8d

# table(occi.date$year)
# occi.date$week=c(seq(from=31, to=46, by=1),  rep(seq(from=1, to=46, by=1), 21))
test=month.day.year(time.occi, c(1,1,1970)) # these are not 8-days apart.... something odd -see manual
occi.date=data.frame(test)
occi.date$F1=paste(occi.date[,3], occi.date[,1], occi.date[,2], sep='-')
occi.date$DOY=as.numeric(strftime(occi.date$F1, '%j'))
ddiff=diff(occi.date$DOY)
occi.date$diff=c(0, ddiff) # see OC-CCI manual in JPSS/calibration folder for list of missing dates, explains why some are not 8d

table(occi.date$year)
occi.date$week=c(seq(from=31, to=46, by=1),  rep(seq(from=1, to=46, by=1), 21))


## Chl_loop over and stack rasters
bb=c(-80, -60, 32, 48)
bb=c(-77, -64, 35, 45)


m2=t(chl.occi[,,1])
# m2=t(m)#[ncol(m):1,] # flip and transpose matrix
occi=raster(m2)
extent(occi)=bb
crs(occi)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
for(i in 2:dim(chl.occi)[3]){
  m2=t(chl.occi[,,i])
  # m2=t(m)#[ncol(m):1,] # flip and transpose matrix
  xx=raster(m2)
  extent(xx)=bb
  crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  # plot(xx)
  occi=stack(occi, xx)
  print(i)
}

save(occi, file='/home/ryan/Git/JPSS/OCCI_v5_8day_stacked_rasters.rda')
save(test, file='/home/ryan/Git/JPSS/OCCI_v5_8day_stacked_raster_dates.rda')


### create raster stacks of data 
chl.av.8d=nc2raster(files.chl1.av[,1], 'CHL1_mean')
chl.gsm.8d=nc2raster(files.chl1.gsm[,1], 'CHL1_mean')
chl.oc5.8d=nc2raster(files.oc5[,1], 'CHL-OC5_mean')

# create time series by box for ecomon strata
# m.av=extract_calc(chl.av.8d[[1:1005]], nes.five)
w.gsm.nes=extract_calc(chl.gsm.8d[[1:1005]], nes.five)
# m.oc5=extract_calc(chl.oc5[[1:1005]], nes.five)
kf.gsm=dates.8d %>% select(M, Y, DOY1, week, DOY2)
kf.gsm$nes=w.gsm.nes[1,]
kf.gsm$gbk=w.gsm.nes[2,]
kf.gsm$gome=w.gsm.nes[3,]
kf.gsm$gomw=w.gsm.nes[4,]
kf.gsm$mabn=w.gsm.nes[5,]
kf.gsm$mabs=w.gsm.nes[6,]
save(kf.gsm, file="hermes_8d_extracted_chl.RData")

### log chl prior to aggregation, then take exponent of aggregate
lg.gsm.8d=log(chl.gsm.8d)
lg.gsm.nes=extract_calc(lg.gsm.8d[[1:1005]], nes.five)
kf.lg.gsm=dates.8d %>% select(M, Y, DOY1, week, DOY2)
kf.lg.gsm$nes=exp(lg.gsm.nes[1,])
kf.lg.gsm$gbk=exp(lg.gsm.nes[2,])
kf.lg.gsm$gome=exp(lg.gsm.nes[3,])
kf.lg.gsm$gomw=exp(lg.gsm.nes[4,])
kf.lg.gsm$mabn=exp(lg.gsm.nes[5,])
kf.lg.gsm$mabs=exp(lg.gsm.nes[6,])
save(kf.lg.gsm, file="hermes_8d_log_aggregate_extracted_chl.RData")

### log normalize data, take mean, exponetiate back to values for OCCCI v5
w.occi.nes=extract_calc_logmean(occi.v5[[1:dim(occi.v5)[3]]], nes.five)
kf.occi.v5.logmean=occi.date #%>% select(month, day, year, F1, DOY, diff)
kf.occi.v5.logmean$nes=w.occi.nes[1,]
kf.occi.v5.logmean$gbk=w.occi.nes[2,]
kf.occi.v5.logmean$gome=w.occi.nes[3,]
kf.occi.v5.logmean$gomw=w.occi.nes[4,]
kf.occi.v5.logmean$mabn=w.occi.nes[5,]
kf.occi.v5.logmean$mabs=w.occi.nes[6,]
save(kf.occi.v5.logmean, file="/home/ryan/Git/JPSS/occci_v5_8d_5area_logmean_extracted_chl.RData")

w.occi.nes=extract_calc_geo(occi[[1:dim(occi)[3]]], nes.five)
kf.occi=occi.date %>% select(month, day, year, F1, DOY, diff)
# kf.occi=occi.date
w.occi.nes=extract_calc(occi[[1:1028]], nes.five)
kf.occi=occi.date %>% select(month, day, year, F1, DOY, diff)
kf.occi$nes=w.occi.nes[1,]
kf.occi$gbk=w.occi.nes[2,]
kf.occi$gome=w.occi.nes[3,]
kf.occi$gomw=w.occi.nes[4,]
kf.occi$mabn=w.occi.nes[5,]
kf.occi$mabs=w.occi.nes[6,]
# kf.occi$year=year(kf.occi$test)
# kf.occi$month=month(kf.occi$test)
# kf.occi$day=day(kf.occi$test)
# kf.occi$DOY=yday(kf.occi$test)

save(kf.occi, file="/home/ryan/Git/JPSS/occi_v5__8d_5area_extracted_chl.RData")
save(kf.occi.geo, file="/home/ryan/Git/JPSS/occi_v5_8d_5area_geomean_extracted_chl.RData")

## plot OCCI for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.occi, month, year, gomw) %>% group_by(year) %>% summarise(mean=gm_mean(gomw))
plot(tt, type='b', main=paste("v4 GOMw OCCI annual Chl"))
tt.sp=dplyr::select(kf.occi, month, year, gomw) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(gomw))
plot(tt.sp, type='b',main=paste("v4 GOMw OCCI Jan-Jun Chl"))
tt.fl=dplyr::select(kf.occi, month, year, gomw) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(gomw))
plot(tt.fl, type='b', main=paste("v4 GOMw OCCI Jul-Dec Chl"))

## plot v5 OCCI for NES, annual, first 6, last 6 months
tt=dplyr::select(new, month, year, gomw) %>% group_by(year) %>% summarise(mean=mean(gomw))
plot(tt, type='b', main=paste("v5 GOMw OCCI annual Chl"))
tt.sp=dplyr::select(new, month, year, gomw) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(gomw))
plot(tt.sp, type='b',main=paste("v5 GOMw OCCI Jan-Jun Chl"))
tt.fl=dplyr::select(new, month, year, gomw) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(gomw))
plot(tt.fl, type='b', main=paste("v5 GOMw OCCI Jul-Dec Chl"))

## plot Hermes GSM for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.gsm, M, Y, nes) %>% group_by(Y) %>% summarise(mean=gmean(nes))
save(kf.occi, file="occi_8d_extracted_chl.RData")

## plot OCCI for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.occi, month, year, nes) %>% group_by(year) %>% summarise(mean=mean(nes))
plot(tt, type='b', main=paste("NES OCCI annual Chl"))
tt.sp=dplyr::select(kf.occi, month, year, nes) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(nes))
plot(tt.sp, type='b',main=paste("NES OCCI Jan-Jun Chl"))
tt.fl=dplyr::select(kf.occi, month, year, nes) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(nes))
plot(tt.fl, type='b', main=paste("NES OCCI Jul-Dec Chl"))

## plot Hermes GSM for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.gsm, M, Y, nes) %>% group_by(Y) %>% summarise(mean=mean(nes))
plot(tt, type='b', main=paste("GSM NES annual Chl"))
tt.sp=dplyr::select(kf.gsm, M, Y, nes) %>% filter(M<=6) %>% group_by(Y) %>% summarise(mean=mean(nes))
plot(tt.sp, type='b',main=paste("GSM NES Jan-Jun Chl"))
tt.fl=dplyr::select(kf.gsm, M, Y, nes) %>% filter(M>6) %>% group_by(Y) %>% summarise(mean=mean(nes))
plot(tt.fl, type='b', main=paste("GSM NES Jul-Dec Chl"))

## plot log-transform aggregrated and then exponentiated Hermes data for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.lg.gsm, M, Y, nes) %>% group_by(Y) %>% summarise(mean=mean(nes))
plot(tt, type='b', main=paste("GSM log NES annual Chl"))
tt.sp=dplyr::select(kf.lg.gsm, M, Y, nes) %>% filter(M<=6) %>% group_by(Y) %>% summarise(mean=mean(nes))
plot(tt.sp, type='b',main=paste("GSM log NES Jan-Jun Chl"))
tt.fl=dplyr::select(kf.lg.gsm, M, Y, nes) %>% filter(M>6) %>% group_by(Y) %>% summarise(mean=mean(nes))
plot(tt.fl, type='b', main=paste("GSM log NES Jul-Dec Chl"))

### plot GBK mean extracted values, averaged with geomean by month
tt=dplyr::select(kf.occi.v5, month, year, gbk) %>% group_by(month, year) %>% summarise(mean=gm_mean(gbk))
tt$date=as.Date(paste0("15/", paste(tt$month, tt$year, sep='/')),format = "%d/%m/%Y")
tt=tt[order(tt$date),]
plot(tt$mean ~ tt$date, type='l', col='red', ylab='GBK 8-d mean Chl', las=1)
## plot OCCI for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.occi.v4, month, year, gbk) %>% group_by(month, year) %>% summarise(mean=gm_mean(gbk))
tt$date=as.Date(paste0("15/", paste(tt$month, tt$year, sep='/')),format = "%d/%m/%Y")
tt=tt[order(tt$date),]
lines(tt$date, tt$mean, col='blue')
legend('topleft', legend=c("v 5.0", "v 4.2"),
       col=c("red", "blue"), lty=c(1:2), box.lty=0)


### Open SAHFOS CPR data
cprgome=readxl::read_excel('/home/ryan/Desktop/RyanMorse_CPRExtract/CPR_Extract_Data_nes.xlsx'); loc='NES'
cprgometx=readxl::read_excel('/home/ryan/Desktop/RyanMorse_CPRExtract/CPR_Extract_TaxaList_nes.xlsx')

cprgome=readxl::read_excel('/home/ryan/Desktop/RyanMorse_CPRExtract/CPR_Extract_Data_gome.xlsx'); loc='GOMe'
cprgometx=readxl::read_excel('/home/ryan/Desktop/RyanMorse_CPRExtract/CPR_Extract_TaxaList_gome.xlsx')
library(maps)
library(mapdata)
library(marmap)
library(mapproj)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
data(stateMapEnv)
map('state', fill = F, add=T) # add state lines
points(cprgome$Longitude[which(cprgome$Year==1991)], cprgome$Latitude[which(cprgome$Year==1991)], pch=16)
points(cprgome$Longitude, cprgome$Latitude, pch=16)


#cprx=dplyr::select(cprgome, Month, Year) %>% group_by(Year) %>% summarise(ct=sum(Month))
#cprx=dplyr::select(cprgome, Month, Year, Latitude, Longitude, PCI) %>% group_by(Year, Latitude, Longitude) %>% summarise(ct=sum(Month), mn=mean(PCI))
cprx=dplyr::select(cprgome, Month, Year, Latitude, Longitude, PCI) %>% group_by(Year, Month) %>% summarise(ct=sum(Month), mn=mean(PCI))
cprx=dplyr::select(cprgome, Month, Year, Latitude, Longitude, PCI) %>% group_by(Year) %>% summarise(ct=sum(Month), mn=mean(PCI))
# table(cprx$Year)
cprx$D=rep(15, length(cprx$Year))
cprx=cprx %>% mutate(date = make_date(Year, Month, D))

cprgome=readxl::read_excel('/home/ryan/Desktop/RyanMorse_CPRExtract/CPR_Extract_Data_mabn.xlsx'); loc='MAB'
cprgometx=readxl::read_excel('/home/ryan/Desktop/RyanMorse_CPRExtract/CPR_Extract_TaxaList_mabn.xlsx')

##yearly
cprx=dplyr::select(cprgome, Month, Year, Latitude, Longitude, PCI) %>% group_by(Year) %>% summarise(ct=sum(Month), mn=mean(PCI))
par(mfrow=c(1,2))
plot(cprx$mn ~cprx$Year, type='b', main='Annual mean')
barplot(cprx$ct, names.arg=cprx$Year, main='sample count')
## first 6 months
cprx=dplyr::select(cprgome, Month, Year, Latitude, Longitude, PCI) %>% subset(Month <6) %>% group_by(Year) %>% summarise(ct=sum(Month), mn=mean(PCI))
par(mfrow=c(1,2))
plot(cprx$mn ~cprx$Year, type='b', main='Jan-May mean')
barplot(cprx$ct, names.arg=cprx$Year, main='sample count')
## last 6 months
cprx=dplyr::select(cprgome, Month, Year, Latitude, Longitude, PCI) %>% subset(Month >5) %>% group_by(Year) %>% summarise(ct=sum(Month), mn=mean(PCI))
par(mfrow=c(1,2))
plot(cprx$mn ~cprx$Year, type='b', main='Jun-Dec mean')
barplot(cprx$ct, names.arg=cprx$Year, main='sample count')

test=cprgometx%>% subset(Is_Diatom==1)
cprx=cprgome%>% dplyr::select( Year, one_of(as.character(test$Taxon_Id)))
cprx$diatoms=rowSums(cprx[2:dim(cprx)[2]], na.rm=T)
t=cprx %>% group_by(Year) %>% summarise(mn=mean(diatoms, na.rm=T), sm=sum(diatoms, na.rm=T), sd=sd(diatoms, na.rm=T), ct=n())
plot(t$mn~t$Year, type='b', ylab='diatom mean count', xlab='', main=loc) 
sdup=t$mn+t$sd
sddn=t$mn-t$sd
lines(t$Year, sdup, col='gray90')
barplot(t$ct, names.arg=t$Year, main='sample count', las=2)



### Read in CZCS data
# https://oceandata.sci.gsfc.nasa.gov/directaccess/CZCS/Mapped/Monthly_Climatology/9km/chlor_a/
nc1=nc_open('/home/ryan/Desktop/Z/CZCS/monthly_clim/C19782741985304.L3m_MC_CHL_chlor_a_9km.nc')
nc1=nc_open('/home/ryan/Downloads/C19783051978334.L3m_MO_CHL_chlor_a_9km.nc')
nc1=nc_open('/home/ryan/Desktop/Z/CZCS/monthly/C19792131979243.L3m_MO_CHL_chlor_a_9km.nc')
lon.czcs=ncvar_get(nc1, 'lon')
lat.czcs=ncvar_get(nc1, 'lat')
chl.czcs=ncvar_get(nc1, 'chlor_a')
lg.chl.czcs=log10(chl.czcs)
time.czcs=ncvar_get(nc1, 'time') #days since Jan 1, 1970
dim(chl.czcs)
colnames(chl.czcs)=as.numeric(lat.czcs)
rownames(chl.czcs)=as.numeric(lon.czcs)
nc_close(nc1)
m2=t(chl.czcs)
xx=raster(m2)
plot(xx)

# subset
bb=c(-80, -60, 32, 48)
LonIdx <- which( nc1$dim$lon$vals > -80 & nc1$dim$lon$vals < -60)
LatIdx <- which( nc1$dim$lat$vals > 32 & nc1$dim$lat$vals < 48)
chl.czcs <- ncvar_get( nc1, 'chlor_a')[ LonIdx, LatIdx]
lon.czcs=ncvar_get(nc1, 'lon')[LonIdx]
lat.czcs=ncvar_get(nc1, 'lat')[LatIdx]
colnames(chl.czcs)=as.numeric(lat.czcs)
rownames(chl.czcs)=as.numeric(lon.czcs)
m2=t(chl.czcs)
xx=raster(m2)
extent(xx)=bb
crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
plot(xx, zlim=c(0,20))

# Monthly Climatology - take geometric mean of values in box
wd='/home/ryan/Desktop/Z/CZCS/monthly_clim/'
files.chl1=list.files(wd, pattern=('_9km.nc'))
bb=c(-80, -60, 32, 48)
LonIdx <- which( nc1$dim$lon$vals > -80 & nc1$dim$lon$vals < -60)
LatIdx <- which( nc1$dim$lat$vals > 32 & nc1$dim$lat$vals < 48)
ttxm=matrix(NA,nrow = 30, ncol=12)
for(i in c(seq(4,12,1), seq(1,3,1))){
  nc1=nc_open(paste(wd,files.chl1[i], sep=''))
  chl.czcs <- ncvar_get( nc1, 'chlor_a')[ LonIdx, LatIdx]
  lon.czcs=ncvar_get(nc1, 'lon')[LonIdx]
  lat.czcs=ncvar_get(nc1, 'lat')[LatIdx]
  colnames(chl.czcs)=as.numeric(lat.czcs)
  rownames(chl.czcs)=as.numeric(lon.czcs)
  m2=t(chl.czcs)
  xx=raster(m2)
  extent(xx)=bb
  crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  # tt=extract_calc_geo(xx, neus.shp)
  tt=extract_calc(xx, neus.shp, fun = 'median', na.option = T)
  ttxm[,i]=tt
}
## NOW CORRECT FOR ORDER ISSUE IN NEUS.SHP
test=neus.shp$BOX_ID
f=as.numeric(levels(test))[test] #box numbers of ttx
ttxm=ttxm[order(f+1),] # winner
# save(ttxm, file='/home/ryan/Desktop/Z/CZCS/czcs_monthly_extracted.rdata')
save(ttxm, file='/home/ryan/Desktop/Z/CZCS/czcs_monthly_extracted_median.rdata')

# Monthly data - take geometric mean of values in box
wd='/home/ryan/Desktop/Z/CZCS/monthly/'
files.chl1=list.files(wd, pattern=('_9km.nc'))
bb=c(-80, -60, 32, 48)
LonIdx <- which( nc1$dim$lon$vals > -80 & nc1$dim$lon$vals < -60)
LatIdx <- which( nc1$dim$lat$vals > 32 & nc1$dim$lat$vals < 48)
ttx=matrix(NA,nrow = 30, ncol=length(files.chl1))
for(i in 1:length(files.chl1)){
  nc1=nc_open(paste(wd,files.chl1[i], sep=''))
  chl.czcs <- ncvar_get( nc1, 'chlor_a')[ LonIdx, LatIdx]
  lon.czcs=ncvar_get(nc1, 'lon')[LonIdx]
  lat.czcs=ncvar_get(nc1, 'lat')[LatIdx]
  colnames(chl.czcs)=as.numeric(lat.czcs)
  rownames(chl.czcs)=as.numeric(lon.czcs)
  m2=t(chl.czcs)
  xx=raster(m2)
  extent(xx)=bb
  crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  tt=extract_calc(xx, neus.shp, fun = 'median', na.option = T)
  ttx[,i]=tt
}
## NOW CORRECT FOR ORDER ISSUE IN NEUS.SHP
test=neus.shp$BOX_ID
f=as.numeric(levels(test))[test] #box numbers of ttx
ttx2=ttx[order(f+1),] # winner
save(ttx2, file='/home/ryan/Desktop/Z/CZCS/czcs_monthly_extracted_median.rdata')

table(czcs.dt$X1)

s02 <- smooth.spline(ttx[1,], spar = 0.2)
xx  <- seq(1,length(ttx[1,]), len=365)
sp=predict(s02, xx)
sp$x=seq(1,365, 1)

# check on method comparability -> geomean gives 1.00000 when others yield NA
i=49 # hardly any data
i=46 #complete coverage
tes=extract_calc(xx, neus.shp, fun='median', na.option = F)
tes1=extract_calc(xx, neus.shp, fun='median', na.option = T)
tes2=extract_calc(xx, neus.shp, fun='mean', na.option = T)
tes3=extract_calc_geo(xx, neus.shp)
tes4=extract_calc_logmean(xx, neus.shp)
tes49=matrix(rbind(tes, tes1, tes2, tes3, tes4), ncol=5, nrow=30)
t49=tes49[order(f+1),] # winner
tes46=matrix(rbind(tes, tes1, tes2, tes3, tes4), ncol=5, nrow=30)
t46=tes46[order(f+1),]
colnames(t46)=c('medNA','med','mean','geo','logmn')
colnames(t49)=c('medNA','med','mean','geo','logmn')
plot(t46[,2], type='l')
lines(t46[,3], col='red')
lines(t46[,4], col='green')
lines(t46[,5], col='yellow')
legend('topleft', lty=c(1,1,1,1),col=c('black', 'red', 'green', 'yellow'), 
       legend=c('med', 'mean', 'geomn', 'logmn'), bty='n')
plot(t49[,2], type='l')
lines(t49[,3], col='red')
lines(t49[,4], col='green')
lines(t49[,5], col='yellow')
legend('topleft', lty=c(1,1,1,1),col=c('black', 'red', 'green', 'yellow'), 
       legend=c('med', 'mean', 'geomn', 'logmn'), bty='n')

## create spline for days in between points
##-- artificial example
# y18 <- c(1:3,5,4,7:3,2*(2:5),rep(10,4))
# xx  <- seq(1,length(y18), len=201)
# (s2  <- smooth.spline(y18)) # GCV
# (s02 <- smooth.spline(y18, spar = 0.2))
# plot(y18, main=deparse(s2$call), col.main=2)
# lines(s2, col = "gray"); lines(predict(s2, xx), col = 2)
# lines(predict(s02, xx), col = 3); mtext(deparse(s02$call), col = 3)

#drop bad data 1978-Oct, add in 2 missing months
ttx2=ttx2[2:91,]
newdt=newdt[2:91,]
#now rbind date to add missing 8/1981, 5/1982 to newdt
newdt2=rbind(newdt[1:33,],newdt[33,],newdt[34:42,],newdt[42,],newdt[43:90,])
newdt$ord=c(seq(1,33,1), seq(35,43,1), seq(45,92,1))
newdt2=rbind(newdt, c(1981, 214,244,229,230, 227,8,34))
newdt2=rbind(newdt2, c(1982,182,212,197,200,196,7,44))
# colnames(ttx2)=newdt$ord
bad=data.frame(matrix(NA,ncol=2, nrow=30))
satmon=data.frame(ttx2, bad)
satmon2=satmon[,order(newdt2$ord)] # add in missing months with NAs in proper order
czcs.dt.final=newdt2[order(newdt2$ord),]

## fill NA with climatology by row (i) and column(j):
satmon3=satmon2
for (i in 1:dim(satmon2)[1]){
  for(j in 1:dim(satmon2)[2]){
    if(bad[i,j]==1){
      satmon3[i,j]=ttxm[i,czcs.dt.final$M[j]]
    }
  }
}

# s1978=satmon3[,1:2]
s1979=data.frame(satmon3[,3:14])
s1980=satmon3[,15:26]
s1981=satmon3[,27:38]
s1982=satmon3[,39:50]
s1983=satmon3[,51:62]
s1984=satmon3[,63:74]
s1985=satmon3[,75:86]
# s1986=satmon3[,87:92]

interpnoise=function(d, name=deparse(substitute(d))){
ttxd=matrix(NA, ncol=365, nrow=30) # daily interpolated
ttxdn=matrix(NA, ncol=365, nrow=30) # daily interp plus white noise
for(i in 1:30){
  s02 <- smooth.spline(as.numeric(d[i,]), spar = 0.2)
  xx  <- seq(1,length(d[i,]), len=365)
  sp=predict(s02, xx)
  spb <- noise(sp$y,'white',0.5)
  ttxd[i,]=sp$y
  ttxdn[i,]=spb
}
# replace where values are below 0
min(ttxdn)
ttxdn[ttxdn<0]=0.001
save(ttxdn, file=paste('/home/ryan/Desktop/Z/CZCS/',name,'noise.Rdata'))
save(ttxd, file=paste('/home/ryan/Desktop/Z/CZCS/',name,'.Rdata'))
}

interpnoise(s1979)
interpnoise(s1980)
interpnoise(s1981)
interpnoise(s1982)
interpnoise(s1983)
interpnoise(s1984)
interpnoise(s1985)

czcsclim=ttxm
interpnoise(czcsclim)

#visualize chl per box as year long time series
plot(ttxdn[2,], type='l', ylim=c(0,20)) #max(ttxd)))
for(i in 3:23){
  lines(ttxdn[i,])
}

DOY=c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # Day 15 of month generic year
plot(ttx[1,], type='p')
s02 <- smooth.spline(ttx[1,], spar = 0.2)
xx  <- seq(1,length(ttx[1,]), len=365)
sp=predict(s02, xx)
sp$x=seq(1,365, 1)

plot(ttx[1,]~DOY, type='p')
lines(sp, col='red')

# func = splinefun(ttx[1,], y=seq(1,12,1), method="natural",  ties = mean)  
# x.daily=func(seq(1, 12, 0.0302)) #365 days

### Add noise to signal
#Generate a HS patterns
library(rpatrec)
# a <- generator()
# #now add white noise with a standdard deviation of 10
# b <- noise(a,'white',10)
# plot(a, type='b')
# lines(b, col='red')
spb <- noise(sp$y,'white',0.5)
plot(ttx[1,]~DOY, type='p', ylim=c(0,7))
lines(spb, col='green')

# now fill monthly to daily with interpolated + noise
ttxd=matrix(NA, ncol=365, nrow=30) # daily interpolated
ttxdn=matrix(NA, ncol=365, nrow=30) # daily interp plus white noise
for(i in 1:30){
s02 <- smooth.spline(ttx[i,], spar = 0.2)
xx  <- seq(1,length(ttx[i,]), len=365)
sp=predict(s02, xx)
spb <- noise(sp$y,'white',0.5)
ttxd[i,]=sp$y
ttxdn[i,]=spb
}
# replace where values are below 0
min(ttxdn)
ttxdn[ttxdn<0]=0.001
save(ttxdn, file='/home/ryan/Desktop/Z/CZCS/DailyInterpWithNoiseCZCSclimatology.Rdata')
save(ttxd, file='/home/ryan/Desktop/Z/CZCS/DailyInterpCZCSclimatology.Rdata')

CZCSclimnoise=loadRData("/home/ryan/Desktop/Z/CZCS/DailyInterpWithNoiseCZCSclimatology.Rdata")



# take geometric mean of values in box
wd='/home/ryan/Desktop/Z/CZCS/monthly/'
files.chl1=list.files(wd, pattern=('_9km.nc'))

# C1978 274 1978 304.
# czcs.dt.mon=list()
czcs.dt.mon=matrix(NA, nrow=length(files.chl1), ncol=1)
for (i in 1:length(files.chl1)){
  czcs.dt.mon[i]=strsplit(files.chl1[i],split=".", fixed=TRUE)[[1]][1]
}
czcs.dt=matrix(NA, nrow=length(files.chl1), ncol=3)
czcs.dt=data.frame(czcs.dt)
czcs.dt$X1=as.numeric(substr(czcs.dt.mon[,1], 2,5)) #YYYY
czcs.dt$X2=as.numeric(substr(czcs.dt.mon[,1], 6,8)) #DOYstart
czcs.dt$X3=as.numeric(substr(czcs.dt.mon[,1], 13,15)) #DOYend
czcs.dt$med=floor(apply(czcs.dt[,2:3], 1, FUN='median'))
doy=data.frame(DOY)
doy$M=c(seq(1,12,1))
doy$R=round(doy$DOY,-1)

newdt=left_join(czcs.dt, doy, by='R')
sum(is.na(newdt$M))
table(newdt$R[(is.na(newdt$M))])
newdt$M[newdt$R==40]=2
newdt$M[newdt$R==80]=3
newdt$M[newdt$R==110]=4


bb=c(-80, -60, 32, 48)
LonIdx <- which( nc1$dim$lon$vals > -80 & nc1$dim$lon$vals < -60)
LatIdx <- which( nc1$dim$lat$vals > 32 & nc1$dim$lat$vals < 48)
ttx=matrix(NA,ncol = 30, nrow=length(czcs.dt.mon))
for (i in 1:length(files.chl1)){
  nc1=nc_open(paste(wd,files.chl1[i], sep=''))
  chl.czcs <- ncvar_get( nc1, 'chlor_a')[ LonIdx, LatIdx]
  lon.czcs=ncvar_get(nc1, 'lon')[LonIdx]
  lat.czcs=ncvar_get(nc1, 'lat')[LatIdx]
  colnames(chl.czcs)=as.numeric(lat.czcs)
  rownames(chl.czcs)=as.numeric(lon.czcs)
  m2=t(chl.czcs)
  xx=raster(m2)
  extent(xx)=bb
  crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  tt=extract_calc_geo(xx, neus.shp)
  ttx[i,]=tt
}
## neus.shp is not in order, affects order of final matrix, need to sort
test=neus.shp$BOX_ID
f=as.numeric(levels(test))[test] #box numbers of ttx
f2=ttx[,order(f+1)] # winner
# f2=ttx[,(f+1)]
# f2=ttx[,f]
# f2=ttx[,sort(f+1)]

#testing extraction of certain boxes
plot(xx)
polygon(neus.shp@polygons[[1]]@Polygons[[1]]@coords) # draw box 2
extract_calc_geo(xx, neus.shp@polygons[[1]]@Polygons[[1]]@coords)

#  Plot images to check for missind values
# pdf(file='/home/ryan/Desktop/Z/CZCS/CZCSimages.pdf')
# pdf(file='/home/ryan/Desktop/Z/CZCS/CZCSclimatologyImages.pdf')
# for (i in 1:length(files.chl1)){
#   nc1=nc_open(paste(wd,files.chl1[i], sep=''))
#   chl.czcs <- ncvar_get( nc1, 'chlor_a')[ LonIdx, LatIdx]
#   lon.czcs=ncvar_get(nc1, 'lon')[LonIdx]
#   lat.czcs=ncvar_get(nc1, 'lat')[LatIdx]
#   colnames(chl.czcs)=as.numeric(lat.czcs)
#   rownames(chl.czcs)=as.numeric(lon.czcs)
#   m2=t(chl.czcs)
#   xx=raster(m2)
#   extent(xx)=bb
#   crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#   plot(xx, zlim=c(0,10), main=paste('Y:',newdt$X1[i], ' M:',newdt$M[i], ' i:',i, sep=' '))
#   lines(neus.shp)
# }
# dev.off()


czcs.mon=heatmap.2(ttxm[2:23,], Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                   labCol = seq(1,12,1), labRow =seq(1,22,1), col=my_palette,
                   tracecol = 'black', main='CZCS Monthly')
czcs.mon=heatmap.2(ttxm, Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                   labCol = seq(1,12,1), labRow =seq(0,29,1), col=my_palette,
                   tracecol = 'black', main='CZCS Monthly')

czcs.d=heatmap.2(t(f2[,2:23]), Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                   labRow = seq(1,22,1), labCol =czcs.dt$X1, col=my_palette,
                   tracecol = 'black', main='CZCS')

### heatmaps
library(gplots)
some_col_func <- colorspace::diverge_hcl
my_palette <- colorRampPalette(c("white", "green"))(n = 7)
my_palette2 <- colorRampPalette(c("blue", "yellow", "green"))(n = 20)
st.anom=apply(dataT, 2, function(x) (x-mean(x))/std(x)) # calc standardized anomaly of data

## monthly Chl value for MARMAP
my_palette <- colorRampPalette(c("white", "red"))(n = 7)
mm.mon.anom=heatmap.2(as.matrix(MM.boxes.mon[,2:13]), breaks=c(0,1,2,3,4,5,6,7),Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                      labCol = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), labRow = MM.boxes.mon[,1], col=my_palette,
                      tracecol = 'black', main='Mean Monthly Chl 1977-1987')
czcs.clim=heatmap.2(ttx, Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                      labCol = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), labRow = seq(0,29,1), col=my_palette2,
                      tracecol = 'black', main='CZCS Climatology')
czcs.clim=heatmap.2(ttx[2:23,], Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    labCol = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), labRow = seq(1,22,1), col=my_palette2,
                    tracecol = 'black', main='CZCS Climatology')

czcs.mon=heatmap.2(f2, Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    labCol = seq(0,29,1), labRow =czcs.dt$X1, col=my_palette2,
                    tracecol = 'black', main='CZCS Monthly')
czcs.mon=heatmap.2(t(f2[,2:23]), Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                   labRow = seq(1,22,1), labCol =czcs.dt$X1, col=my_palette,
                   tracecol = 'black', main='CZCS Monthly')

czcs.clim=heatmap.2(ttx[,2:23], Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    col=my_palette2, labCol = seq(1,22,1),tracecol = 'black')

## testing MARMAP coverage by month
# test=MM.boxes.nano.all
test=MM.boxes.ts
test$D=15
dt=paste(test$Y,test$M, test$D, sep='-')
test$dt=as.Date(dt)
test2=test[,c('box', 'x', 'dt')]
test3=test2 %>% tidyr::spread(dt, x, fill = NA, convert = FALSE)
t=gplots::heatmap.2(as.matrix(test3[,2:96]), Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    col=my_palette2, tracecol = 'black')

# create matrix with monthly values from 1977-1987 to merge values with
t1=matrix(NA, nrow=30*12*11, ncol=4)
## boxes
t1a=rep(1,12)
for(n in 1:29){
  i=n+1
  t1b=rep(i,12)
  t1a=c(t1a,t1b)
}
# t1[,1]=rep(seq(1,30,1),132) #box 1-30
t1[,1]=rep(t1a,11) # boxes 1-30
t1[,2]=rep(seq(1,12,1),330) #M
t1y=c(rep(1977,360),rep(1978,360),rep(1979,360),rep(1980,360),rep(1981,360),rep(1982,360),rep(1983,360),rep(1984,360),rep(1985,360),rep(1986,360),rep(1987,360))#Y
t1[,3]=t1y # Years
colnames(t1)=c('box', 'M', 'Y','x')
t1=data.frame(t1)
MM.b.clim=dplyr::left_join(t1,MM.boxes.ts, by = c("box", "M", "Y"))

test=MM.b.clim
test$D=15
dt=paste(test$Y,test$M, test$D, sep='-')
test$dt=as.Date(dt)
test2=test[,c('box', 'x.y', 'dt')]
MM.box.mon.ts=test2 %>% tidyr::spread(dt, x.y, fill = NA, convert = FALSE)
t=gplots::heatmap.2(as.matrix(test3[,2:96]), Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    col=my_palette2, tracecol = 'black')



### create regimes:
## plot 1 year from 1 box, repeat for long term base
test=ttxdn[13,]
t2=rep(test, 20)
plot(t2, type='l')


