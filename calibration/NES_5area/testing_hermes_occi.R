library(lubridate)
library(raster)
library(dplyr)

### Function to extract data using a shapefile
extract_calc=function(x, shp){
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

test=month.day.year(time.occi, c(1,1,1970)) # these are not 8-days apart.... something odd -> see manual
occi.date=data.frame(test)
occi.date$F1=paste(occi.date[,3], occi.date[,1], occi.date[,2], sep='-')
occi.date$DOY=as.numeric(strftime(occi.date$F1, '%j'))
ddiff=diff(occi.date$DOY)
occi.date$diff=c(0, ddiff) # see OC-CCI manual in JPSS/calibration folder for list of missing dates, explains why some are not 8d

table(occi.date$year)
occi.date$week=c(seq(from=31, to=46, by=1),  rep(seq(from=1, to=46, by=1), 21))


## Chl_loop over and stack rasters
bb=c(-80, -60, 32, 48)
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

w.occi.nes=extract_calc(occi[[1:1028]], nes.five)
kf.occi=occi.date %>% select(month, day, year, F1, DOY, diff)
kf.occi$nes=w.occi.nes[1,]
kf.occi$gbk=w.occi.nes[2,]
kf.occi$gome=w.occi.nes[3,]
kf.occi$gomw=w.occi.nes[4,]
kf.occi$mabn=w.occi.nes[5,]
kf.occi$mabs=w.occi.nes[6,]
save(kf.occi, file="occi_8d_extracted_chl.RData")

## plot OCCI for NES, annual, first 6, last 6 months
# tt=dplyr::select(kf.occi, month, year, nes) %>% group_by(year) %>% summarise(mean=mean(nes))
tt=dplyr::select(kf.occi, month, year, nes) %>% group_by(year,month) %>% summarise(mean=mean(nes)) %>%
  group_by(year)%>% summarise(mean=mean(mean))
plot(tt, type='b', main=paste("NES annual Chl"))
tt.sp=dplyr::select(kf.occi, month, year, nes) %>% filter(month<=6) %>% group_by(year, month) %>% summarise(mean=mean(nes))%>%
  group_by(year)%>% summarise(mean=mean(mean))
plot(tt.sp, type='b',main=paste("NES Jan-Jun Chl"))
tt.fl=dplyr::select(kf.occi, month, year, nes) %>% filter(month>6) %>% group_by(year, month) %>% summarise(mean=mean(nes))%>%
  group_by(year)%>% summarise(mean=mean(mean))
plot(tt.fl, type='b', main=paste("NES Jul-Dec Chl"))

## plot Hermes GSM for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.gsm, M, Y, nes) %>% group_by(Y, M) %>% summarise(mean=mean(nes))%>%
  group_by(Y)%>% summarise(mean=mean(mean))
plot(tt, type='b', main=paste("GSM NES annual Chl"))
tt.sp=dplyr::select(kf.gsm, M, Y, nes) %>% filter(M<=6) %>% group_by(Y, M) %>% summarise(mean=mean(nes))%>%
  group_by(Y)%>% summarise(mean=mean(mean))
plot(tt.sp, type='b',main=paste("GSM NES Jan-Jun Chl"))
tt.fl=dplyr::select(kf.gsm, M, Y, nes) %>% filter(M>6) %>% group_by(Y, M) %>% summarise(mean=mean(nes))%>%
  group_by(Y)%>% summarise(mean=mean(mean))
plot(tt.fl, type='b', main=paste("GSM NES Jul-Dec Chl"))

## plot log-transform aggregrated and then exponentiated Hermes data for NES, annual, first 6, last 6 months
tt=dplyr::select(kf.lg.gsm, M, Y, nes) %>% group_by(Y, M) %>% summarise(mean=mean(nes))%>%
  group_by(Y)%>% summarise(mean=mean(mean))
plot(tt, type='b', main=paste("GSM log NES annual Chl"))
tt.sp=dplyr::select(kf.lg.gsm, M, Y, nes) %>% filter(M<=6) %>% group_by(Y, M) %>% summarise(mean=mean(nes))%>%
  group_by(Y)%>% summarise(mean=mean(mean))
plot(tt.sp, type='b',main=paste("GSM log NES Jan-Jun Chl"))
tt.fl=dplyr::select(kf.lg.gsm, M, Y, nes) %>% filter(M>6) %>% group_by(Y, M) %>% summarise(mean=mean(nes))%>%
  group_by(Y)%>% summarise(mean=mean(mean))
plot(tt.fl, type='b', main=paste("GSM log NES Jul-Dec Chl"))
