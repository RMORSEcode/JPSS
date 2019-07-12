 # NOTE 20190710
 # These are not yet working -- use ftp_download.py in dropbox folder to download HERMES GlobColour data
#  require(RCurl)
# roms_file='C:/Users/ryan.morse/Desktop/NEUS Atl files/HYDRO/roms2010all.nc'
# 
# file_db <- bind_rows(lapply(roms_file, function(x) {
#   nc <- NetCDF(x)
#   tlen <- filter(nc$dimension, name == "ocean_time")$len
#   tibble(fullname = rep(x, tlen), band_level = seq_len(tlen))
# }))


# url = "ftp://ftp_hermes:hermes%@ftp.hermes.acri.fr/345576024" #ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1297/suppl/"
# filenames = getURL(url, ftp.use.epsv = T, dirlistonly = TRUE)
# filenames <- strsplit(filenames, "\r\n")
# filenames = unlist(filenames)
# filenames
# [1] "filelist.txt"    "GSE1297_RAW.tar"
# for (filename in filenames) {
#   download.file(paste(url, filename, sep = ""), paste(getwd(), "/", filename,
#                                                       sep = ""))
# }
# 
# url<- "ftp://ftp.hermes.acri.fr/345576024"
# filenames <- getURL(url, userpwd="ftp_hermes:hermes%", ftp.use.epsv = FALSE, dirlistonly = TRUE) #reading filenames from ftp-server
# destnames <- filenames <-  strsplit(filenames, "\r*\n")[[1]] # destfiles = origin file names
# con <-  getCurlHandle( ftp.use.epsv = FALSE, userpwd="ftp_hermes:hermes%")
# setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES_daily')
# mapply(function(x,y) writeBin(getBinaryURL(x, curl = con, dirlistonly = FALSE), y), x = filenames, y = paste("C:\\temp\\",destnames, sep = "")) #writing all zipped files in one directory
# 
# 
# devtools::install_github("skgrange/threadr")
# 
# library(threadr)
# download_ftp_file(url, destnames, credentials = "ftp_hermes:hermes%",
#                   curl = T, verbose = FALSE, progress = "none")

# OCCI https://esa-oceancolour-cci.org/
# NCSS Request URL
# https://rsg.pml.ac.uk/thredds/ncss/grid/CCI_ALL-v4.0-8DAY/dataset.html
# /thredds/ncss/CCI_ALL-v4.0-8DAY
# https://rsg.pml.ac.uk/thredds/ncss/CCI_ALL-v4.0-8DAY?var=chlor_a&north=48&west=-80&east=-60&south=32&horizStride=1&time_start=1997-09-04T00%3A00%3A00Z&time_end=2018-12-27T00%3A00%3A00Z&timeStride=1&addLatLon=true

library(raster)
library(mgcv)
library(sp)
library(maptools)
library(marmap)
library(rgeos)
library(ncdf4)
library(abind)

### Function to load ncdf files as stacked raster
nc2raster=function(x, varn){
  s=stack()
  for (i in 1:length(x)){
    r <- raster(x[i],  varname = varn)
    s=stack(s, r)
  }
  return(s)
}


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








# ## load shapefiles
# setwd("G:/1 RM/KINGSTON/transfer/shapefiles/epu_shapes")
setwd("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/KINGSTON/transfer/shapefiles/epu_shapes")
wd=getwd()
gbk=rgdal::readOGR("EPU_GBKPoly.shp")
gom=rgdal::readOGR("EPU_GOMPoly.shp")
mab=rgdal::readOGR("EPU_MABPoly.shp")
scs=rgdal::readOGR("EPU_SCSPoly.shp")
#combine shapefiles GOM and GBK
gom.gbk.shp=gUnion(gom, gbk, byid=F, id=NULL)
gom.gbk.shp=gUnion(gom, gbk, byid=F, id=NULL)
gom.scs.shp=gUnion(gom, scs, byid=F, id=NULL)
mab.gbk.shp=gUnion(mab, gbk, byid=F, id=NULL)
NES.shp=gUnion(mab.gbk.shp, gom.scs.shp, byid=F, id=NULL)

#### Load Chl data
# AVW: weighted average of single-sensor Level 2 CHL1 products
# GSM: GSM merging of single sensor L3 NRRS
# The CHL1 algorithms are applicable for "case 1" waters
# OC5 Gohin, F., 2011
# CHL2 is the chlorophyll concentration (mg/m3) for Case 2 waters (see section validity); L3 merge: AV; sensors: MER, OLA; Doerffer and Schiller (2007)
# CHL2 uses the a Neural Network algorithm;The product is valid for case 2 waters, i.e. waters where inorganic particles dominate over phytoplankton (typically in coastal waters).

# setwd('G:/1 RM/3 gridded data/HERMES merged CHL 25km')
setwd("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES_monthly")
wd=getwd()

## HERMES chl1 product for class 1 ocean waters, inlcudes GSM, AV, and AVW files
files.chl1=list.files(wd, pattern=('_CHL1_MO'))
files.chl1.ave=files.chl1[grep(files.chl1, pattern=('_AV-'))] #average SeaWiFS only
files.chl1.avw=files.chl1[grep(files.chl1, pattern=('_AVW-'))] #weighted average merge of multiple satellite data
## Merge the AV and AVW chl1 product filenames
test=data.frame(files.chl1.ave,stringsAsFactors = FALSE); colnames(test)='chl1'
test2=data.frame(files.chl1.avw,stringsAsFactors = FALSE); colnames(test2)='chl1'
files.chl1.av=rbind(test, test2); rm(test); rm(test2)
## Select just GSM chl1 product filenames
files.chl1.gsm=data.frame(files.chl1[grep(files.chl1, pattern=('_GSM-'))]);colnames(files.chl1.gsm)='gsm' #GSM merge of multiple satellite data

## HERMES chl2 product files (coastal, limited data)
files.chl2=data.frame(list.files(wd, pattern=('_CHL2_MO')));colnames(files.chl2)='chl2' #HERMES chl2, only for MER and OLA sats (limited data)

## HERMES OC5 product filenames
files.oc5=data.frame(list.files(wd, pattern=('CHL-OC5_')));colnames(files.oc5)='oc5' #HERMES oc5 algorithm for coastal waters

### sort data lists to make sure it is in chronological order
test1=sort(files.chl1.av[,1])
test2=sort(files.chl1.gsm[,1])
test3=sort(files.oc5[,1])


### create raster stacks of data 
chl.av=nc2raster(files.chl1.av[,1], 'CHL1_mean')
chl.gsm=nc2raster(files.chl1.gsm[,1], 'CHL1_mean')
chl.oc5=nc2raster(files.oc5[,1], 'CHL1_mean')



#get date details from first nc file...
# MM=rep(NA, length(files.GSM))
# YY=rep(NA, length(files.GSM))
# # CHL=array(NA, dim=c())
# for (i in 1:length(files.GSM)){
#   nc1=nc_open(files.GSM[i])
#   v=strsplit(nc1$filename,split="_", fixed=TRUE)
#   dt=strsplit(v[[1]][2],split='-',fixed=T)[[1]][1]
#   MM[i]=as.numeric(substr(dt,5,6)) #month
#   YY[i]=as.numeric(substr(dt,1,4))
#   lon=ncvar_get(nc1, 'lon')
#   lat=ncvar_get(nc1, 'lat')
#   chl=ncvar_get(nc1, 'CHL1_mean')
#   # dim(chl)
#   colnames(chl)=lat
#   rownames(chl)=lon
#   nc_close(nc1)
#   # CHL[,,i]=chl
# }
# lat.val=as.numeric(lat)
# lon.val=as.numeric(colnames(CHL[,,1]))
# month.val=MM
# year.val=YY

## fake some year data to make old code work:
sdat=as.data.frame(matrix(NA, 41, 1))
sdat$year=seq(from=1977, to=2017, by=1)
# #continue...
datalab='Chl' #'CHLrange' #'Chl'
yr.lst=unique(year.val)
yrlist1=unique(sdat$year)
yy1=unique(year.val)
yy2=yrlist1[yrlist1 %in% yy1]
yrlist=data.frame(yy2, yy2)


## spring yearly means stacked raster - these are already monthly means, so just take mean for season and stack
mn.lst=c(2,3,4); season='Spring'
files.02GSM=list.files(wd, pattern='0201');files.02GSM=grep(files.02GSM, pattern='_GSM-', inv=F, value=T)
files.02GSM=files.02GSM[c(1:4,6:length(files.02GSM))] #drop 20020131-20020228
files.03GSM=list.files(wd, pattern='0301');files.03GSM=grep(files.03GSM, pattern='_GSM-', inv=F, value=T)
files.03GSM=files.03GSM[c(1:5,7:length(files.03GSM))] #drop 20030131-20030228
files.04GSM=list.files(wd, pattern='0401');files.04GSM=grep(files.04GSM, pattern='_GSM-', inv=F, value=T)
files.04GSM=files.04GSM[c(1:6,8:length(files.04GSM))] #drop 20040131-20040228
r1 <- raster(files.02GSM[1],  varname = "CHL1_mean")
r2 <- raster(files.03GSM[1],  varname = "CHL1_mean")
r3 <- raster(files.04GSM[1],  varname = "CHL1_mean")
s <- stack(r1, r2, r3)
x <- reclassify(s, cbind(0, NA))
# shp.dat= mean(x, na.rm=TRUE) # mean of pixels
# fun=function(x) { (range(x, na.rm=T)[2]-range(x, na.rm=T)[1])+1} # inclusive difference between high and low of range
fun=function(x) { if (is.na(x[1])){ NA } else {(range(x, na.rm=T)[2]-range(x, na.rm=T)[1])+1}}
shp.dat= calc(s, fun)
for (i in 2:length(files.02GSM)){
  r1 <- raster(files.02GSM[i],  varname = "CHL1_mean")
  r2 <- raster(files.03GSM[i],  varname = "CHL1_mean")
  r3 <- raster(files.04GSM[i],  varname = "CHL1_mean")
  s <- stack(r1, r2, r3)
  x <- reclassify(s, cbind(0, NA))
  # shp.dat2=mean(x, na.rm=TRUE) # mean of pixels
  shp.dat2= calc(s, fun) # function above (range)
  shp.dat=stack(shp.dat, shp.dat2)
}
if (season == 'Spring'){
  yrlist2=yrlist[-(1),] # drop 1997 for spring (no data)
}
# if (season == 'Fall'){
# yrlist=yrlist[-(19),] # drop 2015 for fall (data not downloaded yet)
# }
# save(shp.dat, file='Spr_25km_chl_gsm.rdata')
wd3= "G:/1 RM/2 Plankton Spatial Plots/data"
filename=paste(season, datalab, yrlist2[1,1],yrlist[length(yrlist[,1]),1],"shp.dat.rdata", sep="_")
mypath=file.path(wd3, filename)
save(shp.dat, file=mypath)










#load Chl from HERMES merged product - older set used for plankton analysis (see below for update)
# load("G:/1 RM/2 Plankton Spatial Plots/data/1997-2015 - chl/Spring_Chl_1998_2015_shp.dat.rdata")
# plot(shp.dat[[1]])
# lines(test@polygons[[1]]@Polygons[[1]]@coords, col='blue')
# 
# v2=list()
# for(i in 1:dim(shp.dat)[3]){
# v=extract(shp.dat[[i]], neus.shp)
# v1=lapply(v, function(x) mean(x, na.rm=T))
# v2[i]=list(v1)
# # names(v1[i])=seq(0:29)
# }
# m=matrix(unlist(v2), ncol=18, nrow=30) # box 0-29 =rows, years =cols
# 
# v1=extract(shp.dat[[1]], neus.shp)
# v1m=lapply(v1, function(x) mean(x, na.rm=T))
# 
# ### Spring time series trends by box (can use for calibration later)
# year=seq(1998, 2015, 1)
# for (i in 1:nrow(m)){
#   plot(m[i,]~year, type='l', ylim=c(0, 8), col='gray80', ylab='', xlab='')
#   par(new=T)
# }
# plot(colMeans(m, na.rm=T)~year, type='l', col='blue', ylim=c(0,8), lwd=2)
# par(new=F)

### load monthly mean 25 KM Chl 1997(9-12 only) through 2015(1-4 only) = 212 months (17 years + 8 months) -> (OLD)
#20170610 - added new data from HERMES GlobColour, now all monthly from 1998-2016 included
# load("G:/1 RM/3 gridded data/HERMES merged CHL 25km/monthly_chl_hermes_merged_gsm_25km.rdata")
setwd('G:/1 RM/3 gridded data/HERMES merged CHL 25km/1998_2016')
setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES merged CHL 25km/1998_2016')
wd=getwd()
files.AV=list.files(wd, pattern=('_AV-'))
files.GSM=list.files(wd, pattern=('_GSM-'))

files.01GSM=list.files(wd, pattern='0101');files.01GSM=grep(files.01GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.01GSM)
files.01GSM=files.01GSM[c(1:13,17:22)] #21
ii=data.frame(files.01GSM)
files.02GSM=list.files(wd, pattern='0201');files.02GSM=grep(files.02GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.02GSM)
files.02GSM=files.02GSM[c(1:4,6:20)] #drop 20020131-20020228
files.03GSM=list.files(wd, pattern='0301');files.03GSM=grep(files.03GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.03GSM)
files.03GSM=files.03GSM[c(1:5,7:20)] #drop 20030131-20030228
files.04GSM=list.files(wd, pattern='0401');files.04GSM=grep(files.04GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.04GSM)
files.04GSM=files.04GSM[c(1:6,8:20)] #drop 20040131-20040228
files.05GSM=list.files(wd, pattern='0501');files.05GSM=grep(files.05GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.05GSM)
files.05GSM=files.05GSM[c(1:7,9:20)] 
files.06GSM=list.files(wd, pattern='0601');files.06GSM=grep(files.06GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.06GSM)
files.06GSM=files.06GSM[c(1:8,10:20)] 
files.07GSM=list.files(wd, pattern='0701');files.07GSM=grep(files.07GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.07GSM)
files.07GSM=files.07GSM[c(1:9,11:20)] 
files.08GSM=list.files(wd, pattern='0801');files.08GSM=grep(files.08GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.08GSM)
files.08GSM=files.08GSM[c(1:10,12:20)] 
files.09GSM=list.files(wd, pattern='0901');files.09GSM=grep(files.09GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.09GSM)
files.09GSM=files.09GSM[c(1:11,13:20)] 
files.10GSM=list.files(wd, pattern='1001');files.10GSM=grep(files.10GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.10GSM)
files.10GSM=files.10GSM[c(1:12,14:20)] 
files.11GSM=list.files(wd, pattern='1101');files.11GSM=grep(files.11GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.11GSM)
files.11GSM=files.11GSM[c(1:13,15:20)] 
files.12GSM=list.files(wd, pattern='1201');files.12GSM=grep(files.12GSM, pattern='_GSM-', inv=F, value=T)
ii=data.frame(files.12GSM)
files.12GSM=files.12GSM[c(1:14,16:20)] 

### function to rasterize all files in a list, stack raster and return raster stack
### RM 20170609


# s=stack()
# for (i in 1:length(files.01GSM)){
#   r <- raster(files.01GSM[i],  varname = "CHL1_mean")
#   s=stack(s, r)
# }

### raster stacks of monthly data
r1=nc2raster(files.01GSM)
r2=nc2raster(files.02GSM)
r3=nc2raster(files.03GSM)
r4=nc2raster(files.04GSM)
r5=nc2raster(files.05GSM)
r6=nc2raster(files.06GSM)
r7=nc2raster(files.07GSM)
r8=nc2raster(files.08GSM)
r9=nc2raster(files.09GSM)
r10=nc2raster(files.10GSM)
r11=nc2raster(files.11GSM)
r12=nc2raster(files.12GSM)

### Function to extract data using a shapefile
extractMonths=function(x, shp){
  v2=list()
  for(i in 1:dim(x)[3]){
    v=extract(x[[i]], shp)
    v1=lapply(v, function(xx) mean(xx, na.rm=T))
    v2[i]=list(v1)
  }
  m=matrix(unlist(v2), ncol=dim(x)[3], nrow=30) # box 0-29 =rows, years =cols
  colnames(m)=seq(1998, 2016, by=1)
  rownamse(m)=
    return(m)
}

### returns Mean Chl per box (0-29, rows), by year from 1998-2016 (columns) for month indicated (mg chl a /m3)
Jan.chl=extractMonths(r1, neus.shp)
Feb.chl=extractMonths(r2, neus.shp)
Mar.chl=extractMonths(r3, neus.shp)
Apr.chl=extractMonths(r4, neus.shp)
May.chl=extractMonths(r5, neus.shp)
Jun.chl=extractMonths(r6, neus.shp)
Jul.chl=extractMonths(r7, neus.shp)
Aug.chl=extractMonths(r8, neus.shp)
Sep.chl=extractMonths(r9, neus.shp)
Oct.chl=extractMonths(r10, neus.shp)
Nov.chl=extractMonths(r11, neus.shp)
Dec.chl=extractMonths(r12, neus.shp)

# All.chl=nc2raster(files.GSM)

All.chl.av=nc2raster(files.chl1.av)
All.chl.gsm=nc2raster(files.chl1.gsm)
All.chl.oc5=nc2raster(files.chl1.oc5)

Chl.time=seq(ISOdate(1998,1,15), ISOdate(2016,12,15), "month") # monthly mean values 1998-2016
# dimnames(All.chl[,,3])=Chl.time

NEUSplotChlRaster=function(data, i, maxV){
  rasterX=data[[i]]
  col5=colorRampPalette(c('blue','white','red'))
  max_abolute_value=maxV #what is the maximum absolute value of raster?
  color_levels=20
  br <- seq(0, max_abolute_value, length.out=color_levels+1) 
  rng2=cellStats(rasterX, range)
  rng=c(0, maxV, rng2[2])
  arg=list(at=rng, labels=round(rng,2))
  plot(rasterX, col=col5(length(br) - 1), breaks=br,axis.args=arg, xlim=c(-77,-64),ylim=c(35,45),
       las=1, legend=F, main=Chl.time[[i]])
  map("worldHires", xlim=c(-77,-64),ylim=c(35,45), fill=T,border=0,col="gray", add=T)
  plot(rasterX, legend.only=T, col=col5(length(br) - 1),breaks=br,axis.args=arg, legend.shrink=0.5,
       smallplot=c(0.19,0.21, 0.6,0.80) )
}
NEUSplotChlRaster(All.chl, 10, 4) # data, choose date (1-228), max value

### Mg Chl a /m3 -> Mg Chl a per box (take mg chla/m3 * box area * chl depth (max 50m))
ts1=Jan.chl*bgm.z$chlZ*bgm.z$area
ts2=Feb.chl*bgm.z$chlZ*bgm.z$area
ts3=Mar.chl*bgm.z$chlZ*bgm.z$area
ts4=Apr.chl*bgm.z$chlZ*bgm.z$area
ts5=May.chl*bgm.z$chlZ*bgm.z$area
ts6=Jun.chl*bgm.z$chlZ*bgm.z$area
ts7=Jul.chl*bgm.z$chlZ*bgm.z$area
ts8=Aug.chl*bgm.z$chlZ*bgm.z$area
ts9=Sep.chl*bgm.z$chlZ*bgm.z$area
ts10=Oct.chl*bgm.z$chlZ*bgm.z$area
ts11=Nov.chl*bgm.z$chlZ*bgm.z$area
ts12=Dec.chl*bgm.z$chlZ*bgm.z$area
























