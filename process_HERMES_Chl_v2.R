# 
# url<- "ftp://ftp.hermes.acri.fr/345576024"
# filenames <- getURL(url, userpwd="ftp_hermes:hermes%", ftp.use.epsv = FALSE, dirlistonly = TRUE) #reading filenames from ftp-server

# OCCI https://esa-oceancolour-cci.org/
# https://rsg.pml.ac.uk/thredds/ncss/grid/CCI_ALL-v4.0-8DAY/dataset.html
# NCSS Request URL
# https://rsg.pml.ac.uk/thredds/ncss/grid/CCI_ALL-v4.0-8DAY/dataset.html
# /thredds/ncss/CCI_ALL-v4.0-8DAY
# https://rsg.pml.ac.uk/thredds/ncss/CCI_ALL-v4.0-8DAY?var=chlor_a&north=48&west=-80&east=-60&south=32&horizStride=1&time_start=1997-09-04T00%3A00%3A00Z&time_end=2018-12-27T00%3A00%3A00Z&timeStride=1&addLatLon=true

#get OCCI v4.2
gribfile='http://rsg.pml.ac.uk/thredds/ncss/CCI_ALL-v4.2-8DAY?var=chlor_a&var=chlor_a_log10_bias&var=chlor_a_log10_rmsd&north=48&west=-76&east=-64&south=35&disableProjSubset=on&horizStride=1&time_start=1997-09-04T00%3A00%3A00Z&time_end=2019-12-27T00%3A00%3A00Z&timeStride=1&addLatLon=true'
download.file(gribfile,'junk.nc',mode = "wb")


library(raster)
library(mgcv)
# library(sp)
library(maptools)
library(marmap)
library(maps)
library(mapdata)
library(rgeos)
library(ncdf4)
library(abind)
library(RColorBrewer)
library(chron)

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


### Function to plot a raster for the NES area and set scale
plotChlRaster=function(data, i, maxV, limit=F, dateval){
  rasterX=data[[i]]
    rng2=cellStats(rasterX, range)
  if (limit == 1){
    max_abolute_value=maxV #set limit manually with maxV input
    rng=c(0, max_abolute_value, rng2[2])
  }
  else {
    max_abolute_value=ceiling(rng2[2]) #round up actual max val
    rng=c(0, max_abolute_value)
  }
  color=rev(brewer.pal(11, "Spectral"))
  br <- seq(0, max_abolute_value, length.out=9) 
     arg=list(at=rng, labels=round(rng,1))
  plot(rasterX, col=color, breaks=br,axis.args=arg, xlim=c(-77,-64),ylim=c(35,45),
       las=1, legend=F, main=dateval[[i]])
  map("worldHires", xlim=c(-77,-64),ylim=c(35,45), fill=T,border=0,col="gray", add=T)
  plot(rasterX, legend.only=T, col=color,breaks=br,axis.args=arg, legend.shrink=0.5,
       smallplot=c(0.19,0.21, 0.6,0.80) )
}

# color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
#   scale = (length(lut)-1)/(max-min)
#   
#   dev.new(width=1.75, height=5)
#   plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#   axis(2, ticks, las=1)
#   for (i in 1:(length(lut)-1)) {
#     y = (i-1)/scale + min
#     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#   }
# }
# 
# plot_binned_strata=function(m, i, shp, maxV, limit=F){
#   rng2=cellStats(rasterX, range)
#   if (limit == 1){
#     max_abolute_value=maxV #set limit manually with maxV input
#     rng=c(0, max_abolute_value, rng2[2])
#   }
#   else {
#     max_abolute_value=ceiling(rng2[2]) #round up actual max val
#     rng=c(0, max_abolute_value)
#   }
#   color=rev(brewer.pal(11, "Spectral"))
#   br <- seq(0, max_abolute_value, length.out=length(color)) 
#   m.bin=cut(m, br)
#   arg=list(at=rng, labels=round(rng,1))
#   plot(shp, col=color[m.bin], breaks=br,axis.args=arg, xlim=c(-77,-64),ylim=c(35,45),
#        las=1, legend=F, main=av.dates[[i]])
#   map("worldHires", xlim=c(-77,-64),ylim=c(35,45), fill=T,border=0,col="gray", add=T)
#   color.bar(color, 0,max_abolute_value)
#   # plot(rasterX, legend.only=T, col=color,breaks=br,axis.args=arg, legend.shrink=0.5,
#   #      smallplot=c(0.19,0.21, 0.6,0.80) )
# }

# # testing...
# plotChlRaster(chl.av, 5, 5, limit=T)
# plotChlRaster(chl.av, 5, 20, limit=F)


# ## load shapefiles
# setwd("G:/1 RM/KINGSTON/transfer/shapefiles/epu_shapes")
setwd("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/KINGSTON/transfer/shapefiles/epu_shapes")
d2='/home/ryan/Desktop/shapefiles/epu_shapes'
wd=getwd()
gbk=rgdal::readOGR(paste(d2, '/', "EPU_GBKPoly.shp", sep=''))
gom=rgdal::readOGR(paste(d2, '/', "EPU_GOMPoly.shp", sep=''))
mab=rgdal::readOGR(paste(d2, '/', "EPU_MABPoly.shp", sep=''))
scs=rgdal::readOGR(paste(d2, '/', "EPU_SCSPoly.shp", sep=''))

#combine shapefiles GOM and GBK
gom.gbk.shp=gUnion(gom, gbk, byid=F, id=NULL)
gom.gbk.shp=gUnion(gom, gbk, byid=F, id=NULL)
gom.scs.shp=gUnion(gom, scs, byid=F, id=NULL)
mab.gbk.shp=gUnion(mab, gbk, byid=F, id=NULL)
NES.shp=gUnion(mab.gbk.shp, gom.scs.shp, byid=F, id=NULL)

setwd("C:/Users/ryan.morse/Desktop/NES_5area")
nes.five=rgdal::readOGR('nes_gbk_gome_gomw_mabn_mabsPoly.shp')


## grab ecomon strata and plot
setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/2 Plankton Spatial Plots/shapefiles')
ecomon.strata=rgdal::readOGR("EcoMon_strata.shp")
par(mar = c(0,0,0,0))
par(oma = c(0,0,0,0))
pdf(file='EcoMon_strata.pdf')
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
data(stateMapEnv)
map('state', fill = F, add=T) # add state lines
lines(ecomon.strata)
text(coordinates(ecomon.strata)[,1], coordinates(ecomon.strata)[,2], ecomon.strata$STRATA)
dev.off()

#### Load Chl data
# AVW: weighted average of single-sensor Level 2 CHL1 products
# GSM: GSM merging of single sensor L3 NRRS
# The CHL1 algorithms are applicable for "case 1" waters
# OC5 Gohin, F., 2011
# CHL2 is the chlorophyll concentration (mg/m3) for Case 2 waters (see section validity); L3 merge: AV; sensors: MER, OLA; Doerffer and Schiller (2007)
# CHL2 uses the a Neural Network algorithm;The product is valid for case 2 waters, i.e. waters where inorganic particles dominate over phytoplankton (typically in coastal waters).

### Get 8-day OCCI chlorophyll data
setwd('/media/ryan/Iomega_HDD/1 RM/3 gridded data/OCCI')
setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/OCCI')
setwd('H:/1 RM/3 gridded data/OCCI')


# nc1=nc_open('CCI_ALL-v4.0-8DAY.nc') # just chl
# nc1=nc_open('C:/Users/ryan.morse/Downloads/CCI_ALL-v4.0-8DAY.nc') # new file with error estimates
nc1=nc_open('C:/Users/ryan.morse/Documents/GitHub/JPSS/CCI_ALL-v4.2-8DAY.nc') #udpated with data fix 2019

lon.occi=ncvar_get(nc1, 'lon')
lat.occi=ncvar_get(nc1, 'lat')
chl.occi=ncvar_get(nc1, 'chlor_a')
chl.occi.bias=ncvar_get(nc1, 'chlor_a_log10_bias') # only if using new file in downloads (2GB)
chl.occi.rmsd=ncvar_get(nc1, 'chlor_a_log10_rmsd') # only if using new file in downloads (2GB)
time.occi=ncvar_get(nc1, 'time') #days since Jan 1, 1970
dim(chl.occi)
colnames(chl.occi)=lat.occi
rownames(chl.occi)=lon.occi
nc_close(nc1)

test=month.day.year(time.occi, c(1,1,1970)) # these are not 8-days apart.... something odd
occi.date=data.frame(test)
occi.date$F1=paste(occi.date[,3], occi.date[,1], occi.date[,2], sep='-')
occi.date$DOY=as.numeric(strftime(occi.date$F1, '%j'))
ddiff=diff(occi.date$DOY)
occi.date$diff=c(0, ddiff) # see OC-CCI manual in JPSS/calibration folder for list of missing dates, explains why some are not 8d

table(occi.date$year)
occi.date$week=c(seq(from=31, to=46, by=1),  rep(seq(from=1, to=46, by=1), 21))


# m=(t(chl.occi[,,3]))
# dimnames(m) <- list(lat=as.numeric(lat.occi), lon=as.numeric(lon.occi))
# t=raster(m)
# extent(t)=c(-80, -60, 32, 48)
# crs(t)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
# plot(t)

# occi=brick((chl.occi))
# extent(t)=c(-80, -60, 32, 48)
# crs(t)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

## Chl_loop over and stack
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

## Chl_RMSD loop over and stack
bb=c(-80, -60, 32, 48)
m2=t(chl.occi.rmsd[,,1])
# m2=t(m)#[ncol(m):1,] # flip and transpose matrix
occi.rmsd=raster(m2)
extent(occi.rmsd)=bb
crs(occi.rmsd)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
for(i in 2:dim(chl.occi.rmsd)[3]){
  m2=t(chl.occi.rmsd[,,i])
  # m2=t(m)#[ncol(m):1,] # flip and transpose matrix
  xx=raster(m2)
  extent(xx)=bb
  crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  # plot(xx)
  occi.rmsd=stack(occi.rmsd, xx)
  print(i)
}

## Chl_Bias loop over and stack
bb=c(-80, -60, 32, 48)
m2=t(chl.occi.bias[,,1])
# m2=t(m)#[ncol(m):1,] # flip and transpose matrix
occi.bias=raster(m2)
extent(occi.bias)=bb
crs(occi.bias)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
for(i in 2:dim(chl.occi.bias)[3]){
  m2=t(chl.occi.bias[,,i])
  # m2=t(m)#[ncol(m):1,] # flip and transpose matrix
  xx=raster(m2)
  extent(xx)=bb
  crs(xx)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  # plot(xx)
  occi.bias=stack(occi.bias, xx)
  print(i)
}



### HERMES 8-day composites ###
# setwd('G:/1 RM/3 gridded data/HERMES merged CHL 25km')
setwd("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES_8day")
setwd('/media/ryan/Iomega_HDD/1 RM/3 gridded data/HERMES_8day')
wd=getwd()
## HERMES chl1 product for class 1 ocean waters, inlcudes GSM, AV, and AVW files
files.chl1=list.files(wd, pattern=('_CHL1_8D'))
files.chl1.ave=files.chl1[grep(files.chl1, pattern=('_AV-'))] #average SeaWiFS only
files.chl1.avw=files.chl1[grep(files.chl1, pattern=('_AVW-'))] #weighted average merge of multiple satellite data
## Merge the AV and AVW chl1 product filenames
test=data.frame(files.chl1.ave,stringsAsFactors = FALSE); colnames(test)='chl1'
test2=data.frame(files.chl1.avw,stringsAsFactors = FALSE); colnames(test2)='chl1'
files.chl1.av=rbind(test, test2); rm(test); rm(test2)
## Select just GSM chl1 product filenames
files.chl1.gsm=data.frame(files.chl1[grep(files.chl1, pattern=('_GSM-'))],stringsAsFactors = FALSE);colnames(files.chl1.gsm)='gsm' #GSM merge of multiple satellite data

## HERMES chl2 product files (coastal, limited data)
files.chl2=list.files(wd, pattern=('_CHL2_8D')) #HERMES chl2, only for MER and OLA sats (limited data)
files.chl2=data.frame(files.chl2, stringsAsFactors = FALSE);colnames(files.chl2)='chl2'
## HERMES OC5 product filenames
files.oc5=list.files(wd, pattern=('CHL-OC5_')) #HERMES oc5 algorithm for coastal waters
files.oc5=data.frame(files.oc5, stringsAsFactors = F);colnames(files.oc5)='oc5' #HERMES oc5 algorithm for coastal waters

### sort data lists to make sure it is in chronological order
av.files=sort(files.chl1.av[,1])
gsm.files=sort(files.chl1.gsm[,1])
oc5.files=sort(files.oc5[,1])
chl2.files=sort(files.chl2[,1])

#get dates
av.dates.8d=list()
for (i in 1:length(av.files)){
  av.dates.8d[[i]]=strsplit(av.files[i],split="_", fixed=TRUE)[[1]][2]
}
test=do.call(rbind, strsplit(unlist(av.dates.8d), split="-", fixed=T))
dates.8d=matrix(NA, nrow=length(test[,1]), ncol=3)
dates.8d=data.frame(dates.8d)
dates.8d$X1=as.numeric(substr(test[,1], 1,4)) #YYYY
dates.8d$X2=substr(test[,1], 5,6) #MM
dates.8d$X3=substr(test[,1], 7,8) #DD
dates.8d$Y1=substr(test[,2], 1,4) #YYYY
dates.8d$Y2=substr(test[,2], 5,6) #MM
dates.8d$Y3=substr(test[,2], 7,8) #DD
dates.8d$F1=paste(dates.8d[,1], dates.8d[,2], dates.8d[,3], sep='-')
dates.8d$F2=paste(dates.8d[,4], dates.8d[,5], dates.8d[,6], sep='-')
dates.8d$DOY1=as.numeric(strftime(dates.8d$F1, '%j'))
dates.8d$DOY2=as.numeric(strftime(dates.8d$F2, '%j'))

## add week of year (should be 46 per year for 8-day files)
table(dates.8d$X1)
dates.8d$week=c(seq(from=31, to=46, by=1),  rep(seq(from=1, to=46, by=1), 21), seq(from=1, to=23, by=1))

### create raster stacks of data 
chl.av.8d=nc2raster(files.chl1.av[,1], 'CHL1_mean')
chl.gsm.8d=nc2raster(files.chl1.gsm[,1], 'CHL1_mean')
chl.oc5.8d=nc2raster(files.oc5[,1], 'CHL-OC5_mean')

#testing
plotChlRaster(chl.av.8d, 5, 5, limit=T, av.dates.8d)

### HERMES monthly composites ###
# setwd('G:/1 RM/3 gridded data/HERMES merged CHL 25km')
setwd("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES_monthly")
setwd('/media/ryan/Iomega_HDD/1 RM/3 gridded data/HERMES_monthly')
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
files.chl1.gsm=data.frame(files.chl1[grep(files.chl1, pattern=('_GSM-'))],stringsAsFactors = FALSE);colnames(files.chl1.gsm)='gsm' #GSM merge of multiple satellite data

## HERMES chl2 product files (coastal, limited data)
files.chl2=list.files(wd, pattern=('_CHL2_MO')) #HERMES chl2, only for MER and OLA sats (limited data)
files.chl2=data.frame(files.chl2, stringsAsFactors = FALSE);colnames(files.chl2)='chl2'
## HERMES OC5 product filenames
files.oc5=list.files(wd, pattern=('CHL-OC5_')) #HERMES oc5 algorithm for coastal waters
files.oc5=data.frame(files.oc5, stringsAsFactors = F);colnames(files.oc5)='oc5' #HERMES oc5 algorithm for coastal waters

### sort data lists to make sure it is in chronological order
av.files=sort(files.chl1.av[,1])
gsm.files=sort(files.chl1.gsm[,1])
oc5.files=sort(files.oc5[,1])
chl2.files=sort(files.chl2[,1])

#get dates
av.dates=list()
for (i in 1:length(av.files)){
  av.dates[[i]]=strsplit(av.files[i],split="_", fixed=TRUE)[[1]][2]
}

# nc1=nc_open(av.files[1])
# lon.av=ncvar_get(nc1, 'lon')
# lat.av=ncvar_get(nc1, 'lat')
# chl.av.1=ncvar_get(nc1, 'CHL1_mean')
# time.occi=ncvar_get(nc1, 'time') #days since Jan 1, 1970
# dim(chl.occi)

### create raster stacks of data 
chl.av=nc2raster(files.chl1.av[,1], 'CHL1_mean')
chl.gsm=nc2raster(files.chl1.gsm[,1], 'CHL1_mean')
chl.oc5=nc2raster(files.oc5[,1], 'CHL-OC5_mean')

nc1=nc_open(files.oc5[1,1])

i=20
plot(chl.av[[i]], main=test1[i])

# testing...
plotChlRaster(chl.gsm, 5, 5, limit=T, av.dates)
plotChlRaster(chl.av, 5, 5, limit=T, av.dates)
plotChlRaster(chl.av, 5, 20, limit=F, av.dates)

# create time series by box for ecomon strata
m.av=extract_calc(chl.av[[1:262]], ecomon.strata)
m.gsm=extract_calc(chl.gsm[[1:262]], ecomon.strata)
m.oc5=extract_calc(chl.oc5[[1:262]], ecomon.strata)

#CBay
plot(m.av[6,], type='l', main='Box 6 CBay')
lines(m.gsm[6,], col='red')
#GBK
plot(m.av[30,], type='l', main='Box 30 GBK')
lines(m.gsm[30,], col='red')
#NY
plot(m.av[17,], type='l', main='Box 17 Hudson River')
lines(m.gsm[17,], col='red')

plot((m.av[1,]-m.gsm[1,]), type='l')
abline(h=0, lty=2)
#
plot(m.av[26,], type='l', main='Box 26 GBK outer shelf')
lines(m.gsm[26,], col='red')
#
ii=6
plot(m.oc5[ii,], type='l', main=paste('Box',ii))
lines(m.gsm[ii,], col='red')
lines(m.av[ii,], col='blue')



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

# Chl.time=seq(ISOdate(1998,1,15), ISOdate(2016,12,15), "month") # monthly mean values 1998-2016

### read in calibration files from seaBASS
setwd('C:/Users/ryan.morse/Documents/GitHub/JPSS/calibration')
mrs=read.csv('1563559063629053_chlor_a.csv', skip=26, header=T, stringsAsFactors = F); mrs=mrs[-c(1:2),] # MERIS
swf=read.csv('1563558739820325_chlor_a.csv', skip=26, header=T, stringsAsFactors = F); swf=swf[-c(1:2),] # SeaWiFS
mds=read.csv('1563558871350863_chlor_a.csv', skip=26, header=T, stringsAsFactors = F); mds=mds[-c(1:2),] # Modis
vrs=read.csv('1563558999992368_chlor_a.csv', skip=26, header=T, stringsAsFactors = F); vrs=vrs[-c(1:2),] # Viirs-snpp


mrs$longitude=as.numeric(mrs$longitude); mrs$latitude=as.numeric(mrs$latitude); mrs$insitu_chlor_a=as.numeric(mrs$insitu_chlor_a)
coordinates(mrs)=~longitude+latitude#transform to Spatialpointsdataframe
proj4string(mrs)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #ensure same projection
pointsin=over(mrs, NES.shp) #find which boxes samples belong to
# pointsin=over(mrs, ecomon.strata)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(mrs)

swf$longitude=as.numeric(swf$longitude); swf$latitude=as.numeric(swf$latitude); swf$insitu_chlor_a=as.numeric(swf$insitu_chlor_a)
coordinates(swf)=~longitude+latitude#transform to Spatialpointsdataframe
proj4string(swf)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #ensure same projection
pointsin=over(swf, NES.shp) #find which boxes samples belong to
# pointsin=over(swf, ecomon.strata)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(swf)

mds$longitude=as.numeric(mds$longitude); mds$latitude=as.numeric(mds$latitude); mds$insitu_chlor_a=as.numeric(mds$insitu_chlor_a)
coordinates(mds)=~longitude+latitude#transform to Spatialpointsdataframe
proj4string(mds)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #ensure same projection
pointsin=over(mds, NES.shp) #find which boxes samples belong to
# pointsin=over(mds, ecomon.strata)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(mds)

vrs$longitude=as.numeric(vrs$longitude); vrs$latitude=as.numeric(vrs$latitude); vrs$insitu_chlor_a=as.numeric(vrs$insitu_chlor_a)
coordinates(vrs)=~longitude+latitude#transform to Spatialpointsdataframe
proj4string(vrs)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #ensure same projection
pointsin=over(vrs, NES.shp) #find which boxes samples belong to
# pointsin=over(vrs, ecomon.strata)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(vrs)


### get calibration Chl from WOD
d2='C:/Users/ryan.morse/Documents/GitHub/JPSS/calibration/WOD'
d2='~/Git/JPSS/calibration/WOD'
d2='/home/ryan/Git/JPSS_2/JPSS/calibration/WOD'
ncfiles=list.files(path=d2, pattern='.nc')
nc.str=strsplit(ncfiles, '.nc')

i=1 #CTD hi res
nc1=nc_open(paste(d2, '/', ncfiles[i], sep=''))
# nc1=nc_open(ncfiles[i])
wod.chl=ncvar_get(nc1, 'Chlorophyll')
wod.lat=ncvar_get(nc1, 'lat')
wod.lon=ncvar_get(nc1, 'lon')
wod.time=ncvar_get(nc1, 'time') #days since Jan 1, 1770
wod.z=ncvar_get(nc1, 'z')
nc1$var$time$units
# i=4 #profiler buoy samples
# nc1=nc_open(ncfiles[i])
# wod.pfl.chl=ncvar_get(nc1, 'Chlorophyll')
# wod.pfl.lat=ncvar_get(nc1, 'lat')
# wod.pfl.lon=ncvar_get(nc1, 'lon')
# wod.pfl.time=ncvar_get(nc1, 'time') #days since Jan 1, 1770
# i=2 #glider samples
# nc1=nc_open(ncfiles[i])
# wod.gld.chl=ncvar_get(nc1, 'Chlorophyll')
# wod.gld.lat=ncvar_get(nc1, 'lat')
# wod.gld.lon=ncvar_get(nc1, 'lon')
# wod.gld.time=ncvar_get(nc1, 'time') #days since Jan 1, 1770
# i=3 #bottle samples
# nc1=nc_open(ncfiles[i])
# wod.osd.chl=ncvar_get(nc1, 'Chlorophyll')
# wod.osd.lat=ncvar_get(nc1, 'lat')
# wod.osd.lon=ncvar_get(nc1, 'lon')
# wod.osd.time=ncvar_get(nc1, 'time') #days since Jan 1, 1770

## WOD CTD samples, surface depth
wod.chl.df=data.frame(wod.chl[which(wod.z==0)], wod.lon, wod.lat, wod.z[which(wod.z==0)], wod.time)
test=month.day.year(wod.chl.df$wod.time, c(1,1,1770))
wod.chl.df$month=test$month
wod.chl.df$day=test$day
wod.chl.df$year=test$year
wod.chl.df2=wod.chl.df[which(wod.chl.df$year>1996),]
colnames(wod.chl.df2)=c('chl', 'lon', 'lat', 'z', 'jday', 'month', 'day', 'year')
wod.chl.df2$F1=paste(wod.chl.df2$year, wod.chl.df2$month, wod.chl.df2$day, sep='-')
wod.chl.df2$DOY=as.numeric(strftime(wod.chl.df2$F1, '%j'))
wod.chl.df2=wod.chl.df2[complete.cases(wod.chl.df2),]

save(wod.chl.df2, file='/home/ryan/Git/JPSS_2/JPSS/WOD_CTD_surface_chl.RData')

barplot(table(round(wod.chl.df2$chl, digits=1)))
barplot(table(wod.chl.df2$month))
barplot(table(wod.chl.df2$year))

# coordinates(wod.chl.df2)=~lon+lat #transform to Spatialpointsdataframe
# proj4string(wod.chl.df2)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #ensure same projection
# pointsin=over(wod.chl.df2, NES.shp) #find which boxes samples belong to
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(wod.chl.df2$lon, wod.chl.df2$lat)
points(wod.chl.df2)

### build initial list of time matchups - indicates the dimension of the chl raster stack to extract data from
# yy=unique(wod.chl.df2$year) # unique years
wod.chl.df2$smatch=NA
# wod.chl.df2$DOYmed=NA
wod.chl.df2$sDOY1=NA
wod.chl.df2$sDOY2=NA
for(i in 1:length(wod.chl.df2$lon)){
  ylim=which(dates.8d$X1==wod.chl.df2$year[i])
  xmn=which(dates.8d$DOY1[ylim]<=wod.chl.df2$DOY[i])#[dates.8d$Y1==yj]
  xmx=which(dates.8d$DOY2[ylim]>=wod.chl.df2$DOY[i])#[dates.8d$Y1==yj]
  both=ylim[which(xmn%in%xmx)]
  if(length(both)<1){
    next
  }
  else {
    wod.chl.df2$smatch[i]=both
    wod.chl.df2$sDOY1[i]=dates.8d$DOY1[both]
    wod.chl.df2$sDOY2[i]=dates.8d$DOY2[both]
  }
}




wod.chl.df2$DOYmed=round((wod.chl.df2$sDOY1+wod.chl.df2$sDOY2)/2, digits=0) # median satellite DOY
wod.chl.df2$ddif=wod.chl.df2$DOY-wod.chl.df2$DOYmed # difference from median satellite date


### now extract from raster
# i=2746
# tt=extract(chl.gsm.8d[[140]], wod.chl.df2[i,], method='bilinear', fun='mean', na.rm=T)

wod.chl.df2=wod.chl.df2[complete.cases(wod.chl.df2$smatch),]
WOD=wod.chl.df2
wod.chl.df2$gsm=NA
wod.chl.df2$oc5=NA
wod.chl.df2$av=NA
wod.chl.df2$occi=NA
wod.chl.df2$occiv5=NA

coordinates(wod.chl.df2)=~lon+lat #transform to Spatialpointsdataframe
### method=bilinear #interpolates value from 4 nearest raster cells #simple is for cell only
for(i in 1:length(wod.chl.df2$chl)){
  wod.chl.df2$gsm[i]=raster::extract(chl.gsm.8d[[wod.chl.df2$smatch[i]]], wod.chl.df2[i,], method='bilinear', fun='mean', na.rm=T)
  # wod.chl.df2$av[i]=extract(chl.av.8d[[wod.chl.df2$smatch[i]]], wod.chl.df2[i,], method='bilinear', fun='mean', na.rm=T)
  # wod.chl.df2$oc5[i]=extract(chl.oc5.8d[[wod.chl.df2$smatch[i]]], wod.chl.df2[i,], method='bilinear', fun='mean', na.rm=T)
  wod.chl.df2$occi[i]=raster::extract(occi[[wod.chl.df2$smatch[i]]], wod.chl.df2[i,], method='bilinear', fun='mean', na.rm=T)
  wod.chl.df2$occiv5[i]=raster::extract(occi.v5[[wod.chl.df2$smatch[i]]], wod.chl.df2[i,], method='bilinear', fun='mean', na.rm=T)
  if (i%%100==0){
  }
}

  
## add time to dataframe for satellite matchup
wod.chl.df2$jtime=wod.chl.df2$jday-floor(wod.chl.df2$jday)
wod.chl.df2$time=format(times(wod.chl.df2$jtime))
write.csv(wod.chl.df2, file='WOD_calibration_to_surface_chlorophyll_8day.csv')
write.csv(wod.chl.df2, file='WOD_surf.csv')



WOD$gsm=wod.chl.df2$gsm
WOD$av=wod.chl.df2$av
WOD$oc5=wod.chl.df2$oc5
WOD=WOD[complete.cases(WOD$gsm),]
WOD=WOD[complete.cases(WOD$chl),]


colorpal=viridis::viridis(8)
# plot(log10(WOD$chl)~log10(WOD$gsm), type='n')#, color=colorpal[WOD$ddif+4])
# # points(log10(WOD$chl),log10(WOD$gsm), type='p', col=colorpal[WOD$ddif+4])
# points(log10(WOD$gsm), log10(WOD$chl), type='p', col=colorpal[WOD$ddif+4])
# abline(0,1)

# Non log transformed data
plot(WOD$chl~WOD$oc5, type='n')#, color=colorpal[WOD$ddif+4])
points(WOD$oc5, WOD$chl, type='p', col=colorpal[WOD$ddif+4])
abline(0,1)
reg1=lm(WOD$chl~WOD$oc5)
summary(reg1)
abline(reg1$coefficients[1], reg1$coefficients[2], col='red')

# plot(x=WOD$gsm, y=WOD$chl, log='xy') #, color=colorpal[WOD$ddif+4])
# abline(0,1)
# abline(reg1$coefficients[1], reg1$coefficients[2], col='blue')
# barplot(table(WOD$ddif))

## log transform both X and Y
x=log10(WOD$gsm+0.001)
y=log10(WOD$chl+0.001)
reg1=lm(y~x)
summary(reg1)
xy=data.frame(x,y)
xy=xy[complete.cases(x),]
plot(y~x, type='n')#, log='xy')#, color=colorpal[WOD$ddif+4])
points(log10(WOD$gsm), log10(WOD$chl), type='p', col=colorpal[WOD$ddif+4])
abline(0,1)
abline(reg1$coefficients[1], reg1$coefficients[2], col='red')


library(plotrix)
taylor.diagram(wod.chl.df2$chl, wod.chl.df2$oc5, col='red')
taylor.diagram(wod.chl.df2$chl, wod.chl.df2$gsm, add=T, col='blue')
taylor.diagram(wod.chl.df2$chl, wod.chl.df2$av, add=T, col='green')
taylor.diagram(wod.chl.df2$chl, wod.chl.df2$occi, add=T, col='black')
taylor.diagram(wod.chl.df2$chl, wod.chl.df2$occiv5, add=T, col='red')


#subset to low or high chl
test=wod.chl.df2[which(wod.chl.df2$chl<1),]
taylor.diagram(test$chl, test$oc5, col='red')
taylor.diagram(test$chl, test$gsm, add=T, col='blue')
taylor.diagram(test$chl, test$av, add=T, col='green')
taylor.diagram(test$chl, test$occi, add=T, col='black')
taylor.diagram(test$chl, test$occiv5, add=T, col='red')


RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}
# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}
# Calculate error
error <- actual - predicted


test=WOD[complete.cases(WOD$occciv5),]
RMSE(o=test$chl, m=test$occciv5)
dim(test)
test=WOD[complete.cases(WOD$occciv4),]
RMSE(o=test$chl, m=test$occciv4)
dim(test)
test=WOD[complete.cases(WOD$gsm),]
RMSE(o=test$chl, m=test$gsm)
dim(test)

test=WOD[complete.cases(WOD),] # all vars, equal lenght (much data lost, all remaining in NES)
RMSE(o=test$chl, m=test$occciv5)
RMSE(o=test$chl, m=test$occciv4)
RMSE(o=test$chl, m=test$gsm)


### look at residuals
test$rv4=test$chl-test$occi
test$rv5=test$chl-test$occiv5
test$rgsm=test$chl-test$gsm

mae(test$rgsm)
mae(test$rv4)
mae(test$rv5)
rmse(test$rgsm)
rmse(test$rv4)
rmse(test$rv5)


plot(test$rv4,type='p', ylim=c(-10, 10), main='v4')
plot(test$rv5,type='p', ylim=c(-10, 10), main='v5')
plot(test$rgsm,type='p', ylim=c(-10, 10), main='gsm')
plot(test$rv4~test$occciv4, type='p', main='v4 resid', ylim=c(-5, 15), xlim=c(0, 8))
plot(test$rv5~test$occciv5, type='p', main='v5 resid', ylim=c(-5, 15), xlim=c(0, 8))
plot(test$rgsm~test$gsm, type='p', main='gsm resid', ylim=c(-5, 15), xlim=c(0, 8))
## log scale
plot(test$occciv4 ~ test$chl,type='p', main='v4', log='xy'); abline(a=0, b=1)
plot(test$occciv5 ~ test$chl,type='p', main='v5', log='xy'); abline(a=0, b=1)
plot(test$gsm ~ test$chl,type='p', main='gsm', log='xy'); abline(a=0, b=1)

plot(log(test$occciv4+1) ~ log(test$chl+1) ,type='p', ylim=c(0, 2), xlim=c(0,2), main='v4'); abline(a=0, b=1)
plot(log(test$occciv5+1) ~log(test$chl+1),type='p', ylim=c(0, 2), xlim=c(0,2), main='v5'); abline(a=0, b=1)
plot(log(test$gsm+1) ~log(test$chl+1) ,type='p', ylim=c(0, 2), xlim=c(0,2), main='gsm'); abline(a=0, b=1)


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

tr1=sbmodisa %>% select( latitude, longitude, date_time, insitu_chlor_a, aqua_chlor_a)
tr2=sbmodist %>% select( latitude, longitude, date_time, insitu_chlor_a, terra_chlor_a)
tr3=sbmeris %>% select( latitude, longitude, date_time, insitu_chlor_a, meris_chlor_a)
tr4=sbseawifs %>% select( latitude, longitude, date_time, insitu_chlor_a, seawifs_chlor_a)
tr5=sbviirssnpp %>% select( latitude, longitude, date_time, insitu_chlor_a, viirs_chlor_a)
# tr=tr1 %>% left_join(tr2, by=c('latitude', 'longitude', 'date_time'))

tr=tr1 %>% full_join(tr2, by=c('latitude', 'longitude', 'date_time'), name=NULL)
tr=tr%>% full_join(tr3, by=c('latitude', 'longitude', 'date_time'), name=NULL)
tr=tr%>% full_join(tr4, by=c('latitude', 'longitude', 'date_time'), name=NULL)
tr=tr%>% full_join(tr5, by=c('latitude', 'longitude', 'date_time'), name=NULL)
### merge all in-situ chlorophylls together, drop merged repeats
insitu_chl=tr %>% select(insitu_chlor_a, insitu_chlor_a.x, insitu_chlor_a.y, insitu_chlor_a.x.x, insitu_chlor_a.y.y) %>%
  mutate(ischl=rowMeans(., na.rm=T))
tr$insitu_chl_final=insitu_chl$ischl
tr=tr %>% select(-insitu_chlor_a, -insitu_chlor_a.x, -insitu_chlor_a.y, -insitu_chlor_a.x.x, -insitu_chlor_a.y.y)

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
satv5chl=satshpv5 %>% tidyr::separate(SAT_CHLOR_A_CCI, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% select(v1:v9)
# satv5chl$gm=apply(satv5chl[,1:9], 1, FUN=gm_mean, na.rm=T)
satv5chl$med=apply(satv5chl[,1:9], 1, median, na.rm=T)
satshpv5$SAT_MED_CHLOR=satv5chl$med
tv5=satshpv5[complete.cases(satshpv5$SAT_MED_CHLOR),]
### now do for v4 OCCCI data
satv4chl=satshpv4 %>% tidyr::separate(SAT_CHLOR_A_OCI, sep=";", into=paste("v", 1:9, sep=''), convert=T) %>% select(v1:v9)
satv4chl$med=apply(satv4chl[,1:9], 1, median, na.rm=T)
satshpv4$SAT_MED_CHLOR=satv4chl$med
tv4=satshpv4[complete.cases(satshpv4$SAT_MED_CHLOR),]


### read sat-ship matchup file from Kyle turner
ktsatmat=readxl::read_excel('/home/ryan/Downloads/nes_insitu_satellite_chl_upper10m.xlsx', na="NaN", 
col_types=c('numeric', 'numeric', 'guess', 'numeric', 'text', 'numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
# ktsatmat=ktsatmat %>% mutate_at(vars(`viirs_snpp_chlor_a [mg m-3]`, `snpp_time_dif [hours]`, `viirs_noaa20_chlor_a [mg m-3]`, `noaa20_time_dif [hours]`,
                            # `modis_aqua_chlor_a [mg m-3]`, `modis_time_dif [hours]`, `occci_chlor_a [mg m-3]`), na_if, "NaN") 
# ktsatmat[ktsatmat=="NaN"]=NA #as.numeric(ktsatmat)
t=ktsatmat[complete.cases(ktsatmat$`occci_chlor_a [mg m-3]`),]
plot(log(t$`in_situ_chl [mg m-3]`)~log(t$`occci_chlor_a [mg m-3]`), type='p')
abline(a=0, b=1)
RMSE(o=t$`in_situ_chl [mg m-3]` , m=t$`occci_chlor_a [mg m-3]`)

t2=satshpv5 %>% filter(CRUISE!='SeaBASS')
t3=satshpv5 %>% filter(CRUISE=='SeaBASS')
table(year(t3$SHIP_DATE))
### Turner merge v5 sat-matchups from Kim
ktest5=left_join(tv5,ktsatmat, by=c("date"="datetime", "SHIP_LAT"="lat", "SHIP_LON"="lon"))
ktest50=ktest5[complete.cases(ktest5$`in_situ_chl [mg m-3]`),]
plot(log(ktest50$`in_situ_chl [mg m-3]`) ~log(ktest50$SAT_MED_CHLOR), type='p', main="OCCCI v5.0 median vs in situ chl", xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(ktest50$`in_situ_chl [mg m-3]`) , m=log(ktest50$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(ktest50)[1],sep=''))
### Turner merge v4 sat-matchups from Kim
ktest4=left_join(tv4,ktsatmat, by=c("date"="datetime", "SHIP_LAT"="lat", "SHIP_LON"="lon"))
ktest40=ktest4[complete.cases(ktest4$`in_situ_chl [mg m-3]`),]
plot(log(ktest40$`in_situ_chl [mg m-3]`) ~log(ktest40$SAT_MED_CHLOR), type='p', main="OCCCI v4.2 median vs in situ chl", xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(ktest40$`in_situ_chl [mg m-3]`) , m=log(ktest40$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(ktest40)[1],sep=''))

### now add in tr to see if any additional stations are added... YES
test=full_join(ktest5,tr, by=c("date", "SHIP_LAT"="latitude", "SHIP_LON"="longitude"))
test$insituchl=apply(test[,c(29,43)], 1, mean, na.rm=T)
# sum(is.na(test$`in_situ_chl [mg m-3]`)) #32637
# sum(is.na(test$insitu_chl_final)) #32854
# sum(is.na(test$insituchl)) #32450
test5=test[complete.cases(test$insituchl),] # 916 using median of 3x3 pixels (na.rm=T)
plot(log(test5$insituchl) ~log(test5$SAT_MED_CHLOR), type='p', main="OCCCI v5.0 median vs in situ chl", xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(test5$insituchl) , m=log(test5$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(test5)[1],sep=''))

test=left_join(ktest4,tr, by=c("date", "SHIP_LAT"="latitude", "SHIP_LON"="longitude"))
test$insituchl=apply(test[,c(28,42)], 1, mean, na.rm=T)
sum(is.na(test$`in_situ_chl [mg m-3]`)) #32637
sum(is.na(test$insitu_chl_final)) #32854
sum(is.na(test$insituchl)) #32450
test5=test[complete.cases(test$insituchl),] # 916 using median of 3x3 pixels (na.rm=T)
plot(log(test5$insituchl) ~log(test5$SAT_MED_CHLOR), type='p', main="OCCCI v5.0 median vs in situ chl", xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(test5$insituchl) , m=log(test5$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(test5)[1],sep=''))



### compare time stamps -> need to format for comparison
tr$date_time[1] %in% tv5$SHIP_DATE # not formatted properly
a=tr$date_time[1] #"2009-10-19 15:25:00"
b=satshpv5$SHIP_DATE[40149] #"2009-10-19 15:25:00 UTC"
as.Date(lubridate::intersect(a,b), origin = "1970-01-01") #still not formatted properly
as.Date(lubridate::intersect(as.Date(a),as.Date(b)), origin = "1970-01-01") # works
as.Date(tr$date_time[1]) %in% as.Date(tv5$SHIP_DATE) # works
test=as.Date(tr$date_time) %in% as.Date(tv5$SHIP_DATE); sum(test) # works, number of in situ samples in OCCCI v5 database from Kim
test4=as.Date(tr$date_time) %in% as.Date(tv4$SHIP_DATE); sum(test)

### need POSIXct for times to be recognized
as.POSIXct(tr$date_time[1]) %in% as.POSIXct(tv5$SHIP_DATE) # does not work
as.POSIXct(tr$date_time[1], tz="UTC") %in% as.POSIXct(tv5$SHIP_DATE, tz="UTC") # works

### merge on dates ---> DOES NOT WORK
tr$date=as.Date(tr$date_time, origin = "1970-01-01")
tv5$date=as.Date(tv5$SHIP_DATE, origin = "1970-01-01")
test=left_join(tr, tv5, by="date") # WAY TOO LONG
### USING POSIXct with times ---> WORKS
tr$date=as.POSIXct(tr$date_time, tz="UTC") #origin = "1970-01-01")
tv5$date=as.POSIXct(tv5$SHIP_DATE, tz="UTC") #origin = "1970-01-01")
tv4$date=as.POSIXct(tv4$SHIP_DATE, tz="UTC") #origin = "1970-01-01")
### NOW do merge:
### v5 data using median of 3x3 pixels 
test=left_join(tr, tv5, by=c("date", "latitude"="SHIP_LAT", "longitude"="SHIP_LON")) # 2605 long
test5=test[complete.cases(test$SAT_MED_CHLOR),] # 916 using median of 3x3 pixels (na.rm=T)
plot(log(test5$insitu_chl_final) ~log(test5$SAT_MED_CHLOR), type='p', main="OCCCI v5.0 median vs in situ chl", xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(test5$insitu_chl_final) , m=log(test5$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(test5)[1],sep=''))
### v4 data
test=left_join(tr, tv4, by=c("date", "latitude"="SHIP_LAT", "longitude"="SHIP_LON")) # 2609 long
test4=test[complete.cases(test$SAT_MED_CHLOR),] # 917 using median of 3x3 pixels (na.rm=T)
plot(log(test4$insitu_chl_final) ~log(test4$SAT_MED_CHLOR), type='p', main="OCCCI v4.2 median vs in situ chl",xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(test4$insitu_chl_final) , m=log(test4$SAT_MED_CHLOR))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(test4)[1],sep=''))
### plot non-log transformed
plot(test4$SAT_MED_CHLOR ~ test4$insitu_chl_final, type='p', main="OCCCI v4.2 median vs in situ chl", ylim=c(0,10), xlim=c(0,20)); abline(a=0, b=1)
a=RMSE(o=test4$insitu_chl_final , m=test4$SAT_MED_CHLOR)
text(18,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(18,2, paste('N: ', dim(test4)[1],sep=''))
plot(test5$SAT_MED_CHLOR ~ test5$insitu_chl_final, type='p', main="OCCCI v5.0 median vs in situ chl", ylim=c(0,10), xlim=c(0,20)); abline(a=0, b=1)
a=RMSE(o=test5$insitu_chl_final , m=test5$SAT_MED_CHLOR)
text(18,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(18,2, paste('N: ', dim(test5)[1],sep=''))

### Plot satellite center values log and non log-transformed
## clean Inf values...
test44=test4
test44$SAT_CENTER_CHLOR_A_OCI[which(is.infinite(test44$SAT_CENTER_CHLOR_A_OCI))]=NA
test44=test44[complete.cases(test44$SAT_CENTER_CHLOR_A_OCI),]
test55=test5
test55$SAT_CENTER_CHLOR_A_CCI[which(is.infinite(test55$SAT_CENTER_CHLOR_A_CCI))]=NA
test55=test55[complete.cases(test55$SAT_CENTER_CHLOR_A_CCI),]
plot(test44$SAT_CENTER_CHLOR_A_OCI ~ test44$insitu_chl_final, type='p', main="OCCCI v4.2 center vs in situ chl", ylim=c(0,10), xlim=c(0,20)); abline(a=0, b=1)
a=RMSE(o=test44$insitu_chl_final , m=test44$SAT_CENTER_CHLOR_A_OCI)
text(18,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(18,2, paste('N: ', dim(test44)[1],sep=''))
plot(test55$SAT_CENTER_CHLOR_A_CCI ~ test55$insitu_chl_final, type='p', main="OCCCI v5.0 center vs in situ chl", ylim=c(0,10), xlim=c(0,20)); abline(a=0, b=1)
a=RMSE(o=test55$insitu_chl_final , m=test55$SAT_CENTER_CHLOR_A_CCI)
text(18,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(18,2, paste('N: ', dim(test55)[1],sep=''))

## log transformed
plot(log(test44$insitu_chl_final) ~log(test44$SAT_CENTER_CHLOR_A_OCI), type='p', main="OCCCI v4.2 center vs in situ chl",xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(test44$insitu_chl_final) , m=log(test44$SAT_CENTER_CHLOR_A_OCI))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(test44)[1],sep=''))
plot(log(test55$insitu_chl_final) ~log(test55$SAT_CENTER_CHLOR_A_CCI), type='p', main="OCCCI v5 center vs in situ chl",xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(test55$insitu_chl_final) , m=log(test55$SAT_CENTER_CHLOR_A_CCI))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(test55)[1],sep=''))



## log transformed
plot(log(t44$insitu_chl_final) ~log(t44$SAT_CENTER_CHLOR_A_OCI), type='p', main="OCCCI v4.2 center vs in situ chl",xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(t44$insitu_chl_final) , m=log(t44$SAT_CENTER_CHLOR_A_OCI))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(t44)[1],sep=''))

plot(log(t55$insitu_chl_final) ~log(t55$SAT_CENTER_CHLOR_A_CCI), type='p', main="OCCCI v5 center vs in situ chl",xlim=c(-3,3), ylim=c(-4,4))
abline(a=0, b=1)
a=RMSE(o=log(t55$insitu_chl_final) , m=log(t55$SAT_CENTER_CHLOR_A_CCI))
text(-2,3, paste('RMSE: ', round(a, digits=2),sep=''))
text(-2,2, paste('N: ', dim(t55)[1],sep=''))



#extract just lat/lons for lines
gbk.lonlat =as.data.frame(lapply(slot(gbk, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
gom.lonlat =as.data.frame(lapply(slot(gom, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
mab.lonlat =as.data.frame(lapply(slot(mab, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
scs.lonlat =as.data.frame(lapply(slot(scs, "polygons"), function(x) lapply(slot(x, "Polygons"), function(y) slot(y, "coords"))))
# nes.lonlat =as.data.frame(lapply(slots(), function))
gom.mat=as.matrix(gom.lonlat)
gbk.mat=as.matrix(gbk.lonlat)
mab.mat=as.matrix(mab.lonlat)
scs.mat=as.matrix(scs.lonlat)
m4=as.matrix(wod.chl.df2@coords)
wod.chl.df2$epu=NA
wod.chl.df2$epu[which(in.out(gbk.mat, m4))]='GBK'
wod.chl.df2$epu[which(in.out(gom.mat, m4))]='GOM'
wod.chl.df2$epu[which(in.out(scs.mat, m4))]='SCS'
wod.chl.df2$epu[which(in.out(mab.mat, m4))]='MAB'

#subset to region
limitc='GOM' #GBK SCS MAB
# limitc='All EPU' #GBK SCS MAB
test=wod.chl.df2[which(wod.chl.df2$epu==limitc),]
# test=wod.chl.df2[which(wod.chl.df2$epu=='GOM' | wod.chl.df2$epu=='GBK' | wod.chl.df2$epu=='MAB'),]
taylor.diagram(test$chl, test$oc5, col='red', main=paste(limitc, ' only;',' n=',length(test), sep=''), pos.cor = F)
taylor.diagram(test$chl, test$gsm, add=T, col='blue',pos.cor = F)
taylor.diagram(test$chl, test$av, add=T, col='green',pos.cor = F)
taylor.diagram(test$chl, test$occi, add=T, col='black',pos.cor = F)
barplot(table(test$month), main=limitc)
barplot(table(test$year), main=limitc)

plot(log(test$chl)~log(test$gsm), type='p')
plot(log(test$chl)~log(test$occi), type='p')


map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(test@coords, pch=19)

limitc='NA' #GBK SCS MAB
test=wod.chl.df2[is.na(wod.chl.df2$epu),]
taylor.diagram(test$chl, test$oc5, col='red', main=paste(limitc, ' only;',' n=',length(test), sep=''),pos.cor = F)
taylor.diagram(test$chl, test$gsm, add=T, col='blue')
taylor.diagram(test$chl, test$av, add=T, col='green')
taylor.diagram(test$chl, test$occi, add=T, col='black')

#extract data for shapes:
w.occi.nes=extract_calc(occi[[1:1028]], NES.shp)
occi.date$NESchl=w.occi.nes[1,]
w.occi.gom=extract_calc(occi[[1:1028]], gom)
occi.date$GOMchl=w.occi.gom[1,]
w.occi.gbk=extract_calc(occi[[1:1028]], gbk)
occi.date$GBKchl=w.occi.gbk[1,]
w.occi.mab=extract_calc(occi[[1:1028]], mab)
occi.date$MABchl=w.occi.mab[1,]

w.occi.nes.lg=extract_calc(lg.occi[[1:1028]], NES.shp)
occi.date$NESchllg=w.occi.nes.lg[1,]
w.occi.gom=extract_calc(occi[[1:1028]], gom)
occi.date$GOMchl=w.occi.gom[1,]
w.occi.gbk=extract_calc(occi[[1:1028]], gbk)
occi.date$GBKchl=w.occi.gbk[1,]
w.occi.mab=extract_calc(occi[[1:1028]], mab)
occi.date$MABchl=w.occi.mab[1,]


w.occi.nes.rmsd=extract_calc(occi.rmsd[[1:1028]], NES.shp) #NES
occi.date$NESchlrmsd=w.occi.nes.rmsd[1,]
w.occi.gom.rmsd=extract_calc(occi.rmsd[[1:1028]], gom)
occi.date$GOMchlrmsd=w.occi.gom.rmsd[1,]
w.occi.gbk.rmsd=extract_calc(occi.rmsd[[1:1028]], gbk)
occi.date$GBKchlrmsd=w.occi.gbk.rmsd[1,]
w.occi.mab.rmsd=extract_calc(occi.rmsd[[1:1028]], mab)
occi.date$MABchlrmsd=w.occi.mab.rmsd[1,]

library(dplyr)
tt=select(occi.date, month, year, NESchl) %>% group_by(year) %>% summarise(mean=mean(NESchl))
plot(tt, type='b', main=paste("NES annual Chl"))
tt.sp=select(occi.date, month, year, NESchl) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(NESchl))
plot(tt.sp, type='b',main=paste("NES Jan-Jun Chl"))
tt.fl=select(occi.date, month, year, NESchl) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(NESchl))
plot(tt.fl, type='b', main=paste("NES Jul-Dec Chl"))

tt=select(occi.date, month, year, GBKchl) %>% group_by(year) %>% summarise(mean=mean(GBKchl))
plot(tt, type='b', main=paste("GBK annual Chl"))
tt.sp=select(occi.date, month, year, GBKchl) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(GBKchl))
plot(tt.sp, type='b',main=paste("GBK Jan-Jun Chl"))
tt.fl=select(occi.date, month, year, GBKchl) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(GBKchl))
plot(tt.fl, type='b', main=paste("GBKchl Jul-Dec Chl"))

tt=select(occi.date, month, year, GOMchl) %>% group_by(year) %>% summarise(mean=mean(GOMchl))
plot(tt, type='b', main=paste("GOM annual Chl"))
tt.sp=select(occi.date, month, year, GOMchl) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(GOMchl))
plot(tt.sp, type='b',main=paste("GOM Jan-Jun Chl"))
tt.fl=select(occi.date, month, year, GOMchl) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(GOMchl))
plot(tt.fl, type='b', main=paste("GOM Jul-Dec Chl"))

tt=select(occi.date, month, year, MABchl) %>% group_by(year) %>% summarise(mean=mean(MABchl))
plot(tt, type='b', main=paste("MAB annual Chl"))
tt.sp=select(occi.date, month, year, MABchl) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(MABchl))
plot(tt.sp, type='b',main=paste("MAB Jan-Jun Chl"))
tt.fl=select(occi.date, month, year, MABchl) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(MABchl))
plot(tt.fl, type='b', main=paste("MAB Jul-Dec Chl"))

## now using log chl data, then taking exponent to standardize
tt2=select(occi.date, month, year, NESchllg) %>% group_by(year) %>% summarise(mean=mean(NESchllg)) %>% mutate(mean=10^mean)
plot(tt, type='b', main=paste("NES annual Chl"))
tt.sp=select(occi.date, month, year, NESchl) %>% filter(month<=6) %>% group_by(year) %>% summarise(mean=mean(NESchl))
plot(tt.sp, type='b',main=paste("NES Jan-Jun Chl"))
tt.fl=select(occi.date, month, year, NESchl) %>% filter(month>6) %>% group_by(year) %>% summarise(mean=mean(NESchl))
plot(tt.fl, type='b', main=paste("NES Jul-Dec Chl"))

