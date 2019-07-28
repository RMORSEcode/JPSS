# 
# url<- "ftp://ftp.hermes.acri.fr/345576024"
# filenames <- getURL(url, userpwd="ftp_hermes:hermes%", ftp.use.epsv = FALSE, dirlistonly = TRUE) #reading filenames from ftp-server
# destnames <- filenames <-  strsplit(filenames, "\r*\n")[[1]] # destfiles = origin file names
# con <-  getCurlHandle( ftp.use.epsv = FALSE, userpwd="ftp_hermes:hermes%")
# setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES_daily')
# mapply(function(x,y) writeBin(getBinaryURL(x, curl = con, dirlistonly = FALSE), y), x = filenames, y = paste("C:\\temp\\",destnames, sep = "")) #writing all zipped files in one directory
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
library(maps)
library(mapdata)
library(rgeos)
library(ncdf4)
library(abind)
library(RColorBrewer)

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


### Function to plot a raster for the NES area and set scale
plotChlRaster=function(data, i, maxV, limit=F){
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
       las=1, legend=F, main=av.dates[[i]])
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
nc1=nc_open('CCI_ALL-v4.0-8DAY.nc')
lon=ncvar_get(nc1, 'lon')
lat=ncvar_get(nc1, 'lat')
chl=ncvar_get(nc1, 'chlor_a')
time=ncvar_get(nc1, 'time') #days since Jan 1, 1970
dim(chl)
colnames(chl)=lat
rownames(chl)=lon
nc_close(nc1)


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

### create raster stacks of data 
chl.av=nc2raster(files.chl1.av[,1], 'CHL1_mean')
chl.gsm=nc2raster(files.chl1.gsm[,1], 'CHL1_mean')
chl.oc5=nc2raster(files.oc5[,1], 'CHL-OC5_mean')

nc1=nc_open(files.oc5[1,1])

i=20
plot(chl.av[[i]], main=test1[i])

# testing...
plotChlRaster(chl.av, 5, 5, limit=T)
plotChlRaster(chl.av, 5, 20, limit=F)

# create time series by box for ecomon strata
m.av=extract_calc(chl.av[[1:262]], ecomon.strata)
m.gsm=extract_calc(chl.gsm[[1:262]], ecomon.strata)
#CBay
plot(m.av[6,], type='l')
lines(m.gsm[6,], col='red')
#GBK
plot(m.av[30,], type='l')
lines(m.gsm[30,], col='red')
#NY
plot(m.av[17,], type='l')
lines(m.gsm[17,], col='red')

plot((m.av[1,]-m.gsm[1,]), type='l')
abline(h=0, lty=2)






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
