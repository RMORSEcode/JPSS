setwd('G:/1 RM')
setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM')
MM.chl=read.csv('MARMAPchlorophyll.csv', skip=74, stringsAsFactors = F)
MM.chl=read.csv('/home/ryan/Desktop/Z/MARMAPchlorophyll.csv', skip=74, stringsAsFactors = F)
colnames(MM.chl)=c('CRUISE', 'STA', 'YEAR', 'MON',	'DAY', 'HR',	'MIN',	'LATD',	'LOND',	'DEPTH',	'CHLA',	'PHAE')

MM.surf=MM.chl[which(MM.chl$DEPTH<15),] # limit depth to: surface - 15 m

# MM.chl2=MM.chl[complete.cases(MM.chl),]
MM.chl2=MM.surf[complete.cases(MM.surf),]

coordinates(MM.chl2)=~LOND+LATD #transform to Spatialpointsdataframe
pointsin=over(MM.chl2, neus.shp) #find which boxes samples belong to
BOX_ID=as.numeric(levels(pointsin$BOX_ID))[pointsin$BOX_ID]
DEPTH=as.numeric(levels(pointsin$DEPTH))[pointsin$DEPTH]
MM.boxbio=data.frame(MM.chl2, BOX_ID, DEPTH) #pointsin
### compute monthly means per box
MM.boxes.mon=aggregate(MM.boxbio$CHLA,list('box'=MM.boxbio$BOX_ID, 'M'= MM.boxbio$MON), mean)
MM.boxes.mon=reshape(MM.boxes.mon, idvar='box', timevar = 'M', direction = 'wide')
MM.boxes.mon=MM.boxes.mon[order(MM.boxes.mon['box']),]
### compute means per box by month and year
MM.boxes.ts=aggregate(MM.boxbio$CHLA,list('box'=MM.boxbio$BOX_ID, 'M'= MM.boxbio$MON, 'Y'=MM.boxbio$YEAR), median)

# climatology of box 1-22 (dynamic only boxes)
mm.fill=matrix(NA, ncol = 12, nrow=22)
mm.fill[,]=data.matrix(MM.boxes.mon[2:23,2:13])
t=gplots::heatmap.2(mm.fill, Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    col=my_palette2, tracecol = 'black')


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


# fill holes in mm.fill  to complete climatolgy
#boxes,month: 16-4 -> mean of10-4 and 11-4, 19-4;
#17-3,4; 18-3,4; 22-3,7
mm.fill2=mm.fill
mm.fill2[16,4]=mean(mm.fill[c(10:11,19),4])
mm.fill2[17,3]=mm.fill[10,3]
mm.fill2[17,4]=mean(c(mm.fill[10,3],mm.fill[17,5]))
mm.fill2[18,3]=.5*mm.fill[10,3]
mm.fill2[18,4]=0.75*mean(c(mm.fill[10,3],mm.fill[17,5]))
mm.fill2[22,3]=mm.fill[19,3]
mm.fill2[22,7]=mean(c(mm.fill[22,6],mm.fill[22,8]))
t=gplots::heatmap.2(mm.fill2, Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    col=my_palette2, tracecol = 'black')



## fill NAs in MM.b.clim with climatolgy mm.fill (boxes 1-30 -> 2:23 dynamic)
test=is.na(MM.b.clim$x.y)
MM.b.clim.final=MM.b.clim
for (i in 1:dim(MM.b.clim)[1]){
  if(test[i]==1){
    if(MM.b.clim$box[i] %in% c(1, seq(24,30,1))) {
      next
    }
    MM.b.clim.final$x.y[i]=mm.fill2[(MM.b.clim$box[i]-1),MM.b.clim$M[i]]
  }
}

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

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## fill in boundary boxes with mean climatology from surrounding boxes
mm.fill.full=matrix(NA, ncol=12, nrow=30)
mm.fill.full[2:23,]=mm.fill2
mm.fill.full[1,]=colMeans(mm.fill2[1:2,]) # box zero mean of 1 and 2
mm.fill.full[24:25,]=0
mm.fill.full[26,]=colMeans(mm.fill2[c(18,22),]) # box 25 mean of 18 and 22
mm.fill.full[27,]=colMeans(mm.fill2[21:22,]) # box 26 mean of 21 and 22
mm.fill.full[28,]=colMeans(mm.fill2[14:15,]) # box 27 mean of 14 and 15
mm.fill.full[29,]=colMeans(mm.fill2[c(5,6,9,14),]) # box 28 mean of 5,6,9,14
mm.fill.full[30,]=mm.fill2[3,] # box 29 == box 3 

marmap.clim=mm.fill.full
interpnoise(marmap.clim)


mm.ttxd=loadRData('/home/ryan/Desktop/Z/CZCS/ marmap.clim .Rdata')
mm.ttxdn=loadRData('/home/ryan/Desktop/Z/CZCS/ marmap.clim noise.Rdata')

#visualize chl per box as year long time series
plot(mm.ttxd[2,], type='l', ylim=c(0,20)) #max(mm.ttxd)))
for(i in 3:23){
  lines(mm.ttxd[i,])
}
barplot(table(round(as.numeric(mm.ttxd),1)), main='MARMAP Chl') #

# test without added noise
hirata.diat=apply(mm.ttxd, c(1,2), function(x) ((1.3272 + exp(-3.9828 * log10(x) + 0.1953))^-1))
hirata.micro=apply(mm.ttxd, c(1,2), function(x) ((0.9117 + exp(-2.733 * log10(x) + 0.4003))^-1))
hirata.pico=apply(mm.ttxd, c(1,2), function(x) (-(0.1529+(exp(1.0306*log10(x) -1.5576)))^-1 -1.8597*log10(x) + 2.9954))
# hirata.pico2=apply(mm.ttxd, c(1,2), function(x) ((-(0.1529+(exp((1.0306*log10(x)) -1.5576)))^-1) - (1.8597*log10(x)) + 2.9954)) # testing order of operations
# atl.DF=(hirata.micro-hirata.diat) * mm.ttxdn * 7 # (mg N m-3)
# atl.PL=hirata.diat * mm.ttxdn * 7 # (mg N m-3)
# atl.PS=(1-hirata.micro) * mm.ttxdn * 7 # (mg N m-3)
atl.DF=(hirata.micro-hirata.diat) * mm.ttxd * 7 # (mg N m-3)
atl.PL=hirata.diat * mm.ttxd * 7 # (mg N m-3)
atl.PS=(1-hirata.micro) * mm.ttxd * 7 # (mg N m-3)
barplot(table(round(as.numeric(atl.DF),2)),main='MARMAP DF') 
barplot(table(round(as.numeric(atl.PL),2)),main='MARMAP PL') 
barplot(table(round(as.numeric(atl.PS),2)),main='MARMAP PS')

# wide to long -> matrix to vector to place in netcdf file
# mm.ttxdnw=c(mm.ttxdn)
mm.ttxdw=c(mm.ttxd)

# load non-leap year netcdf file from Joe
jc1969=ncdf4::nc_open('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/Phyto_Forcing_1969.nc')
jcPL=ncvar_get(jc1969, varid = 'Diatom_N')
barplot(table(round(as.numeric(jcPL),1))) #
max(jcPL, na.rm = T)

## wide format as ncdf reads it....
jct1=jcPL[,,1] # z,b,t (day 1)
tznotna=which(!is.na(jct1), arr.ind=TRUE)
jcidx=tznotna[,1][order(tznotna[,2])] # row index for valid Chl (surface) to overwrite
## test assign 
test=jct1
# test=matrix(NA, ncol=30, nrow=5)
# test[jcidx,1:30]=atl.PL[,1]
test[tznotna]=atl.PL[,1]
# ncvar_put(jc1969, varid = 'Diatom_N', start=1,1,1)
## wide
test=matrix(NA, ncol=30, nrow=5)
test[tznotna]=atl.PL[,1]
test2=matrix(NA, ncol=30, nrow=5)
test2[tznotna]=atl.PL[,2]
test3=cbind(test, test2)

### long format 
# jct1=t(jcPL[,,1]) # transpose so boxes go in Y dim, depths in X
# tznotna=which(!is.na(jct1), arr.ind=TRUE)
# jcidx=tznotna[,2][order(tznotna[,1])] # column index for valid Chl (surface) to overwrite
# ## create new data to use in new file:
# PL=matrix(NA, ncol=5, nrow=30)
# PL[tznotna]=atl.PL[,1]
# for(i in 2:dim(atl.PL)[2]){
#   test2=matrix(NA, ncol=5, nrow=30)
#   test2[tznotna]=atl.PL[,i]
#   PL=rbind(PL, test2)
# }
# DF=matrix(NA, ncol=5, nrow=30)
# DF[tznotna]=atl.DF[,1]
# for(i in 2:dim(atl.DF)[2]){
#   test2=matrix(NA, ncol=5, nrow=30)
#   test2[tznotna]=atl.DF[,i]
#   DF=rbind(DF, test2)
# }
# PS=matrix(NA, ncol=5, nrow=30)
# PS[tznotna]=atl.PS[,1]
# for(i in 2:dim(atl.PS)[2]){
#   test2=matrix(NA, ncol=5, nrow=30)
#   test2[tznotna]=atl.PS[,i]
#   PS=rbind(PS, test2)
# }
## create 3-D array want (t,b,z) -> have (b,z,t); Joe: (z,b,t)
# PL=matrix(NA, ncol=5, nrow=30)
# PL[tznotna]=atl.PL[,1]
# for(i in 2:dim(atl.PL)[2]){
#   test2=matrix(NA, ncol=5, nrow=30)
#   test2[tznotna]=atl.PL[,i]
#   PL=abind::abind(PL, test2, along=3)
# }

### try wide instead of long.... matches what Joe has
## create 3-D array want (t,b,z) -> long format gives (b,z,t); Joe: (z,b,t)
PL=matrix(NA, ncol=30, nrow=5)
PL[tznotna]=atl.PL[,1]
for(i in 2:dim(atl.PL)[2]){
  test2=matrix(NA, ncol=30, nrow=5)
  test2[tznotna]=atl.PL[,i]
  PL=abind::abind(PL, test2, along=3)
}
DF=matrix(NA, ncol=30, nrow=5)
DF[tznotna]=atl.DF[,1]
for(i in 2:dim(atl.DF)[2]){
  test2=matrix(NA, ncol=30, nrow=5)
  test2[tznotna]=atl.DF[,i]
  DF=abind::abind(DF, test2, along=3)
}
PS=matrix(NA, ncol=30, nrow=5)
PS[tznotna]=atl.PS[,1]
for(i in 2:dim(atl.PS)[2]){
  test2=matrix(NA, ncol=30, nrow=5)
  test2[tznotna]=atl.PS[,i]
  PS=abind::abind(PS, test2, along=3)
}

### create netcdf forcing file for MARMAP climatology
wd='/home/ryan/Desktop/Z/CZCS/'
library(rbgm) ## read BGM
dz_box=read.csv('/home/ryan/AtlRuns/dz.csv', header=T) # depths for each layer by box
bfile="/home/ryan/Git/atneus_RM/neus_tmerc_RM.bgm" ### new file from Bec, modified faces
bgm <- bgmfile(bfile)
boxes <- boxSpatial(bgm)
box.boxes=bgm$boxes$.bx0 ### NOTE added 1 because index cannot be 0, must remove 1 later w/ ncks

### Define dimensions 
nboxes=30 # length(unique(box_props$.bx0)) #30
ntimes=365 #length(unique(face_props2$band_level))
nlevel=5 #length(unique(face_props2$atlantis_level))
dt=86400 # seconds in one day
# t_stop=max(unique(face_props2$band_level))+ynum ### MAKE SURE TO SELECT ABOVE ynum CORRECTLY
# t_tot=seq(t_start*dt,t_stop*dt,dt)
# jc1964=ncdf4::nc_open('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/Phyto_Forcing_1964.nc')
# jctime=ncvar_get(jc1964, varid = 't')
# t_tot=jctime[1:365]
# filename="/home/ryan/Desktop/Z/CZCS/Atlantis_force/MARMAP_climatology_phyto_forcing.nc"

## get list of years and loop over to extract time to write individual year files with climatology
# yearlist=seq(from=1964, to=1997, by=1)
jclist=list.files('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/', pattern=glob2rx("*Phyto_Forcing_*.nc*"))
for (i in (1:35)){
  jcyear=ncdf4::nc_open(paste('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/', jclist[i], sep=''))
  jctime=ncvar_get(jcyear, varid = 't')
  t_tot=jctime[1:365]
  filename=paste("/home/ryan/Desktop/Z/CZCS/Atlantis_force/MARMAP_climatology_",jclist[i], sep='')
  #define dimensions
  timedim=ncdim_def("t", "", 1:length(t_tot), unlim=T, create_dimvar = F) #as.double(t_tot)
  leveldim=ncdim_def("z", "", 1:nlevel, create_dimvar = F)
  boxesdim=ncdim_def("b", "", 1:nboxes, create_dimvar = F)
  #create variables
  var.time=ncvar_def("t","seconds since 1964-01-01 00:00:00 +10",timedim,prec="double")
  var.PL=ncvar_def("Diatom_N","mg N m-3",list(leveldim, boxesdim, timedim),-999,longname="Diatom_N",prec="double")
  var.PLS=ncvar_def("Diatom_S","mg Si m-3",list(leveldim, boxesdim, timedim),-999,longname="Diatom_S",prec="double")
  var.DF=ncvar_def("DinoFlag_N","mg N m-3",list(leveldim, boxesdim, timedim),-999,longname="DinoFlag_N",prec="double")
  var.PS=ncvar_def("PicoPhytopl_N","mg N m-3",list(leveldim, boxesdim, timedim),-999,longname="PicoPhytopl_N",prec="double")
  nc_varfile=nc_create(filename,list(var.time, var.PL, var.PLS, var.DF, var.PS)) #var.box, var.lev
  #assign global attributes to file
  ncatt_put(nc_varfile,0,"title","MARMAP chlorophyll forcing file NEUS")
  ncatt_put(nc_varfile,0,"geometry","neus_tmerc_RM2.bgm")
  #assign attributes to variables
  ncatt_put(nc_varfile, "t", "dt", prec="double", 86400)
  ncatt_put(nc_varfile,'Diatom_N', '_FillValue', prec="double", -999)
  ncatt_put(nc_varfile,'Diatom_N', 'missing_value', prec="double", -999)
  ncatt_put(nc_varfile,'Diatom_N', 'valid_min', prec="double", -999)
  ncatt_put(nc_varfile,'Diatom_N', 'valid_max', prec="double", 99999)
  ncatt_put(nc_varfile,'Diatom_N', 'units', 'mg N m-3')
  ncatt_put(nc_varfile,'Diatom_N', 'long_name', 'Diatom_N')
  ncatt_put(nc_varfile,'Diatom_S', '_FillValue', prec="double", -999)
  ncatt_put(nc_varfile,'Diatom_S', 'missing_value', prec="double", -999)
  ncatt_put(nc_varfile,'Diatom_S', 'valid_min', prec="double", -999)
  ncatt_put(nc_varfile,'Diatom_S', 'valid_max', prec="double", 99999)
  ncatt_put(nc_varfile,'Diatom_S', 'units', 'mg Si m-3')
  ncatt_put(nc_varfile,'Diatom_S', 'long_name', 'Diatom_S')
  ncatt_put(nc_varfile,'DinoFlag_N', '_FillValue', prec="double", -999)
  ncatt_put(nc_varfile,'DinoFlag_N', 'missing_value', prec="double", -999)
  ncatt_put(nc_varfile,'DinoFlag_N', 'valid_min', prec="double", -999)
  ncatt_put(nc_varfile,'DinoFlag_N', 'valid_max', prec="double", 99999)
  ncatt_put(nc_varfile,'DinoFlag_N', 'units', 'mg N m-3')
  ncatt_put(nc_varfile,'DinoFlag_N', 'long_name', 'DinoFlag_N')
  ncatt_put(nc_varfile,'PicoPhytopl_N', '_FillValue', prec="double", -999)
  ncatt_put(nc_varfile,'PicoPhytopl_N', 'missing_value', prec="double", -999)
  ncatt_put(nc_varfile,'PicoPhytopl_N', 'valid_min', prec="double", -999)
  ncatt_put(nc_varfile,'PicoPhytopl_N', 'valid_max', prec="double", 99999)
  ncatt_put(nc_varfile,'PicoPhytopl_N', 'units', 'mg N m-3')
  ncatt_put(nc_varfile,'PicoPhytopl_N', 'long_name', 'PicoPhytopl_N')
  #assign variables to file
  ncvar_put(nc_varfile, var.time, t_tot, count=ntimes)
  ncvar_put(nc_varfile,var.PL, PL, count=c(nlevel, nboxes, ntimes))
  ncvar_put(nc_varfile,var.PLS, PL, count=c(nlevel, nboxes, ntimes))
  ncvar_put(nc_varfile,var.DF, DF, count=c(nlevel, nboxes, ntimes))
  ncvar_put(nc_varfile,var.PS, PS, count=c(nlevel, nboxes, ntimes))
  nc_close(nc_varfile)
}
##### end force file creation





### actual box numbers, 0 and other boundary boxes omitted ### dynamic 1-22 ONLY
test=MM.b.clim.final
test$D=15
dt=paste(test$Y,test$M, test$D, sep='-')
test$dt=as.Date(dt)
test2=test[,c('box', 'x.y', 'dt')]
MM.box.mon.ts=test2 %>% tidyr::spread(dt, x.y, fill = NA, convert = FALSE)
t=gplots::heatmap.2(as.matrix(MM.box.mon.ts[2:23,2:133]), Rowv=F, Colv=F, dendrogram = 'none', na.color='gray',
                    col=my_palette2, tracecol = 'black')



### create regimes:
## plot 1 year from 1 box, repeat for long term base
test=cz.ttxdn[13,]
t2=rep(test, 20)
plot(t2, type='l')


### FOR CZCS convert from Chl (mg/m3) to mg N for Atlantis forcing
# PS=0.77[1-exp(-0.8/0.13)*CHL] # Lamont et al 2019
# PL=([1.3272 + exp(-3.9828 * log10(CHL) + 0.1953)]^-1)*100 * CHL # Hirata et al. 2011, table 2 modified
# Cm=CHL-PS # microplankton contribution to Chl concentration
# DF=Cm-PL

## all need to be converted to Atlantis units of N
# CHLN=CHL*7
#nano+pico fraction
# test.ps=apply(ttxdn, c(1,2), function(x) 0.77*(1-(exp(-0.8/0.13)*x))) 
# pico only fraction
# test.pico=apply(ttxdn, c(1,2), function(x) 0.13*(1-(exp(-0.8/0.13)*x)))
# nano fraction
# test.nano=test.ps-test.pico
# microplankton CHL
# test.micro=ttxdn-(test.ps*ttxdn) # CHL-CHL_pico+nano
# microplankton fraction
# test.micro=1-test.ps
# test.DF=test.micro-(test.PL*ttxdn)
# test.PL=apply(test.micro, c(1,2), function(x) ((1.3272 + exp(-3.9828 * log10(x) + 0.1953))^-1))

## load CZCS climatology interpolated daily product from 'testing_hermes_occi.r'
cz.ttxdn=loadRData('/home/ryan/Desktop/Z/CZCS/ czcsclim noise.Rdata') # ttxdn
cz.ttxd=loadRData('/home/ryan/Desktop/Z/CZCS/ czcsclim .Rdata') # ttxd

barplot(table(round(as.numeric(cz.ttxd),1)), main="CZCS Chl") #
barplot(table(round(as.numeric(cz.ttxdn),1)), main="CZCS Chl +N") # different structure with noise, more low vals
quantile( cz.ttxd, c(.1, .9 ) )
quantile( cz.ttxdn, c(.1, .9 ))

funt <- function(x){
  quantiles <- quantile( x, c(.05, .90 ) )
  x[ x < quantiles[1] ] <- quantiles[1]
  x[ x > quantiles[2] ] <- quantiles[2]
  x
}
# td=funt(cz.ttxd)
# tdn=funt(cz.ttxdn)
# barplot(table(round(as.numeric(td),1))) 
# barplot(table(round(as.numeric(tdn),1))) 
jc1964=ncdf4::nc_open('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/Phyto_Forcing_1964.nc')
jctime=ncvar_get(jc1964, varid = 'Time')
rmtime=jctime[1:365]






jc1968=ncdf4::nc_open('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/Phyto_Forcing_1968.nc')
## leap years present 1964, 1968
jc1969=ncdf4::nc_open('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/Phyto_Forcing_1969.nc')
jcPL=ncvar_get(jc1969, varid = 'Diatom_N')
jctime=ncvar_get(jc1964, varid = 'Time')
barplot(table(round(as.numeric(jcPL),1))) #
max(jcPL, na.rm = T)

jc2002=ncdf4::nc_open('/home/ryan/Git/neus_atlantis/neus-atlantis/currentVersion/tsfiles/Annual_Files/Phyto_Forcing_2002.nc')
jcPLx=(ncvar_get(jc2002, varid = 'Diatom_N'))
barplot(table(round(as.numeric(jcPLx),1)), main='JCdiatomN') #
max(jcPLx, na.rm = T)
jcPSx=ncvar_get(jc2002, varid = 'PicoPhytopl_N')
barplot(table(round(as.numeric(jcPLx),1)), main='JCpicoN') #

# wide to long -> matrix to vector to place in netcdf file
cz.ttxdnw=c(cz.ttxdn)
cz.ttxdw=c(cz.ttxd)

# NOW just using Hirata et al 2011 global values: fraction * CHL * X_CHLN
hirata.diat=apply(cz.ttxdn, c(1,2), function(x) ((1.3272 + exp(-3.9828 * log10(x) + 0.1953))^-1))
hirata.micro=apply(cz.ttxdn, c(1,2), function(x) ((0.9117 + exp(-2.733 * log10(x) + 0.4003))^-1))
atl.DF=(hirata.micro-hirata.diat) * cz.ttxdn * 7 # (mg N m-3)
atl.PL=hirata.diat * cz.ttxdn * 7 # (mg N m-3)
atl.PS=(1-hirata.micro) * cz.ttxdn * 7 # (mg N m-3)
barplot(table(round(as.numeric(atl.DF),2)),main='CZCS +n DF') 
barplot(table(round(as.numeric(atl.PL),2)),main='CZCS +n PL') 
barplot(table(round(as.numeric(atl.PS),2)),main='CZCS +n PS')

# test without added noise
hirata.diat=apply(cz.ttxd, c(1,2), function(x) ((1.3272 + exp(-3.9828 * log10(x) + 0.1953))^-1))
hirata.micro=apply(cz.ttxd, c(1,2), function(x) ((0.9117 + exp(-2.733 * log10(x) + 0.4003))^-1))
hirata.pico=apply(cz.ttxd, c(1,2), function(x) (-(0.1529+(exp(1.0306*log10(x) -1.5576)))^-1 -1.8597*log10(x) + 2.9954))
# hirata.pico2=apply(cz.ttxd, c(1,2), function(x) ((-(0.1529+(exp((1.0306*log10(x)) -1.5576)))^-1) - (1.8597*log10(x)) + 2.9954)) # testing order of operations
atl.DF=(hirata.micro-hirata.diat) * cz.ttxdn * 7 # (mg N m-3)
atl.PL=hirata.diat * cz.ttxdn * 7 # (mg N m-3)
atl.PS=(1-hirata.micro) * cz.ttxdn * 7 # (mg N m-3)
atl.DF=(hirata.micro-hirata.diat) * cz.ttxd * 7 # (mg N m-3)
atl.PL=hirata.diat * cz.ttxd * 7 # (mg N m-3)
atl.PS=(1-hirata.micro) * cz.ttxd * 7 # (mg N m-3)
barplot(table(round(as.numeric(atl.DF),2)),main='CZCS DF') 
barplot(table(round(as.numeric(atl.PL),2)),main='CZCS PL') 
barplot(table(round(as.numeric(atl.PS),2)),main='CZCS PS')
abline(v=0)
