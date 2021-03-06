
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(foreach)
library(doParallel)
library(signal)
library(robustbase)
library(MODIS)
library(plyr)

library(gimms)
library(ncdf4)
library(rgdal)
library(rasterVis)
library(rgeos)
library(FRutils)
library(signal)
library(scales)

library(MODIStsp)
library(velox)
library(mgcv)
library(rgeos)
library(tmap)
library(data.table)
library(scales)
library(quantreg)
library(FlexParamCurve)
library(signal)
library(zoo)
library(rasterVis)
library(FRutils)
library(phenex)
library(doSNOW)
library(Deriv)

library(vegan)

library(mapview)
library(mapedit)
library(magrittr)

load("C:/Users/rouf1703/Desktop/ram2000_2018.RData")

### This script is for extracting ndvi/evi metrics from RasterStack objects using a set of regions defined by polygons

# RasterStacks may be of the MODIS or GIMMS type
# A raster of julian dates is assumed to accompany each value raster

# y = a year as an integer
# j = the julian day to convert
# j2 = the julian day that marks the beginning of a block


jul2nd<-function(y,j,j2=NULL){
  inc<-0
  if(!is.null(j2)){
    inc<-ifelse(j<j2,1,0) 
  }
  as.integer(as.Date(paste0(y+inc,"-01-01")))+j-1
}

findminmax<-function(x,n=1,beg="06-01",end="11-01",max=TRUE){
  #if(!max)browser()
  stopifnot(!is.null(names(x)))
  d<-substr(names(x),6,10)
  bloc<-d>=beg & d<=end
  run<-rle(bloc)
  # first version
  #l<-Map(":",c(1,head(cumsum(run[[1]]),-1))[run[[2]]],cumsum(run[[1]])[run[[2]]])
  # second version
  l<-Map(":",c(1,head(cumsum(run[[1]])+1,-1)),cumsum(run[[1]]))
  l<-l[run[[2]]]
  #
  r<-rank(ifelse(max,-1,1)*x)
  res<-lapply(l,function(i){
    val<-sort(r[i])[1:n]
    index<-match(val,r)
    index   
  })
  res
}

rangeinc<-function(x,inc=c(0.2,0.1)){
  r<-range(x)
  d<-r[2]-r[1]
  c(r[1]-d*inc[1],r[2]+d*inc[2])
}



#### test logistic

### check package sicegar and associated peerj paper for double sigmoïd curve fitting

### simple logictic
fitLog<-function(x,mmdate=c("12-01","09-15"),plot=FALSE){
  years<-as.integer(unique(substr(names(x),1,4)))
  years<-years[years<=2016]
  l<-lapply(years,function(i){
    if(mmdate[1]>mmdate[2]){
      paste(c(i-1,i),mmdate,sep="-") 
    }else{
      paste(c(i-1,i),mmdate,sep="-")
    }
  })
  peak<-lapply(l,function(i){
    sx<-x[which(names(x)>=i[1] & names(x)<=i[2])]
    d<-data.frame(y=sx,x=as.integer(as.Date(names(sx))))
    min_xmid_up<-as.integer(as.Date(i[1])+120)
    max_xmid_up<-as.integer(as.Date(i[2])-75)
    lo<-list(Asym_up=0,xmid_up=min_xmid_up,scal_up=15,c=0.0)
    up<-list(Asym_up=1,xmid_up=max_xmid_up,scal_up=40,c=0.8)
    m1<-tryCatch(nls(y~Asym_up/(1+exp((xmid_up-x)/scal_up))+c,data=d,start=list(Asym_up=0.5,xmid_up=min_xmid_up+(max_xmid_up-min_xmid_up)/2,scal_up=20,c=0.2),control=list(minFactor=1e-12,maxiter=500),lower=lo,upper=up,algorithm="port"),error=function(j){TRUE})
    if(!isTRUE(m1)){
      se<-seq(min(d$x),max(d$x),by=1) 
      if(plot){  
        #plot(as.Date(d$x),d$y)
        lines(as.Date(se),predict(m1,data.frame(x=se)),col=alpha("grey40",0.5),lwd=4)
      }
      coef(m1)
    }else{
      NA  
    }
  })
  peak
}

### double logistic with Beck's et al (2006) parametrisation
DLog<-function(x,wNDVI,mNDVI,S,A,mA,mS){
  wNDVI+(mNDVI-wNDVI)*((1/(1+exp(-mS*(x-S))))+(1/(1+exp(mA*(x-A))))-1)  
}

DLog1<-Deriv(DLog,"x")

### double logistic
fitDLog<-function(x,mmdate=c("12-01","03-15"),plot=FALSE,type="ndvi",...){
  years<-as.integer(unique(substr(names(x),1,4)))
  #years<-years[years<=2018]
  l<-lapply(years,function(i){
    if(mmdate[1]>mmdate[2]){
      paste(c(i-1,i+1),mmdate,sep="-") 
    }else{
      paste(c(i-1,i+1),mmdate,sep="-")
    }
  })
  peak<-lapply(l,function(i){
    sx<-x[which(names(x)>=i[1] & names(x)<=i[2])]
    d<-data.frame(y=sx,x=as.integer(as.Date(names(sx))))
    
    yy<-as.integer(substr(i[1],1,4))+1
    
    min_S<-as.integer(as.Date(paste0(yy,"-04-01")))
    max_S<-as.integer(as.Date(paste0(yy,"-07-01")))
    
    min_A<-as.integer(as.Date(paste0(yy,"-08-15")))
    max_A<-as.integer(as.Date(paste0(yy,"-12-01")))
    
    ### Beck's et al. (2006) parametrisation 
    
    lo<-list(wNDVI=min(x),S=min_S,mS=0.03,mNDVI=median(x),A=min_A,mA=0.03)
    up<-list(wNDVI=median(x),S=max_S,mS=0.1,mNDVI=max(x),A=max_A,mA=0.1)
    
    if(type=="snow"){

      lo<-list(wNDVI=0,S=min_S,mS=0.03,mNDVI=0.5,A=min_A,mA=0.03)
      up<-list(wNDVI=0.5,S=max_S,mS=0.5,mNDVI=1,A=max_A,mA=0.5)
    }

    start<-mapply(function(i,j){i+j},lo,mapply(function(x,y){(x-y)/2},up,lo,SIMPLIFY=FALSE),SIMPLIFY=FALSE)
    
    if(nrow(d)<10){
      NA
    }else{
      m1<-tryCatch(nls(y~DLog(x,wNDVI,mNDVI,S,A,mA,mS),data=d,start=start,control=list(minFactor=1e-12,maxiter=500),lower=lo,upper=up,algorithm="port"),error=function(j){TRUE})
      
      if(!isTRUE(m1)){
        se<-seq(min(d$x),max(d$x),by=1) 
        if(plot){  
          #plot(as.Date(d$x),d$y)
          lines(as.Date(se),predict(m1,data.frame(x=se)),lwd=4,...)
        }
        coef(m1)
      }else{
        NA  
      }
    }
  })
  names(peak)<-years
  peak
}

# x is a file name ending with "_2002_342.tif"
getDates<-function(x){
    x<-gsub(".tif","",x)
    n<-nchar(x)
    s<-substr(x,n-7,n)
    s<-strsplit(s,"_")
    as.Date(paste0(sapply(s,"[",1),"-01-01"))+as.integer(sapply(s,"[",2))-1
}

#getDate("whatata_2001_324.tif")

#############################################################
##### Get regions of interest ###############################

# A list of SpatialPolygonsDataFrame for which a metric is wanted

### Ram mountain

#x <- mapview() %>% editMap()
#plot(x$finished)
#ram<-as(x$finished,"Spatial")
#writeOGR(ram,dsn="C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram",driver="ESRI Shapefile")




### elevation data comes from tile 083b05: 
# http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec/50k_dem/
e<-raster("S:/NDVI/MODIS/083b05/083b05_0201_deme.dem")
w<-raster("S:/NDVI/MODIS/083b05/083b05_0201_demw.dem")
alt<-merge(e,w)

###  build a tile to crop the MODIS rasters when dowloading data with MODIStsp
#tileram<-as(extent(alt),"SpatialPolygons")
#proj4string(tileram)<-proj4string(alt)
#writeOGR(SpatialPolygonsDataFrame(tileram,data=data.frame(id=1)),dsn="S:/NDVI/MODIS/083b05",layer="tileram",driver="ESRI Shapefile",overwrite_layer=FALSE)
#plot(alt)
#plot(tileram,add=TRUE)


a<-alt
a<-aggregate(alt,2)
av<-1600
a[a<av]<-0
a[a>=av]<-1




ram<-rasterToPolygons(a,fun=function(x){x>=1},dissolve=TRUE)
ram<-disaggregate(ram)
o<-over(ram,SpatialPoints(cbind(c(-115.8549,-115.7963),c(52.38193,52.35439)),proj4string=CRS(proj4string(ram))))
ram<-ram[!is.na(o),]
plot(alt)
plot(ram,add=TRUE)
#plot(ram,add=TRUE)

#x<-list.files("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/test/Gross_PP_8Days_500m_v6/GPP",full.names=TRUE)
#v<-velox(sapply(x,velox))
#r<-v$as.RasterStack()

#plot(r[[1]])
#plot(spTransform(ram,CRS(proj4string(r))),add=TRUE)
#plot(spTransform(patch,CRS(proj4string(r))),add=TRUE)


#ram<-gBuffer(ram,width=-0.015)
#pol<-bbox2pol(c(-72.26,-72.15,45.29,45.40))
#pol<-SpatialPolygonsDataFrame(pol,data.frame(id=1),match.ID=FALSE)

### the ram area could possibly be shrinked

#############################################################
##### Get MODIS raster ######################################

### check which and how to use QA and pixel reliability

### read multiple .tif instead of RData
#path<-"S:/NDVI/MODIS/Ram/VI_16Days_250m_v6/NDVI/"
#x<-list.files(path,pattern=".tif",full.names=TRUE)
#x<-x[-grep(".aux.xml",x)]
#x<-x[order(substr(x,nchar(x)-12,nchar(x)))]
#r2<-stack(x)


### Stack single image files

#l<-list.files("S:/NDVI/MODIS/VI_16Days_500m_v6/NDVI",full.names=TRUE)
#modis<-stack(l)
#l<-list.files("S:/NDVI/MODIS/VI_16Days_500m_v6/DOY",full.names=TRUE)
#modis_jul<-stack(l)

#lr<-lapply(l[1:10],function(i){velox(i,extent=extent(pol))})

#rasterOptions(chunksize = 1e+10)
#rasterOptions(maxmemory = 1e+12)


#path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/MODIS/VI_16Days_250m_v6/Time_Series/RData/"

### loading RData session
load("S:/NDVI/MODIS/VI_16Days_500m_v6/Time_Series/RData/Mixed/NDVI/MOD13A1_MYD13A1_NDVI_49_2000_1_2018_RData.RData")
modis<-raster_ts
r<-raster_ts
#load("S:/NDVI/MODIS/VI_16Days_500m_v6/Time_Series/RData/Mixed/DOY/MOD13A1_MYD13A1_DOY_49_2000_1_2018_RData.RData")
#modis_jul<-raster_ts
#load("S:/NDVI/MODIS/VI_16Days_500m_v6/Time_Series/RData/Mixed/Rely/MOD13A1_MYD13A1_Rely_49_2000_1_2018_RData.RData")
#modis_rely<-raster_ts
#load("S:/NDVI/MODIS/VI_16Days_500m_v6/Time_Series/RData/Mixed/VI_QA/MOD13A1_MYD13A1_VI_QA_49_2000_1_2018_RData.RData")
#modis_QA<-raster_ts
#load("S:/NDVI/MODIS/LAI_8Days_500m_v6/Time_Series/RData/Mixed/Lai/MOD15A2H_MYD15A2H_Lai_49_2000_1_2018_RData.RData")
#modis_lai<-raster_ts
#load("S:/NDVI/MODIS/Snow_Cov_8-Day_500m_v6/Time_Series/RData/Mixed/MAX_SNW/MOD10_A2_MYD10_A2_MAX_SNW_55_2000_1_2018_RData.RData")
##modis_snow<-raster_ts
#load("S:/NDVI/MODIS/Gross_PP_8Days_500m_v6/Time_Series/RData/Mixed/GPP/MOD17A2H_MYD17A2H_GPP_1_2001_1_2018_RData.RData")
#modis_gpp<-raster_ts
#rm(raster_ts)
#l<-strsplit(names(modis),"_")
#modis_doy<-as.Date(paste0(sapply(l,"[",3),"-01-01"))+as.integer(sapply(l,"[",4))-1


pol<-spTransform(ram,CRS(proj4string(modis)))

tmap_mode("view")
tm_shape(pol) +
  tm_polygons(col="red",alpha=0.3) +
  tm_borders(col="red",lwd=5) +
  tm_shape(r[[1]]) +
  tm_raster() +
  tm_layout(basemaps=c("Esri.WorldImagery"))

#plot(test)
#test<-rasterToPolygons(r[[1]])
#ocells<-over(test,spTransform(ram,proj4string(test)))
#plot(test[!is.na(ocells[,1]),],border="red",lwd=2,add=TRUE)


pol<-spTransform(ram,CRS(proj4string(modis)))

#registerDoParallel(8)
#no_cores <- detectCores(logical=FALSE) 
#registerDoParallel(no_cores)
#getDoParWorkers()

lpaths<-c(
"S:/NDVI/MODIS/VI_16Days_500m_v6/DOY",                                          
"S:/NDVI/MODIS/VI_16Days_500m_v6/EVI",                                          
"S:/NDVI/MODIS/VI_16Days_500m_v6/NDVI",                                         
"S:/NDVI/MODIS/VI_16Days_500m_v6/QA_qual",                                      
"S:/NDVI/MODIS/VI_16Days_500m_v6/QA_usef",                                      
"S:/NDVI/MODIS/VI_16Days_500m_v6/Rely", 
"S:/NDVI/MODIS/VI_16Days_500m_v6/VI_QA",
"S:/NDVI/MODIS/Snow_Cov_8-Day_500m_v6/MAX_SNW",                                 
"S:/NDVI/MODIS/Snow_Cov_8-Day_500m_v6/SC_8DAY_bitfield", 
"S:/NDVI/MODIS/LAI_8Days_500m_v6/Fpar",                                         
"S:/NDVI/MODIS/LAI_8Days_500m_v6/Lai",
"S:/NDVI/MODIS/Gross_PP_8Days_500m_v6/GPP",                                     
"S:/NDVI/MODIS/Gross_PP_8Days_500m_v6/PsnNet")

cl <- makeCluster(8)
registerDoSNOW(cl)

lr<-vector(mode="list",length=length(lpaths))
names(lr)<-sapply(strsplit(lpaths,"/"),tail,1)

for(i in 1:length(lpaths)){
  l<-list.files(lpaths[i],full.names=TRUE,pattern=".tif")
  jj<-substr(l,nchar(l)-11,nchar(l)-4)
  l<-l[jj>="2000_001" & jj<="2018_001"] # for a tiny subset######################
  g<-grep("MCD15",l)
  if(any(g)){  
    l<-l[-g]
  }
  lid<-substr(l,nchar(l)-11,nchar(l)-4)
  l<-l[order(lid)]
  lid<-lid[order(lid)]
  cat("\n",paste("",i," / ",length(lpaths)),"\n")
  pb<-txtProgressBar(max=length(l),style=3)
  progress<-function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  v<-foreach(j=1:length(l),.packages=c("velox"),.verbose=FALSE,.options.snow=opts) %dopar% {velox(l[j])}
  v<-velox(v)
  m<-v$extract(pol)
  m<-lapply(m,function(k){
    dimnames(k)[[2]]<-as.character(getDates(l))
    k
  })
  lr[[i]]<-m
}

close(pb)
stopCluster(cl) 




m<-v$extract(pol)
vm<-v$getCoordinates()
vm<-vm[order(vm[,1],-vm[,2]),]
cents<-SpatialPoints(vm,proj4string=CRS(v$crs))
par(mar=c(0,0,0,0))
plot(modis[[1]])
plot(cents,add=TRUE,pch=1,cex=0.6)

ocells<-over(cents,pol)



#modis<-crop(modis,pol)
#snow<-crop(snow,pol)

### keep only data from when Aqua is also present (not sure if this is relevant)
#m<-match("MYD",substr(names(modis),1,3))
#keep<-m:length(names(modis))
#keep<-1:length(names(modis))
#keep<-grep("MOD",names(r))
### subset modis to only get data from when Aqua was also used
#modis<-subset(modis,keep)
#modis_jul<-subset(modis_jul,keep)
#modis_doy<-modis_doy[keep]
#modis_rely<-modis_rely[[keep]]

# Alpha and Terra data are offset by a certain amount in their blocks, but the specific date on which the ndvi is measured for a pixel may be earlier in a later block. Thus, although blocks are ordered, it does not mean that the specific dates are ordered. Ideally, order dates for each pixel one the observations are out of the complete matrix

#for(i in seq_along(dim(modis_jul))){
#  year<-unlist(strsplit(names(modis)[i],"_"))[[3]]
#  dateseq<-as.integer(seq.Date(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),by=1))
#  jul2n<-function(x){
#    dateseq[x]  
#  }
#  res<-calc(modis_jul[,,i],fun=jul2n)
#  modis_jul<-setValues(modis_jul,layer=i,  
#}


#day<-lapply(1980:2017,function(i){
#  seq.Date(as.Date(paste0(i,"-01-01")),as.Date(paste0(i,"-12-31")),by=1)
#})



#rc<-SpatialPointsDataFrame(coordinates(modis),proj4string=CRS(proj4string(modis)),data.frame(id=seq_len(dim(modis)[1]*dim(modis)[2])))

### Rely

#load("C:/Users/rouf1703/Desktop/ram2000_2018.RData")

#pix<-90
#idx<-"EVI"
#plot(lr[[idx]][[1]][pix,])
#w<-which(!lr$Rely[[1]]%in%c(1,2,NA)) # values in here are supposed to be -1,0,1,2,3 not sure why they differ from the user guide
#lr[[idx]][[1]][w]<-NA
#points(lr[[idx]][[1]][pix,],pch=16,col="darkgreen")



#############################################################
##### Get GIMMS raster ######################################

path<-"S:/NDVI/GIMMS/"
#path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/GIMMS/"

x<-list.files(path)
ts<-monthlyIndices(x, version = 1, timestamp = T)
gimms_doy<-as.Date(as.character(ts+round(c(diff(ts),17)/2,0)))
gimms<-rasterizeGimms(x=paste0(path,x),ext=ram,cores=6L) # clipping
names(gimms)<-as.character(gimms_doy)
gimms_jul<-stack(setValues(gimms,as.integer(format(rep(gimms_doy,each=ncell(gimms)),"%j"))))







#############################################################
##### Extract raster values for each region #################

### the following part needs to be run twice, each for gimms and modis and adjust silenced parts

#r<-gimms # determine series to use
#rd<-gimms_jul
#bdoy<-gimms_doy
#divide<-1 # divide values by 1

ndvi<-lr$NDVI # determine series to use
evi<-lr$EVI # determine series to use
doy<-lr$DOY
rely<-lr$Rely
lai<-lr$Lai
snow<-lr$MAX_SNW
gpp<-lr$GPP
fpar<-lr$Fpar
psnnet<-lr$PsnNet
bdoy<-dimnames(doy[[1]])[[2]]
divide<-10000 # divide values by 10000

years<-sort(unique(as.integer(substr(c(gimms_doy,bdoy),1,4))))
years<-years[years<2017]

### mask snow or cloud values
# the ts with removed values does not play well with the SG filter
# and it gives really wild results with NDVI...

#rely<-rely%in%c(-1,2,3,255)
#r<-mask(r,rely,maskvalue=TRUE)
#rd<-mask(rd,rely,maskvalue=TRUE)

#registerDoParallel(6) 
#getDoParWorkers()

### build a list of matrices of each pixel ndvi value for each polygon, the for each is designed for multiple polygons, but here there is only one
#vndvi<-foreach(i=1:length(pol),.packages=c("raster","sp")) %do% {
#  extract(r,pol[i,])[[1]]  
#}

#vdoy<-foreach(i=1:length(pol),.packages=c("raster","sp")) %do% {
#  extract(rd,pol[i,])[[1]]  
#}

#vlai<-foreach(i=1:length(pol),.packages=c("raster","sp")) %do% {
#  extract(lai,spTransform(pol,CRS(proj4string(lai)))[i,])[[1]]  
#}

#vsnow<-foreach(i=1:length(pol),.packages=c("raster","sp")) %do% {
#  extract(snow,spTransform(pol,CRS(proj4string(lai)))[i,])[[1]]  
#}

#vgpp<-foreach(i=1:length(pol),.packages=c("raster","sp")) %do% {
#  extract(gpp,spTransform(pol,CRS(proj4string(gpp)))[i,])[[1]]  
#}


#vgpp<-foreach(i=1:length(pol),.packages=c("raster","sp")) %do% {
#  extract(brick(gpp[[1:10]]),spTransform(pol,CRS(proj4string(gpp)))[i,])[[1]]  
#}

#pol<-spTransform(pol,CRS(proj4string(lai)))

### faster velox version of preceding lines, but transforming to velox actually takes too long
#vx<-velox(r)
#v<-vx$extract(pol)

#vdx<-velox(rd)
#vd<-vdx$extract(pol)


### remove all NA values for some pixels
#v<-lapply(v,function(i){
#  k<-!apply(i,1,function(j){
#    all(is.na(j))  
#  })
#  i[k,]
#})

#vd<-lapply(vd,function(i){
#  k<-!apply(i,1,function(j){
#    all(is.na(j))  
#  })
#  i[k,]
#})


#############################################################
##### Compute metrics #######################################

# for the gimms data, I think the values are obtained using only two pixels, one of which does not overlap a lot with the ram pol
#plot(r[[1]])
#plot(pol,add=TRUE)



ts<-c("ndvi","evi","lai","gpp","snow","psnnet","fpar")
#ts<-c("evi")
cols<-list(ndvi="darkgreen",evi="chartreuse4",lai="chartreuse2",gpp="coral3",snow="azure3",psnnet="aquamarine3",fpar="aquamarine4")
lts<-list(ndvi="Normalized Difference Vegetation Index",evi="Enhanced Vegetation Index",lai="Leaf Area Index",gpp="Gross Primary Productivity",snow="8-Day Snow Cover",psnnet="Net Photosynthesis",fpar="Fraction of Photosynthetically Active Radiation")
#ts<-c("evi")

peak_cell<-lapply(ts,function(i){
    
  lapply(1:nrow(get(i)[[1]]),function(j){
    
    #i<-sample(ts,1)
    
    pl<-FALSE# for plotting or not
    pdfs<-FALSE
    cell<-(-sample(1:nrow(get(i)[[1]]),1))
    
    if(cell>0 && j!=cell){
      return(NA)
    }
    
    type<-i
    
    #################
    ### Values
    val<-get(type)[[1]][j,] # divide by 1 for gimms data and by 10000 for modis data
    #val<-gpp[[i]][j,]/divide # divide by 1 for gimms data and by 10000 for modis data
    if(all(is.na(val))){
      val[seq_along(val)]<-0  
    }
    dates<-dimnames(get(type)[[1]])[[2]]
    if(any(is.na(val))){
      dates<-dates[!is.na(val)]  
      val<-val[!is.na(val)]  
    }
    if(type%in%c("ndvi","evi")){
      val<-val/divide
    }
    if(type=="lai"){ # values over 45 are eliminated (some values in winter are abnormally high)
      dates<-dates[val<45]  
      val<-val[val<45]
    }
    if(type=="fpar"){ # values over 90 are eliminated (some values in winter are abnormally high)
      dates<-dates[val<90]  
      val<-val[val<90]
    }
    if(type=="snow"){
      val[val==50]<-NA
      val[val==200]<-0
      val[val==25]<-1
      dates<-dates[!is.na(val)] 
      val<-val[!is.na(val)]
    }
    jul<-as.integer(format(as.Date(dates),"%j"))
    #jul<-doy[[i]][j,] # for ndvi-evi data
    names(jul)<-dates
    #jul2<-as.integer(sapply(strsplit(names(jul),"_"),"[",4))
    jul2<-as.integer(format(as.Date(dates),"%j"))
    k<-is.na(jul)
    if(any(k)){
      jul[k]<-jul2[k] # missing values in precise julian days are given the beginning of the block  
    }
    name<-unname(as.Date(jul2nd(as.integer(substr(dates,1,4)),jul,jul2),origin="1970-01-01")) # blocks are ordered, but not necessarily precise dates
    names(val)<-name
    o<-order(name)
    val<-val[o]
    
    
    

    ## Savitsky-Golay filter
    #s1<-sgolayfilt(na.spline(val),n=ifelse(length(val)<41,11,41),p=3,m=1)
    #names(s1)<-names(val)
    #pos_up<-unlist(findminmax(s1,n=1,beg="04-01",end="07-01"))
    #pos_do<-unlist(findminmax(s1,n=1,beg="09-01",end="12-01",max=FALSE))
    #sg_up<-as.Date(names(s1)[pos_up])
    #sg_do<-as.Date(names(s1)[pos_do])
    #names(sg_up)<-substr(sg_up,1,4)
    #names(sg_do)<-substr(sg_do,1,4)
    #s0<-sgolayfilt(na.spline(val),n=ifelse(length(val)<21,11,21),p=3,m=0)
    
    if(!j%%20)
      print(paste(i,j))
    if(pl || pdfs){  
      if(j==cell){
        pdf(file.path("C:/Users/rouf1703/Downloads",paste0(i,".pdf")),width=18,height=6,pointsize=4)
      }
      plot(as.Date(names(val)),val,ylim=rangeinc(val),xaxt="n",col=alpha(cols[[type]],0.5),pch=16,ylab=i)
      axis.Date(1,at=as.Date(paste0(substr(names(val),1,4),paste0("-",formatC(1:12,width=2,flag=0),"-01"))),format="%y-%b",las=2)
      axis.Date(3,at=as.Date(paste0(substr(names(val),1,4),paste0("-",formatC(1:12,width=2,flag=0),"-01"))),format="%y-%b",las=2)
      #lines(as.Date(names(val)),s0)
      #points(as.Date(names(val)),s1*7,col="red",cex=0.5)
      abline(0,0)
    }
    
    ### Logistic curve
    mLog<-fitDLog(val[!is.na(val)],plot=pl,type=type,col=alpha(cols[[type]],0.85)) # new up down version (green)
    
    ### not sur eif this is essential
    #y<-setdiff(years,names(mLog)) # we complete the logistic, and because of merge we don't need to complete the SG I think
    #if(any(y)){
    #  mLog<-c(rep(NA,length(y)),mLog)
    #  names(mLog)[1:length(y)]<-y
    #  mLog<-mLog[order(names(mLog))]
    #}
    
    #browser()
    #if(divide==1){
    #  mLog<-mLog[-1] # take out first year for gimms
    #}
    
    ### peaks for the simple 
    #log_up<-as.Date(sapply(mLog,function(k){unname(k["S"])})) # this is only valid for a simple logistic, but not for the double
    #log_do<-as.Date(sapply(mLog,function(k){unname(k["A"])}))
    
    ###
    log_up<-as.Date(sapply(mLog,function(k){if(identical(k,NA)){
      NA
    }else{
      vals<-seq(as.list(k)$S-150,as.list(k)$S+150,by=0.5)
      vals[which.max(do.call("DLog1",c(list(x=vals),as.list(k))))]
    }}))
    log_do<-as.Date(sapply(mLog,function(k){if(identical(k,NA)){
      NA
    }else{
      vals<-seq(as.list(k)$A-150,as.list(k)$A+150,by=0.5)
      vals[which.min(do.call("DLog1",c(list(x=vals),as.list(k))))]
    }}))
    
    
    
    #ans<-rbind(sg_up,log_up,sg_do,log_do)
    #ans<-list(sg_up=sg_up,log_up=log_up,sg_do=sg_do,log_do=log_do)
    ans<-list(log_up=log_up,log_do=log_do)
    ans<-lapply(seq_along(ans),function(k){
      y<-data.frame(year=names(ans[[k]]),date=ans[[k]],stringsAsFactors=FALSE)  
      names(y)[2]<-names(ans)[k]
      y
    })
    ans<-Reduce(function(x,y) merge(x,y,by=c("year"),all=TRUE),ans)

    if(all(val==1)){
      ans[,2:3]<-NA
    }
    #browser()
    if(pl){  
      # check ordering of dates
      #h_up<-sapply(mLog,function(k){k["c"]})+sapply(mLog,function(k){k["Asym_up"]})/2
      #h_do<-sapply(mLog,function(k){k["c"]})+sapply(mLog,function(k){k["Asym_up"]})+sapply(mLog,function(k){k["Asym_do"]})/2
      
      #h_up<-sapply(mLog,function(k){if(identical(k,NA)){NA}else{do.call("DLog",c(x=as.list(k)$S,as.list(k)))}})
      #h_do<-sapply(mLog,function(k){if(identical(k,NA)){NA}else{do.call("DLog",c(x=as.list(k)$A,as.list(k)))}})
      
      h_up<-sapply(seq_along(mLog),function(k){if(identical(mLog[[k]],NA)){NA}else{do.call("DLog",c(x=as.integer(log_up[[k]]),as.list(mLog[[k]])))}})
      h_do<-sapply(seq_along(mLog),function(k){if(identical(mLog[[k]],NA)){NA}else{do.call("DLog",c(x=as.integer(log_do[[k]]),as.list(mLog[[k]])))}})
      
      axis.Date(1,at=log_up,las=2,cex.axis=0.7,col.axis=alpha("green4",0.5),format="%b-%d",line=-3)
      axis.Date(1,at=log_do,las=2,cex.axis=0.7,col.axis=alpha("brown",0.5),format="%b-%d",line=-3)
      #points(ans$sg_up,h_up,col="red",pch=16)
      #points(ans$sg_do,h_do,col="red",pch=16)
      points(ans$log_up,h_up,col="green4",pch=15,cex=1.5)
      points(ans$log_do,h_do,col="brown",pch=15,cex=1.5)
      
      abline(median(val,na.rm=TRUE),0)
      
    }
    
    if(pl || pdfs){  
      if(j==cell){
        dev.off()
      }
    }
    
    ### plotting cells
    #w<-which(!is.na(ocells[,1]))
    #plot(ram)
    #plot(spTransform(cents[w,],CRS(proj4string(ram))),add=TRUE,pch=1)
    #plot(spTransform(cents[w[j],],CRS(proj4string(ram))),add=TRUE,pch=16,col="red",cex=2)
    
    #tmap_mode("plot")
    #tm_shape(pol) +
    #  tm_polygons(col="red",alpha=0.3) +
    #  tm_borders(col="red",lwd=5) +
    #  tm_shape(r[[1]]) +
    #  tm_raster() +
    #  tm_shape(cent[w[j],])+
    #  tm_dots()+
    #  tm_view(basemaps=c("Esri.WorldImagery"))
    
    ans
  })
})
names(peak_cell)<-ts


#metrics<-row.names(peak_cell[[1]][[1]]) # peak_cell is a list (1 for pol) of list of pixels


#peaks<-sapply(metrics,function(k){
#  lapply(peak_cell,function(i){
#    ii<-do.call("rbind",lapply(i,function(j){j[row.names(j)==k,]}))
#    as.Date(colMeans(ii,na.rm=TRUE),origin="1970-01-01")  
#  })
#})
#names(peaks)<-metrics



peaks<-lapply(ts,function(i){
  x<-rbindlist(peak_cell[[match(i,ts)]]) 
  cols<-names(x)[-1]
  x[,(cols):=list(as.Date(log_up),as.Date(log_do))]
  #x<-x[1:5,]
  x<-x[,lapply(.SD,function(j){mean.Date(j,na.rm=TRUE)}),by=year,.SDcols=cols] # problem here!
  names(x)[-1]<-paste(i,names(x)[-1],sep="_")
  x
})
peaks<-Reduce(function(x,y){merge(x,y,all=TRUE)},peaks)

j<-do.call("data.frame",lapply(peaks[,-1],function(i){as.integer(format(i,"%j"))}))
names(j)<-paste0(names(peaks)[-1],"_jul")
peaks2<-cbind(peaks,j)

fwrite(peaks2,"S:/NDVI/MODIS/peaks_ram.csv",row.names=FALSE)



#################################
### graph all time series #######

plot(as.integer(format(peaks[,2],"%j")),peaks$year,xlim=range(as.integer(format(peaks[,-1],"%j")),na.rm=TRUE),type="n",xaxt="n",yaxt="n",xlab="Dates",ylab="Year")

dates<-seq.Date(as.Date("2000-01-01"),as.Date("2001-01-01"),by="day")
dates<-dates[grep("-01|-10|-20",substr(dates,8,10))]
j<-as.integer(format(dates,"%j"))
axis(1,at=j,labels=format(dates,"%b %d"),las=2)
segments(j,-2000,j,4000,lty=3,col=gray(0,0.1))
#segments(j[grep("-01",dates)],-2000,j[grep("-01",dates)],4000,lty=3,col=gray(0,0.25))

invisible(lapply(2:ncol(peaks),function(i){
  tsval<-sapply(strsplit(names(peaks)[i],"_"),"[",1)
  lines(as.integer(format(peaks[,..i],"%j")),peaks$year,type="b",lwd=5,col=alpha(cols[[tsval]],0.75))   
}))

legend("top",cex=1.15,lwd=5,legend=paste(names(cols),"-",lts),col=alpha(unlist(cols),0.75),bty="n",inset=c(0.1,0.00))
axis(2,at=peaks$year,las=2)

##########################################
### map values to rasters ################

w<-which(!is.na(ocells[,1]))
gridi<-v$as.RasterStack()[[1]]

rl<-lapply(peak_cell,function(i){
  years<-i[[1]]$year
  l<-lapply(years,function(j){
    grid<-gridi  
    vals<-rep(NA,ncell(grid))
    vals[w]<-sapply(i,function(k){
      as.integer(format(k[k$year==j,"log_up"],"%j"))
    })
    #vals[w]<-1:length(peak_cell[[1]])
    grid[]<-vals[t(matrix(1:ncell(grid),ncol=ncol(grid)))]
    grid
  })
  res<-stack(l)
  names(res)<-years
  res<-crop(res,pol)
  res<-projectRaster(from=res,crs=CRS(proj4string(ram))) # this changes the values in the cells
  res
})
names(rl)<-names(peak_cell)

levelplot(rl$gpp)

test<-do.call("rbind",peak_cell$evi)
range(as.integer(format(test[,2],"%j")),na.rm=TRUE)
hist(as.integer(format(test[,2],"%j")),na.rm=TRUE)

#tmap_mode("view")
#tm_shape(pol) +
#  tm_polygons(col="red",alpha=0.3) +
#  tm_borders(col="red",lwd=5) +
#  tm_shape(grid) +
#  tm_raster(alpha=0.5) +
#  tm_view(basemaps=c("Esri.WorldImagery"))



###############################
### PCA #######################
up<-as.data.frame(peaks)[,grep("_up",names(peaks))]
row.names(up)<-peaks$year

up[]<-lapply(up,function(i){
  as.integer(format(i,"%j"))  
})
up<-up[apply(up,1,function(i){!any(is.na(i))}),]

cor(up,use="complete.obs")

m<-rda(as.matrix(up),scale=TRUE)
biplot(m,scaling=2)

#
#peak_log_up<-lapply(peak_cell,function(i){
#  ii<-do.call("rbind",lapply(i,function(j){j[1,]}))
#  as.Date(colMeans(ii,na.rm=TRUE),origin="1970-01-01")  
#})
#peak_sg_up<-lapply(peak_cell,function(i){
#  ii<-do.call("rbind",lapply(i,function(j){j[2,]}))
#  as.Date(colMeans(ii,na.rm=TRUE),origin="1970-01-01")  
#})



peak_gimms<-peaks # change here to obtain both series
#peak_modis<-peaks



#tmap_mode("view")
#tm_shape(pol)+tm_borders(lwd=2,col="red")+tm_layout(basemaps=c("Esri.WorldImagery"))


### gimms and modis compare (peak1 peak2)



j_gimms_log<-as.Date(as.integer(format(peak_gimms$log_up,"%j")))
j_gimms_sg<-as.Date(as.integer(format(peak_gimms$sg_up,"%j")))
j_modis_log<-as.Date(as.integer(format(peak_modis$log_up,"%j")))
j_modis_sg<-as.Date(as.integer(format(peak_modis$sg_up,"%j")))
ylim<-range(as.Date(c("1970-04-01","1970-07-01")))


#j1g<-scale(j1g)[,1]
#j2g<-scale(j2g)[,1]
#j1m<-scale(j1m)[,1]
#j2m<-scale(j2m)[,1]
#ylim<-range(c(j1g,j2g))

plot(as.integer(substr(peak_gimms$log_up,1,4)),j_gimms_log,pch=16,col=alpha("red",0.3),ylim=ylim,type="l",las=2,lwd=4,xlim=c(1980,2016))
lines(as.integer(substr(peak_gimms$sg_up,1,4)),j_gimms_sg,pch=16,col=alpha("red",0.3),lwd=4,lty=2)
lines(as.integer(substr(peak_modis$log_up,1,4)),j_modis_log,pch=16,col=alpha("blue",0.3),lwd=4)
lines(as.integer(substr(peak_modis$sg_up,1,4)),j_modis_sg,pch=16,col=alpha("blue",0.3),lwd=4,lty=2)
abline(0,0)


gimmsSG<-peak_gimms$log_up[match(1981:2016,substr(peak_gimms$log_up,1,4))]
gimmsLO<-peak_gimms$sg_up[match(1981:2016,substr(peak_gimms$sg_up,1,4))]
modisSG<-peak_modis$log_up[match(1981:2016,substr(peak_modis$log_up,1,4))]
modisLO<-peak_modis$sg_up[match(1981:2016,substr(peak_modis$sg_up,1,4))]

res<-data.frame(years=1981:2016,gimmsSG,gimmsLO,modisSG,modisLO)
res$modisSG[res$years==2000]<-NA #

res2<-ddply(res,.(years),function(i){format(i[-1],"%j")})
names(res2)[2:ncol(res2)]<-paste0(names(res2)[2:ncol(res2)],"jul")

res<-merge(res,res2)
#fwrite(res,"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc/greenup_ts_all3.csv",row.names=FALSE,sep=";")




####################################################################################################
####################################################################################################
### verifications

x1<-fread("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc/greenup_ts_all.csv",sep=";")
x2<-fread("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc/greenup_ts_all2.csv",sep=";")

plot(x1$years,x1$modisSGjul,type="n",ylim=c(80,160))
lines(x1$years,x1$modisSGjul)
lines(x2$years,x2$modisSGjul)

lines(x1$years,x1$gimmsSGjul)
lines(x2$years,x2$gimmsSGjul)



