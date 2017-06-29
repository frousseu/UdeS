
library(foreach)
library(doParallel)


### This script is for extracting ndvi/evi metrics from RasterStack objects using a set of regions defined by polygons

# RasterStacks may be of the MODIS or GIMMS type
# A raster of julian dates is assumed to accompany each value raster

jul2nd<-function(y,j){
  as.integer(as.Date(paste0(y,"-01-01")))+j-1
}


#############################################################
##### Get regions of interest ###############################

# A list of SpatialPolygonsDataFrame for which a metric is wanted

### Ram mountain
ram<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram")
pol<-gBuffer(ram,width=-0.01)
pol<-SpatialPolygonsDataFrame(pol,data.frame(id=1),match.ID=FALSE)


### Refuges Alberta BC
#z<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/YPoisson/Doc",layer="largezone_BC_Alberta")
#pol<-spTransform(z,CRS(proj4string(r)))


#############################################################
##### Get MODIS raster ######################################

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/MODIS/VI_16Days_250m_v6/Time_Series/RData/"

load(paste0(path,"MOD13Q1_MYD13Q1_NDVI_49_2000_1_2017_RData.RData"))
modis<-raster_ts #r5 is the initial raster used in plot ndvi3
load(paste0(path,"MOD13Q1_MYD13Q1_DOY_49_2000_1_2017_RData.RData"))
modis_jul<-raster_ts
rm(raster_ts)
l<-strsplit(names(modis),"_")
modis_doy<-as.Date(paste0(sapply(l,"[",3),"-01-01"))+as.integer(sapply(l,"[",4))-1

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


#############################################################
##### Get GIMMS raster ######################################

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/GIMMS/"

x<-list.files(path)
ts<-monthlyIndices(x, version = 1, timestamp = T)
gimms_doy<-as.Date(as.character(ts+round(c(diff(ts),17)/2,0)))
gimms<-rasterizeGimms(x=paste0(path,x),ext=pol,cores=6L) # clipping
names(gimms)<-as.character(gimms_doy)
gimms_jul<-stack(setValues(gimms,as.integer(format(rep(gimms_doy,each=ncell(gimms)),"%j"))))







#############################################################
##### Extract raster values for each region #################

r<-modis # determine series to use
rd<-modis_jul
doy<-modis_doy


registerDoParallel(6) 
getDoParWorkers()

# build a list of matrices of each pixel ndvi value for each polygon
v<-foreach(i=1:length(pol),.packages=c("raster","sp")) %dopar% {
  extract(r,pol[i,])[[1]]  
}

vd<-foreach(i=1:length(pol),.packages=c("raster","sp")) %dopar% {
  extract(rd,pol[i,])[[1]]  
}


#############################################################
##### Compute metrics #######################################



peak<-lapply(seq_along(v),function(i){
  sapply(1:nrow(v[[i]]),function(j){
    s0<-sgolayfilt(na.spline(v[[i]][j,]),n=11,p=3,m=0)
    s1<-sgolayfilt(na.spline(v[[i]][j,]),n=11,p=3,m=1)
    names(s1)<-as.integer(as.Date(jul2nd(substr(doy,1,4),vd[[i]][j,]),origin="1970-01-01"))
    #names(s1)<-ifelse(is.na(names(s1)),doy,names(s1))
    pos<-unlist(findminmax(s1,n=1,beg="03-01",end="07-01"))
    names(s1)[pos]
  })
})

peak<-lapply(peak,function(i){
  as.Date(colMeans(i),origin="1970-01-01")  
})


















