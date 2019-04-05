

library(mapview)
library(raster)
library(dplyr)
library(adehabitatHR)
library(sp)
library(rgdal)
library(rgeos)
library(velox)
library(scales)
library(doSNOW)
library(foreach)
library(FRutils)
library(sf)
library(rasterVis)


################################
### Import roo data from 2016
d16<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/LMontana/2016Spatial_181001.csv",stringsAsFactors=FALSE)
d16<-d16[!is.na(d16$X),] # remove NAs
d16<-d16[d16$X>36612 & d16$X<38590 & d16$Y<88929 & d16$Y>88032,] # keep the study area
d16$X<-paste0("4",d16$X) # add first number for east in AGD66
d16$Y<-paste0("56",d16$Y) # add first two numbers for north in AGD66
d16$X<-as.numeric(d16$X) # transform column X in numeric
d16$Y<-as.numeric(d16$Y) # transform column Y in numeric
d16<-d16[!is.na(d16$X),] # create a new dataset without the NAs
coordinates(d16)<-~X+Y # transform object data into a SpatialPoints data.frama
proj4string(d16)<-"+init=epsg:20255" # define the coordinate system of the Spatial.Points data.frame
roos<-d16[sample(1:nrow(d16),5000),]
prom<-gConvexHull(roos)
mapview(spTransform(roos[1:100,],CRS("+init=epsg:4326")) ) # transformation coordinate from AGD66 into latitude-longitude


#########################################
### to extract dates from MODIS data
# x is a file name ending with "_2002_342.tif"
getDates<-function(x){
  x<-gsub(".tif","",x)
  n<-nchar(x)
  s<-substr(x,n-7,n)
  s<-strsplit(s,"_")
  as.Date(paste0(sapply(s,"[",1),"-01-01"))+as.integer(sapply(s,"[",2))-1
}

##########################################
### path to rasters MODIS files
paths<-c("D:/Roos/Gross_PP_8Days_500m_v6/GPP","D:/Roos/VI_16Days_250m_v6/NDVI","D:/Roos/LAI_8Days_500m_v6/Lai")
cols<-c("coral3","darkgreen","chartreuse2")


##########################################
### read and crop MODIS files in parallel
cl <- makeCluster(3)
registerDoSNOW(cl)
lr<-foreach(i=1:length(paths),.packages=c("raster","sp")) %dopar% {
  x<-list.files(paths[i],full.names=TRUE)
  x<-x#[1:100]
  r<-stack(x,quick=TRUE)
  ext<-spTransform(prom,CRS(proj4string(r)))
  cells<-cellFromPolygon(r,ext,weights=TRUE)[[1]]
  r2<-rasterFromCells(r,cells[,"cell"])
  r<-crop(r,extent(r2))
  names(r)<-getDates(x)
  r
}


#########################################################
### plot rasters with locations in sinusoidal projection
par(mfrow=c(2,2),mar=c(1,1,1,1))
sapply(seq_along(lr),function(i){
  r<-lr[[i]]
  plot(r[[1]],col=alpha(colo.scale(100:1,c("white",cols[i])),0.75))
  plot(prom,add=TRUE)
  plot(spTransform(roos,CRS(proj4string(r))),add=TRUE,pch=16,col=gray(0.5,0.25))
  text(SpatialPoints(xyFromCell(r,1:ncell(r))),1:ncell(r))
})


#########################################################
### mapview pixels
pixelsNDVI<-rasterToPolygons(lr[[2]])
pixelsNDVI<-st_as_sf(spTransform(pixelsNDVI,CRS(proj4string(roos))))
pixelsLAI<-rasterToPolygons(lr[[1]])
pixelsLAI<-st_as_sf(spTransform(pixelsLAI,CRS(proj4string(roos))))
roos2<-st_as_sf(roos[1:500,])
#b<-c("Esri.WorldShadedRelief", "OpenStreetMap.DE")
b<-mapviewGetOption("basemaps")
mapview(list(pixelsNDVI,pixelsLAI,roos2),map.types=b)


plot(st_geometry(pixelsLAI),border="darkgreen",lwd=8)
plot(st_geometry(pixelsNDVI),border="green",lwd=2,add=TRUE)
plot(st_geometry(roos2),col=gray(0.5,0.5),lwd=1,pch=16,add=TRUE)


#########################################################
### rasterVis

p.strip <- list(cex=0.1, lines=-1, col="transparent")
levelplot(subset(lr[[2]],order(names(lr[[2]]))),col.regions=rev(colo.scale(1:100,c("white","green","green4","black"))),cuts=99,layout=c(24,24),par.settings=list(strip.background = list(col = "transparent"),strip.border = list(col = 'transparent'),axis.line=list(col="transparent")),scales=list(col="black",tck = c(1,0)),par.strip.text=p.strip)


#########################################################
### extract values in matrices
le<-foreach(i=1:length(lr),.packages=c("raster","sp")) %dopar% {
  r<-lr[[i]]
  e<-extract(r,extent(r))
  cells<-cellFromPolygon(r,prom,weights=TRUE)[[1]]
  dimnames(e)[[2]]<-gsub("X","",gsub("\\.","-",names(r)))
  e<-e[,order(dimnames(e)[[2]])]
  e
}


#########################################################
### plot each ts in a separate graph by using the pixel with the most location
# the 4th graph is for the vegetation data
par(mfrow=c(4,1),mar=c(3,5,1,1))
sapply(seq_along(le),function(i){
  e<-le[[i]]
  
  g<-rasterToPolygons(lr[[i]][[1]])
  roos2<-spTransform(roos,CRS(proj4string(lr[[1]])))
  o<-over(g,roos2,returnList=TRUE)
  cells<-order(sapply(o,nrow),decreasing=TRUE)[1:4] # takes the x cells with the most locations
  cells<-c(1) # to test specific cells
  print(cells)
  
  xaxt<-seq.Date(as.Date("2008-01-01"),as.Date("2018-01-01"),length.out=2)
  plot(xaxt,c(0,0),ylim=quantile(e[cells,],probs=c(0.05,0.98),na.rm=TRUE),xaxt="n",xlab="",type="n",ylab="")
  titre<-toupper(sapply(strsplit(paths,"/"),tail,1))
  legend("topleft",legend=NA,title=titre[i],col=cols[i],cex=2,bty="n",title.col=cols[i])
  axis.Date(1,at=as.Date(paste0(rep(2008:2018,each=11),"-",formatC(2:12,width=2,flag=0),"-01")),format="%b",las=2,cex.axis=0.8)
  axis.Date(1,at=as.Date(paste0(rep(2008:2018,each=1),"-",formatC(1,width=2,flag=0),"-01")),format="%y",las=1,cex.axis=1.5,font=2)
  sapply(cells,function(j){
    lines(as.Date(dimnames(e)[[2]]),e[j,],xaxt="n",col=alpha(cols[i],0.5),pch=16,type="b",lwd=2)
  })
  n<-data.frame(date=rep(as.integer(as.Date(dimnames(e)[[2]])),each=nrow(e[cells,,drop=FALSE])),ndvi=as.vector(e[cells,]))
  
  smooth<-loess(ndvi~date,data=n,span=0.05)
  #smooth<-gam(ndvi~s(date,k=100),data=n)
  val<-as.integer(seq.Date(as.Date(min(dimnames(e)[[2]])),as.Date(max(dimnames(e)[[2]])),by=5))
  s<-predict(smooth,data.frame(date=val))
  lines(as.Date(val,origin="1970-01-01"),s,xaxt="n",col=alpha(cols[i],0.75),pch=16,type="l",lwd=4)
  lines(v,rescale(p$fit,to=c(max(min(s),par("usr")[3]),min(max(s),par("usr")[4]))))
  
})


################################################################################
### add the vegetation data

### this is for sourcing the relevant code in veg.R
#file<-"C:/Users/rouf1703/Documents/UdeS/Consultation/PToni/veg.R"
#veg<-readLines(file)
#g<-grep("### vegetation growth curve",veg)
#source(textConnection(veg[1:g]))

### the actual data
xaxt<-seq.Date(as.Date("2008-01-01"),as.Date("2018-01-01"),length.out=2)
plot(xaxt,c(0,0),ylim=quantile(d$daily,probs=c(0.0,0.99),na.rm=TRUE),xaxt="n",xlab="",type="n",ylab="")
points(d$date,d$daily,col=gray(0.5,0.5),xaxt="n",ylab="Daily Growth",xlab="n")
axis.Date(1,at=as.Date(paste0(rep(2008:2018,each=11),"-",formatC(2:12,width=2,flag=0),"-01")),format="%b",las=2,cex.axis=0.8)
axis.Date(1,at=as.Date(paste0(rep(2008:2018,each=1),"-",formatC(1,width=2,flag=0),"-01")),format="%y",las=1,cex.axis=1.5,font=2)
lines(l$date,l$daily,lwd=2,col=gray(0.5,0.5))
legend("topleft",pch=c(1,NA),lwd=c(NA,2),col=rep(gray(0.5,0.5),2),legend=c("Daily Growth for a quadrat","Mean Daily Growth"),bty="n",inset=c(0.05,0.05))
lines(v,p$fit)
#lines(v,p$fit+1.96*p$se.fit,lty=3)
#lines(v,p$fit-1.96*p$se.fit,lty=3)






##############################
### LAI from MODIS/VIIRS
##############################

#cl <- makeCluster(6)
#registerDoSNOW(cl)
#v<-foreach(j=1:length(x),.packages=c("velox")) %dopar% {
#  velox(x[j])
#}
#v<-velox(v)


untar("C:/Users/rouf1703/Downloads/GTiff.tar.gz",list=TRUE)
untar("C:/Users/rouf1703/Downloads/GTiff.tar.gz",exdir="C:/Users/rouf1703/Downloads/LAI")

l<-sort(list.files("C:/Users/rouf1703/Downloads/LAI",pattern="_Lai_500m",full.names=TRUE))
r<-stack(l[1])

cl <- makeCluster(6)
registerDoSNOW(cl)
v<-foreach(j=1:length(l),.packages=c("velox")) %dopar% {velox(l[j])}
v<-velox(v)
e<-v$extract(as(extent(r),"SpatialPolygons"))[[1]]
lai<-e[1,]
laid<-as.Date(paste0(substr(l,43,46),"-01-01"))+as.integer(substr(l,47,49))-1

points(laid,rescale(lai,to=c(0.001,0.003)))

plot(r[[1]])
b<-as(extent(r),"SpatialPolygons")
proj4string(b)<-proj4string(r)
plot(spTransform(b,CRS("+init=epsg:4326")))
plot(spTransform(prom,CRS(proj4string(r))),add=TRUE)

#e<-extract(r,extent(r))




