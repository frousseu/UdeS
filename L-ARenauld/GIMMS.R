
library(gimms)
library(ncdf4)
library(rgdal)
library(rasterVis)
library(rgeos)
library(FRutils)
library(signal)
library(scales)

gimms_files_v1 <- updateInventory()
tail(gimms_files_v1, 4)

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/GIMMS/"

#gimms_files <- downloadGimms(x = as.Date("2000-01-01"),y = as.Date("2001-09-15"),dsn=path,cores=2L)

## Extent for clipping

# get extent from MODIS raster on ram
pathMODIS<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/MODIS/VI_16Days_250m_v6/Time_Series/RData/"
l<-list.files(pathMODIS)[1]
load(paste0(pathMODIS,l))
bb<-attributes(extent(raster_ts))
ramr<-bbox2pol(c(bb$xmin,bb$xmax,bb$ymin,bb$ymax))
rm(raster_ts)

### ram
shp <- getData("GADM", country = "DEU", level = 0, path = tmpDir())
ram<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram")
plot(ram)
ram<-gBuffer(ram,width=0.25)

### yoanna's region
pathshp<-"C:/Users/rouf1703/Documents/UdeS/Consultation/YPoisson/Doc"
z<-readOGR(pathshp,layer="largezone_BC_Alberta")
z<-spTransform(ram,CRS(proj4string(r)))





## Rasterize without quality control

theme.novpadding <-
  list(layout.heights =
         list(top.padding = 0,
              main.key.padding = 0,
              key.axis.padding = 0,
              axis.xlab.padding = 0,
              xlab.key.padding = 0,
              key.sub.padding = 0,
              bottom.padding = 0),
       layout.widths =
         list(left.padding = 0,
              key.ylab.padding = 0,
              ylab.axis.padding = 0,
              axis.key.padding = 0,
              right.padding = 0))

x<-list.files(path)
ts<-monthlyIndices(x, version = 1, timestamp = T)
doy<-as.Date(as.character(ts+round(c(diff(ts),17)/2,0)))
r <- rasterizeGimms(x = paste0(path,x),ext = z,cores=6L) # clipping
plot(r[[1]])
spplot(subset(r,13:dim(r)[[3]]),layout=c(24,30))
lines(ram)
p.strip <- list(cex=0.1, lines=-1, col="transparent")
levelplot(subset(aggregate(r,10),(13:dim(r)[[3]])[1:(24*34)]),col.regions=rev(terrain.colors(101)),cuts=100,layout=c(24,34),par.settings=list(strip.background = list(col = "transparent"),strip.border = list(col = 'transparent'),axis.line=list(col="transparent")),scales=list(col="black",tck = c(1,0)),par.strip.text=p.strip)


### show all years and dates
png("C:/Users/rouf1703/Documents/UdeS/Consultation/YPoisson/Doc/ndvi.png",width=7,height=10.3,units="in",res=600)
m<-(13:dim(r)[[3]])-12
layout(matrix(m,ncol=24,byrow=FALSE))
r2<-aggregate(r,2)
par(mar=c(0,0,0.025,0),oma=c(0,0.2,0,0))
for(i in m){
  r3<-r2[[i+12]]
  brks <- seq(-0.3,1,by=0.01)
  nbrks <- length(brks)-1
  plot(r3,xaxt="n",yaxt="n",main="",axes=FALSE,legend=FALSE,box=FALSE,breaks=brks,col=rev(terrain.colors(nbrks)),lab.breaks=brks,zlim=c(0,1))
  mtext(doy[i+12],side=3,cex=0.125,line=-0.4,col="white",xpd=TRUE,outer=FALSE)
  if(i==1){plot(z,add=TRUE,lwd=0.05,border=alpha("black",0.25),col=alpha("grey50",0.25))}
}
dev.off()




years<-sort(unique(substr(doy,1,4)))
ee<-extract(r,z)
#e<-ee[[65]]
e<-do.call("rbind",ee)
dimnames(e)[[2]]<-as.character(doy)
plot(doy,e[1,],ylim=c(-0.2,1),type="n",xaxt="n")
axis.Date(1,at=seq(min(doy),max(doy),by="4 month"), format="%Y-%m-%d",las=2,cex.axis=0.5)
abline(0,0)

peak<-invisible(lapply(1:nrow(e[1:min(c(nrow(e),1000)),]),function(i){
  
  s0<-sgolayfilt(e[i,],n=11,p=3,m=0)
  s1<-sgolayfilt(e[i,],n=11,p=3,m=1)
  
  names(s1)<-dimnames(e)[[2]]
  
  #lines(doy,s0,col=alpha("black",0.01))
  #lines(doy,s1,col=alpha("black",0.01))
  
  pos<-unlist(findminmax(s1,n=1,beg="03-01",end="07-01"))
  #peak<-seq.Date(as.Date("2007-01-01"),as.Date("2007-12-31"),by=1)[round(mean(as.integer(format(doy[pos],"%j"))),0)]
  
  doy[pos]
  
  #invisible(lapply(pos,function(i){
  #  lines(rep(doy[i],2),c(-0.2,1),lty=2,col=alpha("red",0.01))
  #}))
  

}))

peak<-as.Date(colMeans(do.call("rbind",peak)),origin="1970-01-01")


## Rasterize with quality control
rq<-rasterizeGimms(x=gimms_files,ext=shp,keep=0)  # quality control
plot(rq[[1]])
lines(shp)
