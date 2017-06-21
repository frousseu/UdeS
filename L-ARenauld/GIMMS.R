
library(gimms)
library(ncdf4)
library(rgdal)
library(rasterVis)
library(rgeos)
library(FRutils)
library(signal)

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


shp <- getData("GADM", country = "DEU", level = 0, path = tmpDir())
ram<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram")
plot(ram)
ram<-gBuffer(ram,width=0.25)





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
r <- rasterizeGimms(x = paste0(path,x),ext = ramr) # clipping
plot(r[[1]])
spplot(subset(r,13:dim(r)[[3]]),layout=c(24,30))
lines(ram)
p.strip <- list(cex=0.1, lines=-2, col="transparent")
levelplot(subset(r,13:dim(r)[[3]]),col.regions=rev(terrain.colors(101)),cuts=100,layout=c(24,34),par.settings=list(strip.background = list(col = "transparent"),strip.border = list(col = 'transparent'),axis.line=list(col="transparent")),scales=list(col="black"),par.strip.text=p.strip)


lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)
levelplot(subset(r,1:12),par.settings=list(axis.line=list(col="transparent")),scales=list(col="black"))


e<-extract(r,ramr)[[1]]
e<-e[,seq_len(ncol(e))[1:828]]
plot(0,0,xlim=c(0,ncol(e)),ylim=c(-0.2,1),type="n")
lapply(1:nrow(e),function(i){
  points(e[i,],col=i) 
  lines(sgolayfilt(e[i,],n=11,p=3),col=i)
})


## Rasterize with quality control
rq <- rasterizeGimms(x = gimms_files,
                                ext = shp, # clipping
                                keep = 0)  # quality control
plot(rq[[1]])
lines(shp)
