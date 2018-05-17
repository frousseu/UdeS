
library(sp)
library(rgdal)
library(scales)
library(gstat)
library(geoR)
library(geostatsp)
library(INLA)
library(rgeos)
library(FRutils)
library(RColorBrewer)
library(brinla)
library(visreg)
library(quantreg)

load("~/UdeS/Consultation/GStetcher/Doc/LLF_size.RData")

d<-llf.size

ds<-d
coordinates(ds)<-~Longitude+Latitude
proj4string(ds)<-"+init=epsg:4326"

prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
ds<-spTransform(ds,CRS(prj))

#plot(ds,col=alpha(ifelse(ds$PA==1,"red","blue"),0.25),pch=16)

d$high_name<-as.factor(d$high_name)
d$logArea<-log(d$Area)

m <- glm (logArea ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = gaussian(link = "identity"), data = d)
m <- rq (logArea ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = d, tau=0.95)
#m <- glm (Area ~ 1, family = gaussian(link = "log"), data = d)
ds$resid<-resid(m)


par(mfrow=c(3,3))
visreg(m,trans=exp)
par(mfrow=c(1,1))

coords <- coordinates(ds)
v<-variog(coords=coords,data=resid(m),breaks=seq(0,200000,by=1000),max.dist=200000)
#v<-variog(coords=coords,data=resid(m))
#v<-variogram(resid~1,data=ds)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b") 


swe <- raster::getData("GADM", country = "SWE", level = 1)
swe<-spTransform(swe,proj4string(ds))