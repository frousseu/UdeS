
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

size<-llf.size

sizes<-size
coordinates(sizes)<-~Longitude+Latitude
proj4string(sizes)<-"+init=epsg:4326"

prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sizes<-spTransform(sizes,CRS(prj))

#plot(ds,col=alpha(ifelse(ds$PA==1,"red","blue"),0.25),pch=16)

size$high_name<-as.factor(size$high_name)
size$logArea<-log(size$Area)

m <- glm (logArea ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = gaussian(link = "identity"), data = size)
m <- rq (logArea ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size, tau=0.95)
#m <- glm (Area ~ 1, family = gaussian(link = "log"), data = d)
sizes$resid<-resid(m)


par(mfrow=c(3,3))
visreg(m)
par(mfrow=c(1,1))

coords <- coordinates(sizes)
v<-variog(coords=coords,data=resid(m),breaks=seq(0,200000,by=1000),max.dist=200000)
#v<-variog(coords=coords,data=resid(m))
#v<-variogram(resid~1,data=ds)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b") 


swe <- raster::getData("GADM", country = "SWE", level = 1)
swe<-spTransform(swe,proj4string(sizes))