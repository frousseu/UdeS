
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
library(visreg)
library(readxl)
library(mgcv)
library(qgam)
library(raster)
library(rasterVis)
library(data.table)
library(plyr)
library(alphahull)
library(concaveman)
library(mapview)
library(fasterize)
library(sf)
library(velox)
library(glmmTMB)
library(viridis)
library(foreach)
library(doParallel)


d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/BD.xlsx"))
d$Long<-d$LongdecSatScan
d$Lat<-d$LatdecSatScan
d$Long<-d$Long+0.3
d$Lat<-d$Lat-0.2
d$Site<-ifelse(nchar(d$Site)==6,sapply(strsplit(d$Site,"(?<=.{3})", perl = TRUE),paste,collapse=" "),d$Site)
#unique(d$Site)

#d$date<-as.Date(d$Day)
d$year<-d$Annee
#d$jul<-as.integer(format(d$date,"%j"))


d<-d[order(d$Site,d$Annee_Week),]

l<-split(d,d$year)

ds<-d
coordinates(ds)<-~Long+Lat
proj4string(ds)<-"+init=epsg:4326"

prj<-"+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
ds<-spTransform(ds,CRS(prj))

#mapview(ds,zcol="Site")

l<-list.files("C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/lcc2000",pattern=".shp",full.names=TRUE)
s<-lapply(l,st_read)
s<-do.call("rbind",s)
s<-st_transform(s,crs=26918)

#https://open.canada.ca/data/en/dataset/97126362-5a85-4fe0-9dc2-915464cfdbb7#
cla<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/lcc2000class.xlsx"))

clan<-sapply(strsplit(cla[,1],"-"),"[",1)
cla<-strsplit(cla[,1],"-") %>% 
  sapply("[",-1) %>% 
  sapply(paste0,collapse="") %>% 
  sapply(function(i){gsub('[[:digit:]]+', '', i)}) %>% 
  sapply(function(i){gsub("  "," ",i)}) %>% 
  sapply(function(i){gsub("  ","",i)}) %>% 
  unname


cla[grep("Conifé|Mixte|Feuillu|Forêt",cla)]<-"Natural"
cla[grep("Cultures|Prairies|agricoles",cla)]<-"Agriculture"
cla[grep("humide",cla)]<-"Wet"
cla[grep("Arbustes|arbustes|herbacées",cla)]<-"Natural"
cla[grep("découvert|Stérile|Bryo",cla)]<-"Exposed"
cla[grep("développées",cla)]<-"Urban"
cla[grep("Ombre",cla)]<-NA

s$cover<-factor(cla[match(s$COVTYPE,clan)])

r<-raster(nrow=2000,ncol=3000,ext=extent(s))

rx<-fasterize(s,r,field="cover",fun="max")
rn<-fasterize(s,r,field="cover",fun="min")

r<-ratify(rx)

levels(s$cover)
hn<-levels(r)[[1]][,1]
cols<-c("lightgoldenrod","skyblue","grey30","forestgreen","grey60","brown")
par(mar=c(0,0,0,0))
plot(r,legend=FALSE,col=cols)
legend("bottomright",legend=levels(s$cover),fill=cols,bty="n",border=NA,cex=2)
points(ds,col="red",pch=16,cex=0.8)


v<-velox(r)
b<-1000
g<-st_buffer(st_as_sf(ds),dist=b)
l<-v$extract(g)


l<-lapply(l,function(i){
  x<-i[,1]
  x<-x[!is.na(x)]
  sum(x==4L)/length(x)
})

ds$natural<-unlist(l)


grid<-st_make_grid(r,cellsize=c(500,500),what="centers")
pr<-raster(SpatialGrid((points2grid(as_Spatial(grid)))))
g<-st_buffer(st_as_sf(xyFromCell(pr,1:ncell(pr),spatial=TRUE)),dist=b)
l<-v$extract(g)

l<-lapply(l,function(i){
  x<-i[,1]
  x<-x[!is.na(x)]
  sum(x==4L)/length(x)
})
pr<-setValues(pr,unlist(l))
plot(pr,col=rev(viridis(100)))

#tmap_mode("view")
#tm_shape(ds) +
#  tm_dots() +
#  tm_layout(basemaps = c("Esri.WorldImagery", "Esri.WorldShadedRelief", "Esri.NatGeoWorldMap"))

plot(ds)
l<-locator()
h<-gConvexHull(SpatialPoints(cbind(l$x,l$y),proj4string=CRS(proj4string(ds))))
plot(h,add=TRUE)
o<-over(ds,h)
d<-d[!is.na(o),]
ds<-ds[!is.na(o),]

#plot(d$jul,log(d$sp+1))

#m<-gam(sp~s(jul),data=d,family=nb())

can<-raster::getData("GADM", country = "CAN", level = 2)
que<-can[can$NAME_1=="Québec",]

rm(rx,rn,s,can);gc();gc()

# "http://forobs.jrc.ec.europa.eu/products/glc2000/legend/GLC2000_Lccs_110604_export.xls"
#lcc<-as.data.frame(read_xls("C:/Users/rouf1703/Downloads/GLC2000_Lccs_110604_export.xls"))
#cov<-raster("C:/Users/rouf1703/Downloads/CAN_cov/CAN_cov.grd")
#cov<-crop(cov,extent(-76,-70,45,49))
#cov<-crop(cov,xs2)
#que<-spTransform(que,proj4string(ds))
#cov<-projectRaster(cov,crs=CRS(proj4string(ds)))
#mapview(cov)+mapview(ds[1:1000,])
#plot(ds)
#plot(que,add=TRUE)
#plot(cov,add=TRUE)



###########################################
#### TEMPORAL

x<-d[d$year==2015,]
x$sp<-x$A23
#$sp<-log(x$sp+1)
x$sp<-x$sp
l<-split(x,x$Site)
plot(x$jul,x$sp,col="white",xlim=c(120,300))
invisible(lapply(l,function(i){
  points(i$jul,i$sp,col=gray(0,0.2))
}))

#m<-gam(sp~s(jul),data=x,family=nb())
#m<-qgam(sp~s(jul),data=x,qu=0.1)
m<-gam(sp~s(jul),data=x,qu=0.1,family=gaussian(link=log))

v<-seq(min(x$jul),max(x$jul),by=0.1)
p<-predict(m,data.frame(jul=v),type="response",se.fit=TRUE)
lines(v,p$fit)
lines(v,p$fit+2*p$se.fit,lty=3)
lines(v,p$fit-2*p$se.fit,lty=3)




#########################################
#### SPATIAL

xs<-ds[ds$year==2014 & ds$Week==27,]
xs$sp<-log(xs$A29+1)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(xs,cex=0.05,pch=16)
plot(que,add=TRUE,border="grey50")
points(xs,pch=16,cex=xs$sp/3,col=gray(0,0.2))
#par(new=FALSE)
#plot(xs$jul,xs$sp)
mg<-glm(sp~1,data=xs@data)
coords <- coordinates(xs)
v<-variog(coords=coords,data=resid(mg),breaks=seq(0,100000,by=2000),max.dist=100000)
#v<-variog(coords=coords,data=resid(m))
plot(v,type="b") 
#par(mfrow=c(1,1))


prdomain <- inla.nonconvex.hull(coordinates(xs),convex=-0.05, resolution = c(100, 100))
mesh<-inla.mesh.2d(loc=coordinates(xs),max.edge=c(3000,40000),offset=c(20000,10000),cutoff=3000,boundary=prdomain)
plot(mesh,asp=1)

#spde<-inla.spde2.matern(mesh,alpha=2)
spde<-inla.spde2.pcmatern(mesh,prior.range=c(50000,0.9),prior.sigma=c(2,0.1),alpha=2)

g<-makegrid(xs,n=10000) # makes sure pixels touching are included too
g<-SpatialPoints(g,proj4string=CRS(proj4string(xs)))
g<-SpatialPixels(g)
o<-over(as(g,"SpatialPolygons"),gBuffer(xs,width=50000))
#g<-g[apply(o,1,function(i){!all(is.na(i))}),]
g<-g[!is.na(o),]

s.index<-inla.spde.make.index(name="spatial",n.spde=spde$n.spde)

model<-sp~-1+intercept+f(spatial,model=spde)

A<-inla.spde.make.A(mesh=mesh,loc=coordinates(xs))
Ap<-inla.spde.make.A(mesh=mesh,loc=coordinates(g))
n<-100

stack.est<-inla.stack(data=list(sp=xs$sp),A=list(A),effects=list(c(s.index,list(intercept=1))),tag="est")
stack.latent<-inla.stack(data=list(sp=NA),A=list(Ap),effects=list(s.index),tag="latent")

full.stack<-inla.stack(stack.est,stack.latent)

index.est<-inla.stack.index(full.stack,tag="est")$data
index.latent<-inla.stack.index(full.stack,tag="latent")$data

m<-inla(model,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),family="gaussian",control.compute=list(dic=TRUE,waic=TRUE))

### map predictions
#par(mfrow=c(1,1),oma=c(0,5,0,5))
plot(xs,type="n",axes=TRUE)
p<-m$summary.fitted.values[index.latent,"mean"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(raster(gp),col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=2,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),legend.args=list(text='logAbundance', side=1, font=2, line=2.3),horizontal=TRUE)
plot(que,add=TRUE,border=gray(0,0.05),lwd=0.01)
points(xs,pch=1,cex=xs$sp/3,col=gray(0,0.15))
#par(mfrow=c(1,1))

#########################################
### GAM

xs<-ds
xs$sp<-log(xs$A9+1)
xs$sp<-xs$A29
xs$lon<-coordinates(xs)[,1]
xs$lat<-coordinates(xs)[,2]
xs$year<-as.factor(xs$Year)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(xs,cex=0.05,pch=16)
plot(que,add=TRUE,border="grey50")
points(xs,pch=16,cex=xs$sp/3,col=gray(0,0.2))
#par(new=FALSE)
#plot(xs$jul,xs$sp)

m<-gam(sp~te(lon,lat,Week)+s(year,bs="re"),data=xs@data,family=nb())

g<-makegrid(xs,n=2000) # makes sure pixels touching are included too
g<-SpatialPoints(g,proj4string=CRS(proj4string(xs)))
g<-SpatialPixels(g)
o<-over(as(g,"SpatialPolygons"),gConvexHull(xs))
#g<-g[apply(o,1,function(i){!all(is.na(i))}),]
g<-g[!is.na(o),]

newdat<-data.frame(lon=coordinates(g)[,1],lat=coordinates(g)[,2],Week=30,year=2015)

p<-predict(m,newdat,type="response")
#p<-exp(p)-1

plot(xs,type="n",axes=TRUE)
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
#brks <- seq(min(p),max(p),by=0.01)
br<-classIntervals(p,n=10,style="kmeans",rtimes = 1)
brks<-br$brks
factor(brks)
cols<-colo.scale(length(brks)-1,rev(brewer.pal(11,"RdYlGn")))

r<-raster(gp)
r<- cut(r,breaks=brks)
r<-as.factor(r)
rat <- levels(r)[[1]]
rat[["ab"]] <- c("land","ocean/lake", "rivers","water bodies")
levels(r) <- rat


plot(raster(gp),col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=2,add=TRUE,axis.args=list(at=brks, labels=round(brks,0)),legend.args=list(text='logAbundance', side=1, font=2, line=2.3),horizontal=TRUE)
plot(que,add=TRUE,border=gray(0,0.25),lwd=0.01)
points(xs,pch=1,cex=xs$sp/3,col=gray(0,0.3))

#########################################
### Rasters

xs<-ds[ds$year==2014,]
xs$sp<-log(xs$A29+1)

#xs<-ds
#xs$sp<-log(xs$A23+1)
#xs$sp<-xs$A9

r <- raster(ncol = 20, nrow = 15, ext = extent(xs))
l<-split(xs,xs$Week)
lr<-lapply(l,function(i){
  rasterize(i,r,field="sp",fun=median,background=NA)
})
r<-stack(lr)
r<-log(r+1)
levelplot(r,col.regions=rev(rasterTheme()$regions$col))

#########################################
### RW1

xs<-ds[ds$year==2015,]
xs$sp<-log(xs$A9+1)

plot(xs$jul,xs$sp)

g<-gam(sp~s(jul),data=xs@data)

v<-120:300
newdat<-data.frame(jul=v)

p<-predict(g,newdat,type="link",se.fit=TRUE)
plot(xs$jul,xs$sp)
lines(v,p$fit)
lines(v,p$fit+2*p$se.fit,lty=3)
lines(v,p$fit-2*p$se.fit,lty=3)


dat<-rbind.fill(xs@data,newdat)
i<-inla(sp~f(jul,model="rw1"),data=dat,control.predictor=list(compute=TRUE))

p<-i$summary.fitted.values[(nrow(xs)+1):(nrow(i$summary.fitted.values)),c("mean","0.025quant","0.975quant")]

lines(v,p$mean,col="red")

######################
### weeks
library(reshape)
library(tidyr)
par(mar=c(0,0,0,3))
d$wy<-paste(d$year,d$Week,sep="_")
d$k<-1
x<-reshape(d[,c("Site","wy","k")], idvar="Site", timevar="wy", direction="wide")
x<-x[,c(names(x)[1],sort(names(x)[-1]))]
x<-raster(as.matrix(x)[,-1])
plot(x)

################################################
### space-time simple from spde tutorial

xs<-ds[ds$year=="2014",]
xs$sp<-log(xs$A29+1)
#xs$sp<-xs$A29
xs<-xs[order(xs$Annee,xs$Week),]

prdomain <- inla.nonconvex.hull(coordinates(xs),convex=-0.05, resolution = c(100, 100))
prmesh1<-inla.mesh.2d(loc=coordinates(xs),max.edge=c(4000,40000),offset=c(20000,10000),cutoff=4000,boundary=prdomain)
plot(prmesh1,asp=1)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
  mesh=prmesh1, alpha=2, ### mesh and smoothness parameter
  prior.range=c(50000, 0.9), ### P(practic.range<0.05)=0.01
  prior.sigma=c(4, 0.1)) ### P(sigma>1)=0.01

## ----rfindex-------------------------------------------------------------
k<-length(unique(xs$Week))
iset <- inla.spde.make.index('i', n.spde=spde$n.spde, n.group=k)

## ----apred---------------------------------------------------------------
A <- inla.spde.make.A(mesh=prmesh1, 
                      loc=coordinates(xs), 
                      group=as.integer(factor(xs$Week))) 

## ----stack---------------------------------------------------------------
sdat<-inla.stack(tag='stdata',data=list(y=xs$sp),A=list(A,1),effects=list(iset,w=xs$Typ)) 

## ----hbeta---------------------------------------------------------------
h.spec <- list(theta=list(prior='pccor1',param=c(0, 0.95)))

## ----remote,echo=FALSE---------------------------------------------------
##inla.setOption(inla.call='remote')

## ----ft------------------------------------------------------------------
formulae <- y ~ 0 + w + 
  f(i, model=spde, group=i.group, 
    control.group=list(model='ar1', hyper=h.spec)) 
prec.prior <- list(prior='pc.prec', param=c(4, 0.01))
m <- inla(formulae,  data=inla.stack.data(sdat), 
            control.predictor=list(compute=TRUE, A=inla.stack.A(sdat),link=1), 
            control.family=list(hyper=list(theta=prec.prior)), 
            control.fixed=list(expand.factor.strategy='inla'),
            family="gaussian")

## ----sbeta---------------------------------------------------------------
round(cbind(observed=tapply(dat$y, dat$w, mean), m$summary.fixed), 4) 

## ----echo=FALSE,fig.width=5.5,fig.height=5.5-----------------------------
par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=2:0)
for (j in 1:4) {
  plot(m$marginals.hyper[[j]], type='l', 
       xlab=names(m$marginals.hyper)[j], ylab='Density')
  abline(v=c(1/sd.y^2, sqrt(8)/params[1], 
             params[2]^0.5, rho)[j], col=2)
}

## ----rfidx---------------------------------------------------------------
str(idat <- inla.stack.index(sdat, 'stdata')$data) 

## ----meanrf--------------------------------------------------------------
cor(xs$sp, m$summary.linear.predictor$mean[idat])

## ----projgrid------------------------------------------------------------
coords<-coordinates(xs)

stepsize <- 400*1/1
nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
projgrid <- inla.mesh.projector(prmesh1, xlim=range(coords[,1]),ylim=range(coords[,2]), dims=nxy,crs=CRS(proj4string(xs)))

## ----projpmean-----------------------------------------------------------
xmean <- list()
for (j in 1:k){
  xmean[[j]] <- inla.mesh.project(projgrid,m$summary.random$i$mean[iset$i.group==j])
}

## ----inout---------------------------------------------------------------
b<-gBuffer(gConvexHull(SpatialPoints(prdomain$loc,p=CRS(proj4string(ds)))),width=100,byid=FALSE)
a<-concaveman(coordinates(xs),2)
a<-gBuffer(spPolygons(a,crs=CRS(proj4string(xs))),width=5000)
#plot(a)
#plot(xs,add=TRUE)
o <- over(SpatialPoints(projgrid$lattice$loc,p=CRS(proj4string(ds))),b)

## ----setNAs---------------------------------------------------------------
for (j in 1:k)   xmean[[j]][is.na(o)] <- NA


r<-stack(lapply(xmean,function(i){
  raster(nrows=nxy[2], ncols=nxy[1], xmn=min(projgrid$x), xmx=max(projgrid$x), ymn=min(projgrid$y), ymx=max(projgrid$y),crs=CRS(proj4string(xs)),vals=as.vector(i[,ncol(i):1])) ## some crazy ordering in INLA output be careful
  #raster(i)
}))
names(r)<-unique(xs$Week)
#r<-exp(r)-1
cols<-colo.scale(200,rev(brewer.pal(11,"RdYlBu")))

# voir argument panel.number ou packets de layer

xxs<-split(xs,xs$Week)
levelplot(r,col.regions=cols,cuts=199)+
  layer(sp.points(xxs[[panel.number()]],col=gray(0,0.1),pch=1,cex=identity(xxs[[panel.number()]]$sp)/3))
  #+#layer(sp.polygons(xxs))

#plot(m$summary.fitted.values$mean[1:nrow(xs)],xs$sp)

#levelplot(exp(r)-1,zscaleLog=10)

####################################
### test color raster
####################################
  
n<-10000
r<-raster(matrix(runif(n,0,1000),ncol=sqrt(n)))

cols<-colo.scale(200,rev(brewer.pal(11,"RdYlBu")))
levelplot(r,col.regions=cols,cuts=199)


breaks <- seq(0, 1000, by=10)
b <- c(0,10,20,50,100,500,1000)

cols <- colorRampPalette(c("red", "yellow", "lightgreen"))(length(breaks)-1)
#cols<-colo.scale(200,rev(brewer.pal(11,"RdYlBu")))
levelplot(r,col.regions=cols)
levelplot(r)



####################################
### nbinom
####################################

x<-rnbinom(n=100000,size=0.01,mu=10)
mean(x)
hist(x)

####################################
### gam
####################################

xs<-ds[ds$year=="2005" & ds$Week==29,]
xs$x<-coordinates(xs)[,1]
xs$y<-coordinates(xs)[,2]
#xs$sp<-log(xs$A9+1)
xs$sp<-xs$A23
xs<-xs[order(xs$Annee,xs$Week),]

points(xs,cex=log(xs$sp+1)/2)

#registerDoParallel(6)
#getDoParWorkers() 

ran<-seq(1000,50000,by=1000)
REML<-foreach(i=seq_along(ran),.packages="mgcv") %dopar% {
  fit<-gam(sp~natural+s(x,y,bs="gp",k=40,m=c(2,ran[i])),data=xs@data,family=tw,method="REML")
  fit$gcv.ubre
}
plot(ran,unlist(REML),type="l")
rho<-ran[which.min(unlist(REML))]
#rho<-9000

fit<-gam(sp~natural+s(x,y,bs="gp",k=40,m=c(2,rho)),data=xs@data,family=tw,method="REML")
#ctrl<-list(nthreads=6)
#fit<-gam(sp~1+s(x,y,bs="ds",k=40,m=c(1,0.5)),data=xs@data,family=tw,method="REML")

pred<-crop(pr,extent(gBuffer(bbox2pol(ds),width=2000)))
ce<-xyFromCell(pred,1:ncell(pred),spatial=TRUE)
proj4string(ce)<-proj4string(xs)
p<-predict(fit,newdata=data.frame(Week=27,x=coordinates(ce)[,1],y=coordinates(ce)[,2],natural=values(pred),year=2014),type="response",se.fit=TRUE)$fit
pred<-setValues(pred,as.vector(p))
plot(pred,col=rev(viridis_pal(option="viridis")(100)))
points(xs,cex=log(xs$sp+1),col=gray(0,0.15),xpd=TRUE)
points(ds[!duplicated(ds$Site),],cex=0.35,pch=3,col=gray(0,0.15))
#text(xs,cex=0.6,label=xs$sp,col="tomato4")



####################################
### glmmTMB
####################################

xs<-ds[ds$year=="2016" & ds$Week==30,]
xs$x<-coordinates(xs)[,1]
xs$y<-coordinates(xs)[,2]
xs$sp<-log(xs$A29+1)#
#xs$sp<-xs$A29
xs<-xs[order(xs$Annee,xs$Week),]

points(xs,cex=log(xs$sp+1)/2)

xs$pos <- numFactor(coordinates(xs)[,1:2])
xs$group <- factor(rep(1, nrow(xs)))

fit<-glmmTMB(sp~1+exp(pos+0|group),data=xs@data)
visreg(fit,scale="response")


pred<-crop(pr,extent(gBuffer(bbox2pol(xs),width=10000)))
pred<-aggregate(pred,5)
ce<-xyFromCell(pred,1:ncell(pred),spatial=TRUE)
proj4string(ce)<-proj4string(xs)
p<-predict(fit,newdata=data.frame(pos=numFactor(coordinates(ce)[,1:2]),natural=values(pred),group=factor(rep(1,ncell(pred)))),type="response", allow.new.levels=TRUE)
pred<-setValues(pred,as.vector(p))
plot(pred,col=rev(viridis_pal(option="viridis")(100)))
points(xs,cex=log(xs$sp+1))
text(xs,cex=0.6,label=xs$sp,col="tomato4")


####################################
### climate
####################################

t<-raster::getData('worldclim', var='tmin', res=0.5, lon=-70, lat=45)
xs2<-spTransform(xs,CRS=CRS(proj4string(t)))
t<-crop(t,xs2)
plot(subset(t,1))
plot(xs2,add=TRUE)



#tmap_mode("view")
#tm_shape(cov) +
#  tm_raster(palette = rev(terrain.colors(100)),n=12) +
#  tm_layout(basemaps = c("Esri.WorldImagery", "Esri.WorldShadedRelief", "Esri.NatGeoWorldMap"))


#######################################
### species counts
#######################################

x<-as.data.table(ds@data)
x<-x[,lapply(.SD,sum),by=Annee,.SDcols=paste0("A",1:29)][order(Annee)]
x


#######################################
### se.fit pattern
#######################################

n<-100
x<-c(rnorm(n,0.5,0.1),rnorm(n,-0.5,0.1))
y<-c(rnorm(n,0.5,0.1),rnorm(n,0.5,0.1))
sp<-x+rnorm(n+n,0,0.2)+10
fit<-gam(sp~s(x,y,k=50,bs="ds",m=c(1,0.5)),data=data.frame(x,y,sp))
#fit<-gam(sp~s(x,y,k=50),data=data.frame(x,y,sp))
r<-raster(nrow=200,ncol=200,ext=extent(-1,1,0,1))
ce<-xyFromCell(r,1:ncell(r),spatial=TRUE)
p<-predict(fit,newdata=data.frame(x=coordinates(ce)[,1],y=coordinates(ce)[,2]),se.fit=TRUE)$se.fit
#plot(ds$Week,ds$sp)
r<-setValues(r,as.vector(p))
#r[is.na(over(ce,h))]<-NA
plot(r,col=rev(viridis_pal(option="viridis")(100)))
points(x,y,cex=10*(max(sp)-sp)/max(sp))



#######################################
### se.fit pattern2
#######################################

n<-500
x<-c(rnorm(n,0.5,0.1))
y<-c(rnorm(n,0.5,0.1))
sp<-2*x+rnorm(n,0,0.2)+10

m1<-gam(sp~s(x,y,k=50),data=data.frame(x,y,sp))
m2<-gam(sp~s(x,y,k=50,bs="ds",m=c(1,0.5)),data=data.frame(x,y,sp))



r<-raster(nrow=200,ncol=200,ext=extent(cbind(x,y)))
ce<-xyFromCell(r,1:ncell(r),spatial=TRUE)

par(mfrow=c(1,2))

p1<-predict(m1,newdata=data.frame(x=coordinates(ce)[,1],y=coordinates(ce)[,2]),se.fit=TRUE)$se.fit
r1<-setValues(r,as.vector(p1))
plot(r1,col=rev(viridis_pal(option="viridis")(100)))
points(x,y,cex=10*(max(sp)-sp)/max(sp))

p2<-predict(m2,newdata=data.frame(x=coordinates(ce)[,1],y=coordinates(ce)[,2]),se.fit=TRUE)$se.fit
r2<-setValues(r,as.vector(p2))
plot(r2,col=rev(viridis_pal(option="viridis")(100)))
points(x,y,cex=10*(max(sp)-sp)/max(sp))















