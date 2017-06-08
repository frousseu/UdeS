
library(splancs)
library(geoR)
library(gstat)
library(sp)
library(scales)
library(ncf)
library(car)
library(data.table)
library(readxl)
library(raster)
library(visreg)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(INLA)
library(velox)
library(geostatsp)
library(FRutils)

# FRutils doit être installé une fois comme cela pour obtenir la fonction colo.scale:

library(devtools)
install_github("frousseu/FRutils")

####################################################################
### function to build waic table from a list of geostatsp models
####################################################################

# x is a named list of geostatsp models 

waictab<-function(x){
  stopifnot(!is.null(names(x)))
  waic<-sapply(x,function(i){
    i$inla$waic$waic  
  })
  p.eff<-sapply(x,function(i){
    i$inla$waic$p.eff
  })
  dwaic<-waic-min(waic)
  w<-sapply(dwaic,function(i){
    exp((-1/2)*i)/sum(exp((-1/2)*dwaic))
  })
  d<-data.frame(model=names(x),p.eff=round(p.eff,2),waic=round(waic,2),dwaic=round(dwaic,2),w=round(w,2))
  d<-d[order(d$dwaic),]
  d
}

################################################
### load data
################################################

d<-as.data.frame(fread('C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/LandEco_Complet_22_05.txt',dec=",",sep="\t"))
d[]<-lapply(d,function(i){
  if(any(grep(",",i))){
    as.numeric(gsub(",",".",i))
  }else{
    i
  }
})

### build shapefile of locations
ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
ds<-d

coordinates(ds)<- ~X+Y
proj4string(ds) <- ll
ds<-spTransform(ds,CRS(prj))

d$X<-coordinates(ds)[,1]
d$Y<-coordinates(ds)[,2]

d$x<-d$X
d$y<-d$Y

d$Cat_Typ<-as.factor(d$Cat_Typ)

# ds is a spatial object that can be plotted


#########################################
### import raster and build pred grid
#########################################

#path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Mature_forest_2015_include_bajos_secondary.tif"
#path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Mature_forest_2015_not_include_bajos_secondary.tif"
path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Landcover_2015_extended.tif"
code<-read.table("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Legend.txt",header=TRUE,sep=",",stringsAsFactors=FALSE)
code<-code[code$Code!=0,]# ignore background values in raster
r <- stack(path)

### build grid to extract pixel values around a buffer where predictions are going to be made
ras <- raster(xmn=bbox(r)[1,1],xmx=bbox(r)[1,2],ymn=bbox(r)[2,1],ymx=bbox(r)[2,2],ncols=100,nrows=100)
ras[] <- runif(ncell(ras))
g<-as(ras,"SpatialPixelsDataFrame")
proj4string(g)<-proj4string(r)
w<-500 # this is the buffer width, values over 2000 start to take a long time to cumpute
p<-gBuffer(SpatialPoints(coordinates(g)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
plot(g)
plot(p,add=TRUE)

### extract values from grid using velox cause it's faster
v<-velox(r)
e<-v$extract(p,fun=function(x){paste(x,collapse="_")})

### compute proportions for different categories in at each location
ans<-lapply(e,function(i){
  x<-table(unlist(strsplit(i,"_")))
  x<-x[names(x)!=0] # ignore background values in raster
  s<-setdiff(code$Code,names(x))
  if(length(s)==0){
    x<-x[order(names(x))]
  }else{
    temp<-rep(0,length(s))
    names(temp)<-s
    x<-c(x,temp)
    x<-x[order(names(x))]
  }
  100*x/sum(x)
})
land<-as.data.frame(do.call("rbind",ans))
names(land)<-gsub(" ","_",code$Category[match(names(land),code$Code)])
g@data<-land
g<-g[!is.na(g@data[,1]),] #remove NA values for grid without only background data


### get proportions at observations
p<-gBuffer(SpatialPoints(coordinates(ds)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
eobs<-v$extract(p,fun=function(x){paste(x,collapse="_")})
ans<-lapply(eobs,function(i){
  x<-table(unlist(strsplit(i,"_")))
  x<-x[names(x)!=0] # ignore background values in raster
  s<-setdiff(code$Code,names(x))
  if(length(s)==0){
    x<-x[order(names(x))]
  }else{
    temp<-rep(0,length(s))
    names(temp)<-s
    x<-c(x,temp)
    x<-x[order(names(x))]
  }
  100*x/sum(x)
})
obs<-as.data.frame(do.call("rbind",ans))
names(obs)<-gsub(" ","_",code$Category[match(names(obs),code$Code)])
d<-cbind(d,obs,stringsAsFactors=FALSE)
names(d)[which(names(d)=="Agriculture\t")]<-"Agriculture"
names(g)[which(names(g)=="Agriculture\t")]<-"Agriculture"
g$Cat_Typ<-"C"

#g is a SpatialPixelsDataFrame

######################################################################
### build model formulas
######################################################################

m0<-Attack~1
m1<-Attack~Cat_Typ
m2<-Attack~Cat_Typ+Secondary
m3<-Attack~Cat_Typ+Secondary+Pasture 
m4<-Attack~Cat_Typ+Secondary+Dist_Road

ml<-list(m0=m0,m1=m1,m2=m2,m3=m3,m4=m4)

# Puisque Dist_Road est manquantr du raster, on va surtout utiliser le m3

######################################################################
### Explore variogram from glm residuals from model3
######################################################################

glm1<-glm(m3,data=d,family="binomial")

coords<-as.matrix(d[,c("X","Y")])
v<-variog(coords=coords,data=resid(glm1),breaks=seq(0,50000,by=500))
fitv<-variofit(v,ini.cov.pars=c(2,5000),cov.model="matern",fix.nugget=FALSE,nugget=0,fix.kappa=TRUE,kappa=1)
plot(v)
lines(fitv)


######################################################################
### Build prediction raster stack with available variables 
######################################################################

covr<-stack(g[,c("Secondary","Pasture")])
glm1<-glm(m3,data=d,family="binomial")
p<-predict(glm1,data.frame(as.matrix(covr),Cat_Typ="C"),type="response")

# this is the predicted attack grid from the simple glm with model4
tempr<-covr[[1]]
tempr[]<-p

plot(covr)

#######################################################################
### evalute models and store them in a list
#######################################################################

# Only variable available in raster format (Secondary, Pasture) are used in predictions

predr<-aggregate(covr,fac=3) ## this produce a coarser grid to make computations faster

plot(predr)

### ceci peut être long à rouler...
#mm<-lapply(ml,function(i){
#                glgm(i, 
#                  data=ds,
#                  grid=predr,
#                  covariates=covr, 
#                  family="binomial", 
#                  buffer=10000,
#                  shape=1,
#                  priorCI=list(sd=c(u=0.2, alpha=0.05),range=c(2000,50000)),
#                  control.compute=list(waic=TRUE,mlik=TRUE))
#})
#names(mm)<-names(ml)

### version sans predictions précises dans le raster, roule plus vite
mm<-lapply(ml,function(i){
  glgm(i, 
       data=ds,
       grid=20,
       covariates=NULL, 
       family="binomial", 
       buffer=10000,
       shape=1,
       priorCI=list(sd=c(u=0.2, alpha=0.05),range=c(2000,50000)),
       control.compute=list(waic=TRUE,mlik=TRUE))
})
names(mm)<-names(ml)

summary(mm[["m4"]]$inla)

##############################################################
### WAIC table and possible model-averaging with INLABMA
##############################################################

### liste des modèles soumise à waictab
l<-lapply(mm,function(i){i$inla})
waictab(mm)


#########################################################################
### Plot predictions
#########################################################################

### first run a model with raster predictions (m3)
mpred<-glgm(m3, 
     data=ds,
     grid=predr,
     covariates=covr, 
     family="binomial", 
     buffer=10000,
     shape=1,
     priorCI=list(sd=c(u=0.2, alpha=0.05),range=c(2000,50000))
)

### visualize predictions (static)
bp<-brewer.pal(5,"RdBu")
bp[3]<-"#4dac26"
cols<-colo.scale(1:100,bp)
par(mfrow=c(1,2))
plot(tempr,col=cols,main="glm")
plot(ds,add=TRUE,pch=1)
plot(mpred$raster[["predict.invlogit"]],col=cols,main="glgm")
plot(ds,add=TRUE,pch=1)

### visualisation prediction (dynamic)
tmap_mode("view")
tm_shape(mpred$raster[["predict.invlogit"]])+tm_raster(alpha=0.6,palette=colo.scale(1:10,bp),n=7)+tm_shape(ds)+tm_dots(col=c("Attack"))








############################################################
############################################################

### Below is exploratory

##############################################################
### check priors and posteriors for matern parameters
##############################################################

plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density', lwd=2,lty=2)
lines(m$parameters$sd$posterior, lty=1, lwd=2)
legend("topright", lty=2:1, legend=c("prior","posterior"))

plot(m$parameters$range$posterior,type="l",xlim = c(0,50*1000),xlab='range (m)', ylab='density', lwd=2,lty=1)
lines(m$parameters$range$prior, lty=2, lwd=2)
legend("topright", lty=2:1, legend=c("prior","posterior"))


###############################################################
### extract matern parameters (work in progress...)
###############################################################

### extract matern parameters and plot variogram function
param=c(shape=1,range=7194,variance=0.06^2)
u=seq(0,4,len=20)
uscale = sqrt(8*param['shape'])* u / param['range']
theMaterns = cbind(
  dist=u, 
  manual=	param['variance']*
    ( 1/(gamma(param['shape']) * 2^(param['shape']-1)  ) ) * 
    uscale^param['shape'] * besselK( uscale , param['shape']),
  geostatsp=matern(u, param=param)
)
head(theMaterns)
matplot(theMaterns[,'dist'], 
        theMaterns[,c('manual','geostatsp')],
        col=c('red','blue'), type='l')
legend('topright', fill=c('red','blue'),
       legend=c('manual','geostatsp'))



mymatern = function(u, phi, kappa) {
  uscale = sqrt(8 * kappa) * u/phi
  res = (1/(gamma(kappa) * 2^(kappa - 1))) * uscale^kappa *
    besselK(uscale, kappa)
  res[u == 0] = 1
  res
}

i<-1.1
u<-seq(0,15000,by=1)
lines(u,i*(1-mymatern(u,7194,1)),type="l",ylim=c(0,i))

v.f <- function(x, ...){3-cov.spatial(x, ...)}
### from geoR
curve(v.f(x, cov.pars = c(3, 10000),cov.model="matern", kappa = 1),from = 0, to = 15000,
      xlab = "distance", ylab = expression(gamma(h)), lty = 2,
      main = "models with equivalent \"practical\" range")
curve(v.f(x, cov.pars = c(1, 0.188), kappa = 1),from = 0, to = 1,
      add = TRUE)      
curve(v.f(x, cov.pars = c(1, 0.14), kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2, lty=2)      
curve(v.f(x, cov.pars = c(1, 0.117), kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2)      
legend("bottomright",
       expression(list(kappa == 0.5, phi == 0.250), 
                  list(kappa == 1, phi == 0.188), list(kappa == 2, phi == 0.140),
                  list(kappa == 3, phi == 0.117)), lty=c(2,1,2,1), lwd=c(1,1,2,2))
# plotting a nested variogram model



