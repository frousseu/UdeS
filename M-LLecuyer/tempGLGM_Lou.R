
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
library(RColorBrewer)
library(tmap)
library(AICcmodavg)

# FRutils doit être installé une fois comme cela pour obtenir la fonction colo.scale:

#library(devtools)
#install_github("frousseu/FRutils")

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

#d<-as.data.frame(fread('C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/LandEco_Complet_22_05.txt',dec=",",sep="\t"))

d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/LandEco_Complet_22_05.xlsx"),stringsAsFactors=FALSE)

#d[]<-lapply(d,function(i){
#  if(any(grep(",",i))){
#    as.numeric(gsub(",",".",i))
#  }else{
#    i
#  }
#})

### test bt adding data
#d<-rbind(d,d,d,d,d)
#d$X<-d$X+rnorm(nrow(d),0,0.05)
#d$Y<-d$Y+rnorm(nrow(d),0,0.05)

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

plot(ds)


#########################################
### import raster and build pred grid
#########################################

### import land cover raster
path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Landcover_2015_extended.tif"
code<-read.table("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Legend.txt",header=TRUE,sep=",",stringsAsFactors=FALSE)
code<-code[code$Code!=0,]# ignore background values in raster
r <- stack(path)

plot(r)
plot(ds,add=TRUE)

### import 


### build grid to extract pixel values around a buffer where predictions are going to be made
ras <- raster(xmn=bbox(r)[1,1],xmx=bbox(r)[1,2],ymn=bbox(r)[2,1],ymx=bbox(r)[2,2],ncols=100,nrows=100)
ras[] <- runif(ncell(ras))
g<-as(ras,"SpatialPixelsDataFrame")
proj4string(g)<-proj4string(r)
w<-500 # this is the buffer width, values over 2000 start to take a long time to cumpute
p<-gBuffer(SpatialPoints(coordinates(g)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)

plot(p,add=TRUE)

### extract values from grid using velox cause it's faster
v<-velox(r)
e<-v$extract(p,fun=function(x){paste(x,collapse="_")})

### compute proportions for different categories in at each location
ans<-lapply(e,function(i){
  x<-table(unlist(strsplit(i,"_")))
  x<-x[names(x)!=0] # ignore background values in raster
  s<-setdiff(code$Code,names(x))
  #s<-setdiff(names(x),names(x))
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

ds@data<-cbind(ds@data,obs)

### add pop index

pop<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="poblacion")


di<-gDistance(pop,ds,byid=TRUE)
f<-sqrt
d$index_pop2<-rowSums(t(apply(di,1,function(i){exp(((log(0.1)/2000))*i)*f(pop$POBTOT)})))
d$index_pop4<-rowSums(t(apply(di,1,function(i){exp(((log(0.1)/4000))*i)*f(pop$POBTOT)})))
d$index_pop8<-rowSums(t(apply(di,1,function(i){exp(((log(0.1)/8000))*i)*f(pop$POBTOT)})))
#d$index_pop16<-rowSums(t(apply(g,1,function(i){exp(((log(0.1)/16000))*i)*f(s$z)})))
#d$index_pop32<-rowSums(t(apply(g,1,function(i){exp(((log(0.1)/32000))*i)*f(s$z)})))

m1<-glm(Attack~Cat_Typ,data=d,family="binomial")
m2<-glm(Attack~Cat_Typ+index_pop2,data=d,family="binomial")
m3<-glm(Attack~Cat_Typ+index_pop4,data=d,family="binomial")
m4<-glm(Attack~Cat_Typ+index_pop8,data=d,family="binomial")
#m5<-glm(Attack~Cat_Typ+index_pop16,data=d,family="binomial")
#m6<-glm(Attack~Cat_Typ+index_pop32,data=d,family="binomial")

ml<-list(m1=m1,m2=m2,m3=m3,m4=m4)#,m5=m5,m6=m6)

aictab(ml)

visreg(m4,"index_pop8",scale="response")

par(mfrow=c(1,2),oma=c(1,1,0,3))
x<-0:20000
plot(x,exp(co*x),type="l",ylim=0:1,xlim=c(0,10000),xlab="distance",ylab="pondération")
abline(0.1,0,lty=2)
plot(r)

### add preidctions

ds2<-ds[1:10,]
ds2$Cat_Typ<-"C"
ds2$Secondary<-seq(0,100,length.out=10)
ds2$Attack<-NA

ds<-rbind(ds,ds2)



######################################################################
### build model formulas
######################################################################

m0<-Attack~1
m1<-Attack~Cat_Typ
m2<-Attack~Cat_Typ+Secondary

ml<-list(m0=m0,m1=m1,m2=m2)

# Puisque Dist_Road est manquantr du raster, on va surtout utiliser le m3

######################################################################
### Explore variogram from glm residuals from model3
######################################################################

glm1<-glm(m2,data=d,family="binomial")

coords<-as.matrix(d[,c("X","Y")])
v<-variog(coords=coords,data=resid(glm1),breaks=seq(0,50000,by=500))
fitv<-variofit(v,ini.cov.pars=c(2,5000),cov.model="matern",fix.nugget=FALSE,nugget=0,fix.kappa=TRUE,kappa=1)
plot(v)
lines(fitv)


######################################################################
### Build prediction raster stack with available variables 
######################################################################

covr<-stack(g[,c("Secondary","Pasture")])
glm1<-glm(m2,data=d,family="binomial")
p<-predict(glm1,data.frame(as.matrix(aggregate(covr,fac=2)),Cat_Typ="C"),type="response")

# this is the predicted attack grid from the simple glm with model4
tempr<-aggregate(covr,fac=2)[[1]]
tempr[]<-p

plot(covr)

#######################################################################
### evalute models and store them in a list
#######################################################################

# Only variable available in raster format (Secondary, Pasture) are used in predictions

predr<-aggregate(covr,fac=2) ## this produce a coarser grid to make computations faster

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

print(ml)

mm<-lapply(ml,function(i){
  glgm(i, 
       data=ds,
       grid=20,
       covariates=NULL, 
       family="binomial", 
       buffer=10000,
       shape=1,
       priorCI=list(sd=c(0.2,4),range=c(2000,50000)),
       control.compute=list(waic=TRUE,mlik=TRUE))#,
       #control.predictor=list(link=c(rep(NA,101),rep(1,10))))
})
names(mm)<-names(ml)

# summary(mm[["m3"]]$inla)

##############################################################
### check priors and posteriors for matern parameters
##############################################################

m<-mm[["m2"]]

par(mfrow=c(1,2))

ma<-max(c(m$parameters$sd$prior[,2],m$parameters$sd$posterior[,2]))
plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density',lty=2,xlim=c(0,5),ylim=c(0,ma))
lines(m$parameters$sd$posterior,lty=1)
#lines(seq(0,10,by=0.1),dlgamma(seq(0,10,by=0.1),m$parameters$sd$params.intern$param[1],m$parameters$sd$params.intern$param[2]),col="blue")
legend("topright", lty=2:1, legend=c("prior","posterior"))

ma<-max(c(m$parameters$range$prior[,2],m$parameters$range$posterior[,2]))
plot(m$parameters$range$posterior,type="l",xlim = c(0,50*1000),xlab='range (m)', ylab='density',lty=1,ylim=c(0,ma))
lines(m$parameters$range$prior,lty=2)
#lines(dgamma(seq(0,50000,by=10),m$parameters$range$params.intern[1],m$parameters$range$params.intern[2]),col="red")
legend("topright", lty=2:1, legend=c("prior","posterior"))


##############################################################
### plot prediction graphs
##############################################################

m$inla$summary.fitted.values


##############################################################
### plot prediction rasters
##############################################################

brks<-seq(0,1,by=0.01)
lab.brks<-ifelse(((brks*100)%%10)==0,brks,"")
cols<-rev(gray((0:100)/100))
par(mfrow=c(1,2))
plot(m$raster$predict.invlogit,col=cols,breaks=brks,lab.breaks=lab.brks)
plot(tempr,col=cols,breaks=brks,lab.breaks=lab.brks)

#plot(stack(m$raster$predict.invlogit,tempr))


### exploratory
mymatern<-function(u,phi,kappa){
  uscale<-sqrt(8*kappa)*(u/phi)
  res<-(1/(gamma(kappa)*2^(kappa-1)))*(uscale^kappa)*besselK(uscale,kappa)
  res[u==0]=1
  res
}
dist<-seq(0,50000,10)
maternfit<-m$param$summary["sd",c("0.5quant")]*(1-mymatern(dist,phi=m$param$summary["range",c("0.5quant")],kappa=1))
maternfit<-mymatern(dist,phi=m$param$summary["range",c("0.5quant")],kappa=1)
plot(maternfit~dist,type="l")
lines(maternfit~dist,col="red")




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
     priorCI=list(sd=c(0.2,10),range=c(2000,50000))
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
tm_shape(mpred$raster[["predict.invlogit"]])+tm_raster(alpha=0.6,palette=colo.scale(1:10,bp),n=7)+tm_shape(ds)+tm_dots("Attack")






###############################################################
### extract matern parameters (work in progress...)
###############################################################

### extract matern parameters and plot variogram function
param=c(shape=1,range=20000,variance=0.06^2)
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


# example with raster
myraster = raster(nrows=40,ncols=60,xmn=150000,xmx=300000,ymn=1950000,ymx=2150000)
param = c(range=200000, shape=2,	anisoRatio=0, 
          anisoAngleDegrees=0,variance=10)
# plot correlation of each cell with the origin
myMatern = matern(myraster, y=c(225000,2050000), param=param)
plot(myMatern, main="anisortopic matern")


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



