
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
#library(brinla)
library(visreg)

load("~/UdeS/Consultation/GStetcher/Doc/LLF_occur.RData")

occ<-na.omit(llf.occur)
occ<-occ[sample(1:nrow(occ),5000),] # sample location to reduce computing time

occ$high_name<-as.factor(occ$high_name)
occ$logPop_2017<-log(occ$Pop_2017+0.2)

occs<-occ
coordinates(occs)<-~Longitude+Latitude
proj4string(occs)<-"+init=epsg:4326"

prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
occs<-spTransform(occs,CRS(prj))

#plot(occs,col=alpha(ifelse(occs$PA==1,"red","blue"),0.25),pch=16)

m1 <- glm (PA ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name, family = binomial(link = "logit"), data = na.omit(occ))
m1 <- glm (PA ~ -1+VEGZONSNA+logPop_2017+trees_age+Road_dens+WtrUrb_km, family = binomial(link = "logit"), data = na.omit(occ))
par(mfrow=c(3,3),mar=c(4,4,3,3))
visreg(m1,scale="response")
par(mfrow=c(1,1))

coords <- coordinates(occs)
v<-variog(coords=coords,data=resid(m1),breaks=seq(0,100000,by=500),max.dist=100000,bin.cloud=TRUE)
#fitv<-variofit(v,ini.cov.pars=c(0.2,30000),cov.model="exponential",fix.nugget=FALSE,nugget=0.125,fix.kappa=FALSE,kappa=0.45)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b") 
#lines(fitv)

swe <- raster::getData("GADM", country = "SWE", level = 0)
swe<-spTransform(swe,proj4string(occs))

#prdomain <- inla.nonconvex.hull(coordinates(occs),convex=-0.02, resolution = c(100, 100))
mesh<-inla.mesh.2d(loc=coordinates(occs),max.edge=c(20000,200000),offset=c(20000,50000),cutoff=20000,boundary=swe)
plot(mesh,asp=1)

#spde<-inla.spde2.matern(mesh,alpha=2)
spde<-inla.spde2.pcmatern(mesh,prior.range=c(100000,0.9),prior.sigma=c(3,0.1))

#m<-inla(model,data=list(y=occ$PA,intercept=rep(1,spde$n.spde),spatial=1:spde$n.spde),control.predictor=list(A=A,compute=TRUE),family="binomial")

g<-makegrid(swe,n=10000) # makes sure pixels touching are included too
g<-SpatialPoints(g,proj4string=CRS(proj4string(occs)))
g<-SpatialPixels(g)
o<-over(as(g,"SpatialPolygons"),swe)
g<-g[apply(o,1,function(i){!all(is.na(i))}),]

s.index<-inla.spde.make.index(name="spatial",n.spde=spde$n.spde)

model<-PA~-1+VEGZONSNA+intercept+logPop_2017+trees_age+Road_dens+WtrUrb_km+f(spatial,model=spde)

A<-inla.spde.make.A(mesh=mesh,loc=coordinates(occs))
Ap<-inla.spde.make.A(mesh=mesh,loc=coordinates(g))
n<-100
Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(312180,6342453),ncol=2)[rep(1,n),,drop=FALSE])

v<-setdiff(all.vars(model),c("PA","intercept","spatial","spde"))
lp<-newdata(x=occ[,v,drop=FALSE],v=v,n=n,fun=median,list=TRUE)
lpmed<-lapply(newdata(x=occ[,v,drop=FALSE],v=v,n=1,fun=median,list=TRUE)[[1]],function(i){rep(i,length(g))})


stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(occ[,v,drop=FALSE])),tag="est")
#stack.latent<-inla.stack(data=list(xi=NA),A=list(Ap),effects=list(s.index),tag="latent")
stack.map<-inla.stack(data=list(PA=NA),A=list(Ap,1),effects=list(c(s.index,list(intercept=1)),lpmed),tag="map")
#stack.map<-inla.stack(data=list(xi=NA),A=list(Ap),effects=list(s.index),tag="map")

full.stack<-inla.stack(stack.est,stack.map)

for(i in seq_along(v)){
  le<-length(lp[[v[i]]][[1]])
  if(le!=n){
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(819006,6545844),ncol=2)[rep(1,le),,drop=FALSE])
  }else{
    AA<-Apn
  }
  stack<-inla.stack(data=list(PA=NA),A=list(AA,1),effects=list(c(s.index,list(intercept=1)),lp[[v[i]]]),tag=v[i])     
  full.stack<-inla.stack(full.stack,stack)
}

index.est<-inla.stack.index(full.stack,tag="est")$data
index.map<-inla.stack.index(full.stack,tag="map")$data
index<-list(est=index.est,map=index.map)
for(i in seq_along(v)){
  index<-c(index,list(inla.stack.index(full.stack,tag=v[i])$data))
}  
names(index)[3:length(index)]<-v

m<-inla(model,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),family="binomial",control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=FALSE))


###########################################
### raw predictions on probability of usage
par(mfrow=c(1,4),oma=c(0,5,0,5))
plot(swe,border=gray(0,0.25),lwd=0.01)
points(occs,col=alpha(ifelse(occ$PA==1,"red","blue"),0.4),pch=ifelse(occ$PA==1,16,16),cex=ifelse(occ$PA==1,0.15,0.15))
p<-m$summary.fitted.values[index[["map"]],"mean"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
#cols<-colo.scale(length(brks)-1,rev(brewer.pal(11,"RdYlGn")))
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
#plot(gp,breaks=brks,col=cols,at=pretty(brks,10))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Probability of use', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)
###############################
### prediction on relative risk
p<-m$summary.fitted.values[index[["map"]],"mean"]
p<-p/(sum(occ$PA)/nrow(occ))
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
#cols<-colo.scale(length(brks)-1,rev(brewer.pal(11,"RdYlGn")))
cols<-colo.scale(brks-1,rev(brewer.pal(11,"RdBu")))
#plot(gp,breaks=brks,col=cols,at=pretty(brks,10))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Relative risk', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)
#par(mfrow=c(1,1))
##################################
### sd
p<-m$summary.fitted.values[index[["map"]],"sd"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Sd of probability of use', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)


#####################################
### graphs predictions
par(mfrow=c(2,3),mar=c(4,4,3,3),oma=c(0,10,0,0))
for(i in seq_along(v)){
  p<-m$summary.fitted.values[index[[v[i]]],c("0.025quant","0.5quant","0.975quant")]
  plot(lp[[v[i]]][[1]],p[,2],type="l",ylim=c(0,1),xlab=v[i],font=2,ylab="",lty=1)
  if(nrow(p)==n){
    lines(lp[[v[i]]][[1]],p[,1],lty=3)
    lines(lp[[v[i]]][[1]],p[,3],lty=3)
  }else{
    segments(x0=as.integer(lp[[v[i]]][[1]]),x1=as.integer(lp[[v[i]]][[1]]),y0=p[,1],y1=p[,3],lty=3)
  }
  points(occ[,v[i]],jitter(occ$PA,amount=0.035),pch=16,col=gray(0,0.15))
}
mtext("Probability of being an actual fire",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)





#image(inla.mesh.project(mesh,field=m$summary.fitted.values[inla.stack.index(full.stack,tag="latent")$data,"mean"]),dims=c(10,10))

projgrid <- inla.mesh.projector(mesh, dims=c(200,200))
xmean <- inla.mesh.project(projgrid, m$summary.random$s$mean)
image(xmean,asp=2)

res<-inla.spde2.result(m,"spatial",spde)

plot(res[["marginals.range.nominal"]][[1]], type = "l",main = "Nominal range, posterior density")


par(mfrow=c(1,2))

ma<-max(c(m$parameters$sd$prior[,2],m$parameters$sd$posterior[,2]))
plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density',lty=2,xlim=c(0,5),ylim=c(0,ma))
lines(m$parameters$sd$posterior,lty=1)
legend("topright", lty=2:1, legend=c("prior","posterior"))

ma<-max(c(m$parameters$range$prior[,2],m$parameters$range$posterior[,2]))
plot(m$parameters$range$posterior,type="l",xlim = c(0,50*1000),xlab='range (m)', ylab='density',lty=1,ylim=c(0,ma))
lines(m$parameters$range$prior,lty=2)
legend("topright", lty=2:1, legend=c("prior","posterior"))



#######################################
### glgm

fit<-glgm(PA~VEGZONSNA+Populati_2+Road_dens+MDC+dom,
          data=occs,
          grid=20,
          covariates=NULL, 
          family="binomial", 
          buffer=10000,
          shape=1,
          priorCI=list(sd=c(0.1,4),range=c(5000,100000)),
          control.compute=list(waic=TRUE,dic=TRUE,mlik=TRUE),#,config=TRUE),
          #control.predictor=list(compute=TRUE,link=1)
)

summary(fit$inla)

plot(fit$raster$predict.invlogit)

#summary(inla.rerun(fit$inla))

m<-fit

par(mfrow=c(1,2))

ma<-max(c(m$parameters$sd$prior[,2],m$parameters$sd$posterior[,2]))
plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density',lty=2,xlim=c(0,5),ylim=c(0,ma))
lines(m$parameters$sd$posterior,lty=1)
#lines(seq(0,10,by=0.1),dlgamma(seq(0,10,by=0.1),m$parameters$sd$params.intern$param[1],m$parameters$sd$params.intern$param[2]),col="blue")
legend("topright", lty=2:1, legend=c("prior","posterior"))

ma<-max(c(m$parameters$range$prior[,2],m$parameters$range$posterior[,2]))
plot(m$parameters$range$posterior,type="l",xlim = c(0,500*1000),xlab='range (m)', ylab='density',lty=1,ylim=c(0,ma))
lines(m$parameters$range$prior,lty=2)
#lines(dgamma(seq(0,50000,by=10),m$parameters$range$params.intern[1],m$parameters$range$params.intern[2]),col="red")
legend("topright", lty=2:1, legend=c("prior","posterior"))


################################
### lgcp




# is it better to do a lgcp for a point process

ds2<-occs[occs$PA==1,]

covList<-with(ds2@data,list(Populati_2=log(Populati_2),Road_dens=log(Road_dens),MDC=log(MDC)))

ds2$Road_dens2<-log(ds2$Road_dens)

r <- raster(ncol = 100, nrow = 200, ext = extent(ds2))
r <- rasterize(ds2, r, field = 1, fun = "count", background = 0)
plot(r)

roads<-rasterize(ds2[,"Road_dens2"],r, field = "Road_dens2", fun = mean, background = 0)
plot(roads)

fit<-lgcp(formula=~Road_dens2,
          data=ds2[1:50,],
          grid=20,
          covariates=list(Road_dens2=roads),
          family="binomial", 
          buffer=100000,
          shape=1,
          priorCI=list(sd=c(0.1,4000),range=c(20000,100000)),#,
          verbose=TRUE
          #control.compute=list(waic=TRUE,dic=TRUE,mlik=TRUE)#,config=TRUE),
          #control.predictor=list(compute=TRUE,link=1)
)




# a log-Gaussian Cox process example

library(inlabru)

data(gorillas)
# Use tutorial setting and thus empirical Bayes for faster inference
init.tutorial()
# Plot the Gorilla nests, the mesh and the survey boundary
ggplot() +
  gg(gorillas$mesh) +
  gg(gorillas$nests) +
  gg(gorillas$boundary) +
  coord_fixed()
# Define SPDE prior
matern <- inla.spde2.pcmatern(gorillas$mesh,
                              prior.sigma = c(0.1, 0.01),
                              prior.range = c(5, 0.01))
# Define domain of the LGCP as well as the model components (spatial SPDE effect and Intercept)
cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept
# Fit the model
fit <- inlabru:::lgcp(cmp, gorillas$nests, samplers = gorillas$boundary)
# Predict the spatial intensity surface
lambda <- predict(fit, pixels(gorillas$mesh), ~ exp(mySmooth + Intercept))
# Plot the intensity
ggplot() +
  gg(lambda) +
  gg(gorillas$mesh) +
  gg(gorillas$nests) +
  gg(gorillas$boundary) +
  coord_fixed()
## End(Not r


library(raster)
library(spatstat)

sweden <- raster:::getData("GADM", country = "SWE", level = 1)  
sweden<-spTransform(sweden,CRS(proj4string(occs)))

plot(sweden)
plot(ds2,add=TRUE)

o<-owin(xrange=bbox(ds2)[1,],yrange=bbox(ds2)[2,])
X<-ppp(coordinates(ds2)[,1],coordinates(ds2)[,2],window=o)



# inhomogeneous pattern of maples
#X <- unmark(split(lansing)$maple)

# (1) intensity function estimated by model-fitting
# Fit spatial trend: polynomial in x and y coordinates
fit <- ppm(X, ~ polynom(x,y,4), Poisson(),data=ds2@data)
# (a) predict intensity values at points themselves,
#     obtaining a vector of lambda values
lambda <- predict(fit, locations=X, type="trend")
# inhomogeneous K function
Ki <- Kinhom(X, lambda, nlarge=2000)
Ki <- Ginhom(X, lambda, nlarge=2000)
#Ki <- Kest(X, lambda)
plot(Ki,xlim=c(0,5000))
  
  

# method for point patterns
kppm(redwood, ~1, "Thomas")
# method for formulas
kppm(redwood ~ 1, "Thomas")

kppm(redwood ~ 1, "Thomas", method="c")
kppm(redwood ~ 1, "Thomas", method="p")

kppm(redwood ~ x, "MatClust") 
kppm(redwood ~ x, "MatClust", statistic="pcf", statargs=list(stoyan=0.2)) 
kppm(redwood ~ x, cluster="Cauchy", statistic="K")
kppm(redwood, cluster="VarGamma", nu = 0.5, statistic="pcf")

# LGCP models
kppm(redwood ~ 1, "LGCP", statistic="pcf")
if(require("RandomFields")) {
  k<-kppm(redwood ~ x, "LGCP", statistic="pcf",
       model="matern", nu=0.3,
       control=list(maxit=10))
}
  
  
  



