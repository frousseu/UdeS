
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

load("~/UdeS/Consultation/GStetcher/Doc/glgm_non-LFY.RData")
load("~/UdeS/Consultation/GStetcher/Doc/glgm_LFY.RData")

d<-glm.sp

newdat<-data.frame()

ds<-d
coordinates(ds)<-~Ostkoordin+Norrkoordi
proj4string(ds)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(ds)<-"+init=epsg:3006"

plot(ds,col=alpha(ifelse(ds$PA==1,"red","blue"),0.25),pch=16)

m<-glm(PA~VEGZONSNA+Populati_2+Road_dens+MDC+dom,data=d,family=binomial)

coords <- coordinates(ds)
v<-variog(coords=coords,data=resid(m),breaks=seq(0,250000,by=500))
fitv<-variofit(v,ini.cov.pars=c(0.2,30000),cov.model="exponential",fix.nugget=FALSE,nugget=0.125,fix.kappa=FALSE,kappa=0.45)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)") 
lines(fitv)

hist(d$Year)

swe <- raster::getData("GADM", country = "SWE", level = 1)
swe<-spTransform(swe,proj4string(ds))
#plot(swe)
#plot(ds,add=TRUE)

prdomain <- inla.nonconvex.hull(coordinates(ds),convex=-0.02, resolution = c(100, 100))

mesh<-inla.mesh.2d(loc=coordinates(ds),max.edge=c(100000,200000),offset=c(20000,50000),cutoff=15000,boundary=prdomain)
plot(mesh,asp=1)
points(ds,pch=16,cex=0.25,col="red")

#spde<-inla.spde2.matern(mesh,alpha=2)
spde<-inla.spde2.pcmatern(mesh,prior.range=c(100000,0.9),prior.sigma=c(3,0.1))

#m<-inla(model,data=list(y=d$PA,intercept=rep(1,spde$n.spde),spatial=1:spde$n.spde),control.predictor=list(A=A,compute=TRUE),family="binomial")

g<-makegrid(ds,n=3000)
g<-SpatialPoints(g,proj4string=CRS(proj4string(ds)))
g<-SpatialPixels(g)
gbuff<-gBuffer(gBuffer(gBuffer(ds,width=1000),width=5000),width=40000)
o<-over(g,gbuff)
g<-g[!is.na(o),]

plot(gbuff)
plot(g,add=TRUE)
plot(ds,add=TRUE)


s.index<-inla.spde.make.index(name="spatial",n.spde=spde$n.spde)

A<-inla.spde.make.A(mesh=mesh,loc=coordinates(ds))
Ap<-inla.spde.make.A(mesh=mesh,loc=coordinates(g))

stack.est<-inla.stack(data=list(y=d$PA),A=list(A),effects=list(c(s.index,list(intercept=1))),tag="est")
stack.latent<-inla.stack(data=list(xi=NA),A=list(Ap),effects=list(s.index),tag="latent")
stack.pred<-inla.stack(data=list(y=NA),A=list(Ap),effects=list(c(s.index,list(intercept=1))),tag="pred")

full.stack<-inla.stack(stack.est,stack.latent,stack.pred)

model<-y~-1+intercept+f(spatial,model=spde)

m<-inla(model,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),family="binomial",control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE))

index.est<-inla.stack.index(full.stack,tag="est")$data
index.latent<-inla.stack.index(full.stack,tag="latent")$data
index.pred<-inla.stack.index(full.stack,tag="pred")$data

p<-m$summary.fitted.values[index.pred,"0.5quant"]
p<-inla.link.invlogit(m$summary.linear.predictor[index.pred,"mean"])

gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))

brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(length(brks)-1,rev(brewer.pal(11,"RdYlGn")))
#plot(gp,breaks=brks,col=cols,at=pretty(brks,10))

plot(gp,col=cols)
points(ds,col=alpha(ifelse(d$PA==1,"red","black"),0.35),pch=16,cex=0.3)

image(inla.mesh.project(mesh,field=m$summary.fitted.values[inla.stack.index(full.stack,tag="latent")$data,"mean"]),dims=c(10,10))

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
          data=ds,
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

ds2<-ds[ds$PA==1,]

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
sweden<-spTransform(sweden,CRS(proj4string(ds)))

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
  
  
  



