
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
library(DHARMa)
library(ROCR)
library(doParallel)
library(foreach)
library(fields)
library(viridisLite)
library(MASS)

#######################################
#######################################
#######################################
#######################################



#############
### load data
load("~/UdeS/Consultation/GStetcher/Doc/LLF_size.RData")
size<-llf.size
size<-size[sample(1:nrow(size),2000),]

#####################################################################
### source newdata and toseq

source("https://raw.githubusercontent.com/frousseu/UdeS/master/GStecher/newdata.R")

#####################################################################
### convert high_name to factor and log transform the population size
size$high_name<-as.factor(size$high_name)
size$logPop_2017<-log(size$Pop_2017+0.2)

##########################################################################################
### Use a boxcox transformation determined from the most complete linear model

BoxCox<-function(x,lambda=0){
  if(lambda==0){
    log(x)
  }else{
    ((x^lambda)-1)/lambda  
  }
}

BoxCoxI<-function(x,lambda=0){
  if(lambda==0){
    exp(x)
  }else{
    ((x*lambda)+1)^(1/lambda)  
  }
}

#trans<-BoxCox
#transI<-BoxCoxI

trans<-identity
transI<-identity

bc<-boxcox(lm (Area ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size))
lambda<-bc$x[which.max(bc$y)]
size$tArea<-trans(size$Area,lambda)
size$tArea<-size$Area

#############################################
### build a spatial object with the locations
sizes<-size
coordinates(sizes)<-~Longitude+Latitude
proj4string(sizes)<-"+init=epsg:4326"
prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sizes<-spTransform(sizes,CRS(prj))

#########################################
### build a simple glm to compare results
m1 <- lm (tArea ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI,data = na.omit(size))
#m1 <- glm (Area ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = Gamma(link = "log"), data = na.omit(size))
par(mfrow=c(3,3),mar=c(4,4,3,3))
visreg(m1,scale="response",trans=function(i){i^(1/lambda)},ylim=c(0,2))
par(mfrow=c(1,1))

####################################################################
### look at the variogram with the residuals from the previous model
coords <- coordinates(sizes)
v<-variog(coords=coords,data=resid(m1),breaks=seq(0,100000,by=1000),max.dist=100000,bin.cloud=TRUE)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire sizeurrence)",type="b") 

##############
### get sweden
swe <- raster::getData("GADM", country = "SWE", level = 0)
swe<-spTransform(swe,proj4string(sizes))

#############################################
### build a mesh and use sweden as a boundary
# smaller values put in here will make a more precise grid, but will take longer to run
mesh<-inla.mesh.2d(loc=coordinates(sizes),max.edge=c(5000,15000),offset=c(5000,15000),cutoff=5000,boundary=swe)
plot(mesh,asp=1)

#############################
### build spde with pc priors
spde<-inla.spde2.pcmatern(mesh,prior.range=c(100000,0.9),prior.sigma=c(3,0.1))

###########################################################
### build the raster/grid that will be used for predictions
g<-makegrid(swe,n=20000)
g<-SpatialPoints(g,proj4string=CRS(proj4string(sizes)))
#o<-over(as(g,"SpatialPolygons"),swe) # makes sure pixels touching are included too, does not change much when the grid gets small
o<-over(g,swe)
g<-SpatialPixels(g)
g<-g[apply(o,1,function(i){!all(is.na(i))}),]

############################################
### make spatial index to link with the spde
s.index<-inla.spde.make.index(name="spatial",n.spde=spde$n.spde)

###########################
### make observation matrix
A<-inla.spde.make.A(mesh=mesh,loc=coordinates(sizes))


#############
### model set

modell<-list(
  tArea ~ 0 + intercept + FWI + f(spatial,model=spde),
  tArea ~ 0 + intercept + VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI + f(spatial,model=spde)
)

###################
### model selection

# this runs every model in the model set

registerDoParallel(min(detectCores()-1,length(modell))) 
getDoParWorkers()

#ml<-vector(mode="list",length=length(modell))

ml<-foreach(i=seq_along(modell),.packages=c("stats","INLA"),.verbose=TRUE) %dopar% {
  # list of variables
  v<-setdiff(all.vars(modell[[i]]),c("tArea","intercept","spatial","spde"))
  # build the data stack
  spde<-spde # this only to make sure spde is exported to the nodes
  stack.est<-inla.stack(data=list(tArea=size$tArea),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(size[,v,drop=FALSE])),tag="est")
  # run the model with the eb strategy for faster runs (more approximate)
  m<-inla(modell[[i]],data=inla.stack.data(stack.est),control.predictor=list(A=inla.stack.A(stack.est)),control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE,return.marginals=FALSE),control.inla=list(strategy='gaussian',int.strategy="eb"),num.threads=1)
  # print iterations
  m
  #print(paste(" ",i,"/",length(ml)," "))
}

#######################################
### compute dic and waic values
dic<-sapply(ml,function(i){i$dic$dic})
waic<-sapply(ml,function(i){i$waic$waic})

#################################################
### order models according to one of the criteria
sel<-data.frame(waic,model=sapply(modell,function(i){format(deparse(i,width.cutoff=400))}))
sel<-sel[order(sel[,1]),]
head(sel)

####################################
### predictions

#######################
### find the best model
b<-which.min(waic)
modell[[b]]
bmodel<-modell[[b]]

###################################################################
### build prediction matrices for the map and the prediction graphs
Ap<-inla.spde.make.A(mesh=mesh,loc=coordinates(g))
n<-50 # number of divisions in generated values for the focus variable
Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(312180,6342453),ncol=2)[rep(1,n),,drop=FALSE]) # the graphs are build using a random points in the area

################################################
### build newdata with variable values to submit
v<-setdiff(all.vars(bmodel),c("tArea","intercept","spatial","spde"))
lp<-newdata(x=size[,v,drop=FALSE],v=v,n=n,fun=median,list=TRUE)
lpmed<-lapply(newdata(x=size[,v,drop=FALSE],v=v,n=1,fun=median,list=TRUE)[[1]],function(i){rep(i,length(g))})

########################################################
### bind the data stack for the estimate and for the map
stack.est<-inla.stack(data=list(tArea=size$tArea),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(size[,v,drop=FALSE])),tag="est")
stack.map<-inla.stack(data=list(tArea=NA),A=list(Ap,1),effects=list(c(s.index,list(intercept=1)),lpmed),tag="map")
full.stack<-inla.stack(stack.est,stack.map)

#######################################
### add a stack for each focus variable
for(i in seq_along(v)){
  le<-length(lp[[v[i]]][[1]])
  if(le!=n){
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(819006,6545844),ncol=2)[rep(1,le),,drop=FALSE]) # for categorical variables
  }else{
    AA<-Apn # for numerical variables
  }
  stack<-inla.stack(data=list(tArea=NA),A=list(AA,1),effects=list(c(s.index,list(intercept=1)),lp[[v[i]]]),tag=v[i])     
  full.stack<-inla.stack(full.stack,stack)
}

##########################################
### extract index for each stack
index.est<-inla.stack.index(full.stack,tag="est")$data
index.map<-inla.stack.index(full.stack,tag="map")$data
index<-list(est=index.est,map=index.map)
for(i in seq_along(v)){
  index<-c(index,list(inla.stack.index(full.stack,tag=v[i])$data))
}  
names(index)[3:length(index)]<-v

##################################################
### rerun best model with each variable to predict
#m<-inla(bmodel,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.inla=list(strategy='gaussian',int.strategy="eb"),num.threads=6)

m<-inla(bmodel,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.inla=list(strategy='gaussian',int.strategy="eb"),family="gp",control.family=list(control.link=list(quantile=0.98)),num.threads=6)

##################################
### build a relative frequency map

### raw predictions on probability of usage
par(mfrow=c(1,4),oma=c(0,5,0,5))
plot(swe,border=gray(0,0.25),lwd=0.01)
points(sizes,col=alpha("blue",0.2),pch=16,cex=30*size$Area/max(size$Area))
mtext("Fire size",side=4,font=2)

plot(swe,border=gray(0,0.25),lwd=0.01)
points(sizes,col=alpha("blue",0.2),pch=16,cex=log(size$Area))
mtext("log Fire size",side=4,font=2)

p<-transI(m$summary.fitted.values[index[["map"]],"mean"]) # the lambda is to back-transform on the original scale)
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Predicted fire size', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)

### sd
p<-m$summary.fitted.values[index[["map"]],"sd"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='sd of predicted fire size', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)


#####################################
### graphical predictions

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
for(i in seq_along(v)){
  p<-m$summary.fitted.values[index[[v[i]]],c("0.025quant","0.5quant","0.975quant")]
  p[]<-lapply(p,transI)
  plot(lp[[v[i]]][[1]],p[,2],type="l",ylim=c(0,100),xlab=v[i],font=2,ylab="",lty=1)
  if(nrow(p)==n){
    lines(lp[[v[i]]][[1]],p[,1],lty=3)
    lines(lp[[v[i]]][[1]],p[,3],lty=3)
  }else{
    segments(x0=as.integer(lp[[v[i]]][[1]]),x1=as.integer(lp[[v[i]]][[1]]),y0=p[,1],y1=p[,3],lty=3)
  }
  points(size[,v[i]],transI(size$tArea),pch=16,col=gray(0,0.15))
}
mtext("Fire size in ha",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)

#rq<-rq(tArea~FWI,data=size,tau=0.98)
#p<-predict(rq,data.frame(FWI=toseq(size$FWI,100)))
#lines(toseq(size$FWI,100),p,col="blue")



######################################################
### model checking and posterior predictive checks

# from bayesian regression modeling with INLA, Faraway frpom google books
# the same thing is shown on Blangardio et Cameletti on p. 168
# something does not work cause response is 1 or 0 and there are NaN in predictions
post.predicted.pval<-vector(mode="numeric",length=nrow(size))
for(i in 1:nrow(size)){
  post.predicted.pval[i]<-inla.pmarginal(q=size$tArea[i],marginal=m$marginals.fitted.values[[i]])
}
hist(post.predicted.pval,main="",breaks=10,xlab="Posterior predictive p-value")


### from haakon bakka, BTopic112
samples<-inla.posterior.sample(500,m)
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor"
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-t(matrix(unlist(samples.effect),byrow=T,nrow=length(samples.effect)))

### check with inla model
prob<-m$summary.fitted.values[index[["est"]],"0.5quant"]
matprob<-apply(s.eff,2,function(i){
  rnorm(length(i),size=1,prob=inla.link.invlogit(i))  
})
o<-createDHARMa(simulatedResponse=s.eff,observedResponse=size$tArea,fittedPredictedResponse=prob,integerResponse=TRUE)
par(mfrow=c(2,2))
plot(o,quantreg=TRUE)
#hist(o$scaledResiduals)

### check same for glm using the two methods
mod1 <- lm(tArea ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size)

mod1 <- glm(Area ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size, family=Gamma(link="log"))
mod1 <- lm(Area^(0.1) ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size)
#mod1 <- glm(Area ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size, family=tweedie(var.power = 2.5, link.power = 1))
simulationOutput <- simulateResiduals(fittedModel = mod1)
plot(simulationOutput,quantreg=TRUE)
s<-simulateResiduals(mod1)
o<-createDHARMa(simulatedResponse = s[["simulatedResponse"]],observedResponse = size$Area^(0.1), fittedPredictedResponse = fitted(mod1),integerResponse = F)
plot(o,quantreg=TRUE)
hist(o$scaledResiduals)


##########################################
### confusion matrix with training data

##########################################
### confusion matrix on hold-out data


##########################################
### confusion matrix on hold-out data

### code adapted from Myer et al. 2017 (spatiotemporal mosquitoes, https://doi.org/10.1002/ecs2.1854) supplementary material

### but see Boyce et al. 2002 pour le use-availability design

##########################################
### cross-validation/cpo/pit measures ?

#################################
###

#image(inla.mesh.project(mesh,field=m$summary.fitted.values[inla.stack.index(full.stack,tag="latent")$data,"mean"]),dims=c(10,10))
projgrid <- inla.mesh.projector(mesh, dims=c(500,500))
xmean <- inla.mesh.project(projgrid, m$summary.random$spatial$mean)
xsd <- inla.mesh.project(projgrid, m$summary.random$spatial$sd)
image(xmean,asp=2,col=heat.colors(100))
res<-inla.spde2.result(m,"spatial",spde)
plot(res[["marginals.range.nominal"]][[1]], type = "l",main = "Nominal range, posterior density")










#######################################
### PLAY
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################


load("~/UdeS/Consultation/GStetcher/Doc/LLF_size.RData")

size<-llf.size

sizes<-size
coordinates(sizes)<-~Longitude+Latitude
proj4string(sizes)<-"+init=epsg:4326"

prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sizes<-spTransform(sizes,CRS(prj))


#plot(ds,col=alpha(ifelse(ds$Area==1,"red","blue"),0.25),pch=16)

bc<-function(y,div=50){
  default<-par()
  lambda<-seq(-1,0.1,length.out=div)
  par(mfrow=c(ceiling(sqrt(div)),ceiling(sqrt(div))),mar=c(0,0,0,0))
  invisible(lapply(lambda,function(i){
    hist((y^i-1)/i,main="",xaxt="n",yaxt="n")
    mtext(paste("lambda",round(i,2)),side=3,line=-2,col=gray(0,0.5))
  }))
  par(default)
}
bc(size$Area)

size$high_name<-as.factor(size$high_name)
size$logArea<-log(size$Area)

m <- glm (Area ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = Gamma(link = "log"), data = size)
m <- glm (Area ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = Gamma(link = "log"), data = size)
#m <- rq (logArea ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size, tau=0.95)
#m <- glm (Area ~ 1, family = gaussian(link = "log"), data = d)
sizes$resid<-resid(m)

par(mfrow=c(3,3))
visreg(m,scale="response")
par(mfrow=c(1,1))

coords <- coordinates(sizes)
v<-variog(coords=coords,data=resid(m),breaks=seq(0,200000,by=1000),max.dist=200000)
#v<-variog(coords=coords,data=resid(m))
#v<-variogram(resid~1,data=ds)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b") 

simulationOutput <- simulateResiduals(fittedModel = m, n = 250, integerResiduals=FALSE)
plot(simulationOutput)
testDispersion(simulationOutput)

swe <- raster::getData("GADM", country = "SWE", level = 1)
swe<-spTransform(swe,proj4string(sizes))





















f<-size$Area
f2<-log(f)
brks<-exp(seq(log(min(f)-0.00001),log(max(f)+0.00001),length.out=20))
cl<-cut(f2,brks)
f2[is.na(cl)]
tab<-table(cl)
barplot(tab,names.arg=names(tab),las=2)



bc<-boxcox(lm (Area ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size))
lambda<-bc$x[which.max(bc$y)]

#lambda<--1000.0
f<-BoxCox(size$Area,lambda)
hist(f)

visreg(mod1,trans=function(i){BoxCoxI(i,lambda)})



### build the raster/grid that will be used for predictions
g<-makegrid(swe,n=200)
g<-SpatialPoints(g,proj4string=CRS(proj4string(sizes)))
o<-over(g,swe)
g<-SpatialPixels(g)
g<-g[apply(o,1,function(i){!all(is.na(i))}),]
plot(swe)
plot(g,add=TRUE)


cfire<-function(x,show=FALSE,n=20,mode=TRUE){ # gives characteristics fire zei according to  Lehsten et al 2013
  x2<-log(x)
  brks<-exp(seq(log(min(x)-0.00001),log(max(x)+0.00001),length.out=n))
  cl<-cut(x2,brks)
  x2[is.na(cl)]
  tab<-table(cl)
  if(show){
    barplot(tab,names.arg=names(tab),las=2)
  }
  if(mode){
    sapply(strsplit(gsub("\\(|\\]","",names(tab[which.max(tab)])),","),function(i){mean(as.numeric(i))})
  }else{
    tab  
  }
}

o<-over(sizes,g)
size$o<-o

s<-data.table(size)
s<-s[,.(cf=cfire(Area)),by=o]
hist(s$cf)


barplot(cfire(size$Area,mode=FALSE))


hist(BoxCox(size$Area+0.1,lambda))



s














