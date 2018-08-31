
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

#######################################
#######################################
#######################################
#######################################

######################################################################
### function to create a sequence from the range of values in a vector
toseq<-function(x,n=100){
  if(is.numeric(x)){
    r<-range(x,na.rm=TRUE)
    seq(r[1],r[2],length.out=n)
  }else{
    sort(unique(x))  
  }
}

#####################################################################################
### function to create a data.frame for predictions for each variable in a data.frame
newdata<-function(x,v=names(x),n=100,fun=mean,list=FALSE){
  
  # returns the mean or the levels
  mm<-function(y){
    if(is.numeric(y)){
      fun(y)
    }else{
      names(rev(sort(table(y))))[1]
    }
  }
  ## possibly a bug with the function did not save the last part
  ans<-lapply(v,function(i){
    if(n==1L){
      val<-mm(x[,i])
    }else{
      val<-toseq(x[,i],n=n)
    }
    l<-lapply(x[,setdiff(names(x),i),drop=FALSE],mm)
    res<-data.frame(val,as.data.frame(l),stringsAsFactors=FALSE)
    names(res)[1]<-i
    if(list){ # inefficient, should not be turned to data.frame if list=TRUE
      as.list(res)
    }else{
      res
    } 
  })  
  names(ans)<-v
  
  if(length(v)==1L){
    unlist(ans,recursive=FALSE)
  }else{
    ans
  }
  
}

#############
### load data
load("~/UdeS/Consultation/GStetcher/Doc/LLF_size.RData")

size<-llf.size

#####################################################################
### convert high_name to factor and log transform the population size
size$high_name<-as.factor(size$high_name)
size$logPop_2017<-log(size$Pop_2017+0.2)

#############################################
### build a spatial object with the locations
sizes<-size
coordinates(sizes)<-~Longitude+Latitude
proj4string(sizes)<-"+init=epsg:4326"
prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sizes<-spTransform(sizes,CRS(prj))

#########################################
### build a simple glm to compare results
m1 <- glm (Area ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = Gamma(link = "log"), data = na.omit(size))
par(mfrow=c(3,3),mar=c(4,4,3,3))
visreg(m1,scale="response")
par(mfrow=c(1,1))

####################################################################
### look at the variogram with the residuals from the previous model
coords <- coordinates(sizes)
v<-variog(coords=coords,data=resid(m1),breaks=seq(0,200000,by=1000),max.dist=200000,bin.cloud=TRUE)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire sizeurrence)",type="b") 

##############
### get sweden
swe <- raster::getData("GADM", country = "SWE", level = 0)
swe<-spTransform(swe,proj4string(sizes))

#############################################
### build a mesh and use sweden as a boundary
# smaller values put in here will make a more precise grid, but will take longer to run
mesh<-inla.mesh.2d(loc=coordinates(sizes),max.edge=c(20000,200000),offset=c(20000,50000),cutoff=20000,boundary=swe)
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
  Area ~ -1 + intercept + FWI + f(spatial,model=spde),
  Area ~ -1 + intercept + VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI + f(spatial,model=spde)
)

###################
### model selection

# this runs every model in the model set

ml<-vector(mode="list",length=length(modell))

for(i in seq_along(modell)){
  # list of variables
  v<-setdiff(all.vars(modell[[i]]),c("Area","intercept","spatial","spde"))
  # build the data stack
  stack.est<-inla.stack(data=list(Area=size$Area),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(size[,v,drop=FALSE])),tag="est")
  # run the model with the eb strategy for faster runs (more approximate)
  ml[[i]]<-inla(modell[[i]],data=inla.stack.data(stack.est),control.predictor=list(A=inla.stack.A(stack.est)),family="gamma",control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE,return.marginals=FALSE),control.inla=list(strategy='gaussian',int.strategy="eb"))
  # print iterations
  print(paste(" ",i,"/",length(ml)," "))
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
n<-100 # number of divisions in generated values for the focus variable
Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(312180,6342453),ncol=2)[rep(1,n),,drop=FALSE]) # the graphs are build using a random points in the area

################################################
### build newdata with variable values to submit
v<-setdiff(all.vars(bmodel),c("Area","intercept","spatial","spde"))
lp<-newdata(x=size[,v,drop=FALSE],v=v,n=n,fun=median,list=TRUE)
lpmed<-lapply(newdata(x=size[,v,drop=FALSE],v=v,n=1,fun=median,list=TRUE)[[1]],function(i){rep(i,length(g))})

########################################################
### bind the data stack for the estimate and for the map
stack.est<-inla.stack(data=list(Area=size$Area),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(size[,v,drop=FALSE])),tag="est")
stack.map<-inla.stack(data=list(Area=NA),A=list(Ap,1),effects=list(c(s.index,list(intercept=1)),lpmed),tag="map")
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
  stack<-inla.stack(data=list(Area=NA),A=list(AA,1),effects=list(c(s.index,list(intercept=1)),lp[[v[i]]]),tag=v[i])     
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
m<-inla(bmodel,Ntrials=1,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),family="gamma",control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.inla=list(strategy='gaussian',int.strategy="eb"))


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
p<-m$summary.fitted.values[index[["map"]],"mean"]
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
  plot(lp[[v[i]]][[1]],p[,2],type="l",ylim=c(0,50),xlab=v[i],font=2,ylab="",lty=1)
  if(nrow(p)==n){
    lines(lp[[v[i]]][[1]],p[,1],lty=3)
    lines(lp[[v[i]]][[1]],p[,3],lty=3)
  }else{
    segments(x0=as.integer(lp[[v[i]]][[1]]),x1=as.integer(lp[[v[i]]][[1]]),y0=p[,1],y1=p[,3],lty=3)
  }
  points(size[,v[i]],size$Area,pch=16,col=gray(0,0.15))
}
mtext("Fire size in ha",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)


######################################################
### model checking and posterior predictive checks

# from bayesian regression modeling with INLA, Faraway frpom google books
# the same thing is shown on Blangardio et Cameletti on p. 168
# something does not work cause response is 1 or 0 and there are NaN in predictions
post.predicted.pval<-vector(mode="numeric",length=nrow(size))
for(i in 1:nrow(size)){
  post.predicted.pval[i]<-inla.pmarginal(q=size$Area[i],marginal=m$marginals.fitted.values[[i]])
}
hist(post.predicted.pval,main="",breaks=10,xlab="Posterior predictive p-value")


### from haakon bakka, BTopic112
samples<-inla.posterior.sample(2000,m)
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
  rbinom(length(i),size=1,prob=inla.link.invlogit(i))  
})
o<-createDHARMa(simulatedResponse=matprob,observedResponse=size$Area,fittedPredictedResponse=prob,integerResponse=TRUE)
par(mfrow=c(2,2))
plot(o,quantreg=TRUE)
#hist(o$scaledResiduals)

### check same for glm using the two methods
mod1 <- glm(Area ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size, family=Gamma(link="log"))
simulationOutput <- simulateResiduals(fittedModel = mod1)
plot(simulationOutput,quantreg=TRUE)
s<-simulateResiduals(mod1)
o<-createDHARMa(simulatedResponse = s[["simulatedResponse"]],observedResponse = size$Area, fittedPredictedResponse = fitted(mod1),integerResponse = T)
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

m <- glm (logArea ~ VEGZONSNA + WtrUrb_km + Pop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = gaussian(link = "identity"), data = size)
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