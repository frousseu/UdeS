
library(sp)
library(rgdal)
library(scales)
library(gstat)
library(geoR)
library(INLA)
library(rgeos)
library(FRutils)
library(RColorBrewer)
library(visreg)
library(quantreg)
library(DHARMa)
library(ROCR)
library(doParallel)
library(foreach)
library(fields)
library(viridisLite)
library(MASS)
library(maptools)
library(ggeffects)
library(gridExtra)
library(raster)

#######################################
#######################################
#######################################
#######################################

cat("\014")

#############
### load data
load("~/UdeS/Consultation/GStetcher/Doc/MSB_size.RData")
#size<-size
#size<-llf.size
size<-size[sample(1:nrow(size),1000),]

vars<-c("Total","Road1k","Pp_1000","urbwtr1k","frbreak1k","Ag_1000","h__1000","VEGZONS","FWI")

# some names seem to have been abbreviated by ArcGIS in the database
size$Road1k<-size$Rod_dns
size$urbwtr1k<-size$urb_wtr_onesR
size$frbreak1k<-size$firebrk_onesR
size$Total<-size$Total/10000 # put fire size in ha

plot(size[,vars])

#####################################################################
### source newdata and toseq

source("https://raw.githubusercontent.com/frousseu/UdeS/master/GStecher/newdata.R")

#####################################################################
### convert factors and log transform the population size

#size$urbwtr1k<-as.factor(size$urbwtr1k)
#size$frbreak1k<-as.factor(size$frbreak1k)
size$h__1000<-as.factor(size$h__1000)
size$VEGZONS<-as.factor(gsub(" ","_",size$VEGZONS))
size$Pp_1000<-log(size$Pp_1000+1)

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

bc<-boxcox(lm (Total ~ Road1k + Pp_1000 + urbwtr1k + frbreak1k + Ag_1000 + h__1000 + VEGZONS + FWI, data = size))
lambda<-bc$x[which.max(bc$y)]
size$ttTotal<-BoxCox(size$Total,lambda) # the final version won't use any transformation I think
size$tTotal<-size$Total

#############################################
### build a spatial object with the locations
sizes<-size
coordinates(sizes)<-~LongRT90+LatRT90
proj4string(sizes)<-"+init=epsg:4326"
prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sizes<-spTransform(sizes,CRS(prj))

#########################################
### build a simple glm to compare results
m1 <- lm (ttTotal ~ Road1k + Pp_1000 + urbwtr1k + frbreak1k + Ag_1000 + h__1000 + VEGZONS + FWI,data = size)
#m1 <- glm (Area ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, family = Gamma(link = "log"), data = na.omit(size))
par(mfrow=c(3,3),mar=c(4,4,3,3))
visreg(m1,scale="response",trans=BoxCoxI,rug=FALSE)
par(mfrow=c(1,1))

####################################################################
### look at the variogram with the residuals from the previous model
coords <- coordinates(sizes)
v<-variog(coords=coords,data=resid(m1),breaks=seq(0,100000,by=500),max.dist=100000,bin.cloud=TRUE)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire sizeurrence)",type="b",xlim=c(0,100000)) 

##############
### get sweden
swe <- raster::getData("GADM", country = "SWE", level = 0)
swe<-spTransform(swe,proj4string(sizes))

bb<-bbox(swe)
swediv<-elide(swe,bb=bb,scale=TRUE) # this scales the data to [0,1] using the bb while maintaining the aspect
sizesdiv<-elide(sizes,bb=bb,scale=TRUE)

#############################################
### build a mesh and use sweden as a boundary
# smaller values put in here will make a more precise grid, but will take longer to run

w<-which.max(apply(bb,1,max))
div<-diff(abs(bb[w,])) # this is a scale factor to reduce absolute values to smaller values by division
me<-10000/div # this is the max edge in meters divided by the value to scale this to [0,1]
out_fac<-3 # this the max edge in the buffer area

mesh<-inla.mesh.2d(loc=coordinates(sizesdiv),max.edge=c(me,me*out_fac),offset=c(me,me*out_fac),cutoff=me,boundary=swediv)
plot(swediv,col="tomato3",axes=TRUE)
plot(mesh,asp=1,add=TRUE)

#############################
### build spde with pc priors
spde<-inla.spde2.pcmatern(mesh,prior.range=c(50000/div,0.5),prior.sigma=c(3,0.1))

####################################################################
### sensible priors on factors and intercept

fac<-sapply(size,is.factor)
fac<-names(fac[fac])
fac<-fac[fac%in%vars]
fac<-sapply(fac,function(i){levels(size[,i])})
fac<-unlist(sapply(seq_along(fac),function(i){paste0(names(fac)[i],fac[[i]])}))
vals<-rep(1/5^2,length(fac))
names(vals)<-fac
vals<-as.list(vals)
vals<-c(vals,list(intercept=1/5^2,default=0.001)) 
control.fixed<-list(prec=vals,mean=list(intercept=0,default=0),expand.factor.strategy = "model.matrix")

###########################################################
### build the raster/grid that will be used for predictions
g<-makegrid(swediv,n=5000)
g<-SpatialPoints(g,proj4string=CRS(proj4string(sizesdiv)))
#o<-over(as(g,"SpatialPolygons"),swe) # makes sure pixels touching are included too, does not change much when the grid gets small
o<-over(g,swediv)
g<-SpatialPixels(g)
g<-g[apply(o,1,function(i){!all(is.na(i))}),]

############################################
### make spatial index to link with the spde
s.index<-inla.spde.make.index(name="spatial",n.spde=spde$n.spde)

###########################
### make observation matrix
A<-inla.spde.make.A(mesh=mesh,loc=coordinates(sizesdiv))


#############
### model set

modell<-list(
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + FWI + f(spatial,model=spde),
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + VEGZONS + FWI + f(spatial,model=spde),
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + Ag_1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + urbwtr1k + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + h__1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + h__1000 + VEGZONS + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + h__1000 + Ag_1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + frbreak1k + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + h__1000 + urbwtr1k + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + h__1000 + frbreak1k + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + VEGZONS + Ag_1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + urbwtr1k + Ag_1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + frbreak1k + Ag_1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + frbreak1k + Ag_1000 + h__1000 + FWI + f(spatial,model=spde),
  #Total ~ -1 + intercept + Road1k + Pp_1000 + frbreak1k + Ag_1000 + VEGZONS + FWI + f(spatial,model=spde),
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + urbwtr1k + Ag_1000 + VEGZONS + FWI + f(spatial,model=spde),
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + urbwtr1k + Ag_1000 + h__1000 + FWI + f(spatial,model=spde),
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + urbwtr1k + Ag_1000 + h__1000 + VEGZONS + FWI + f(spatial,model=spde),
  tTotal ~ -1 + intercept + Road1k + Pp_1000 + urbwtr1k + frbreak1k + Ag_1000 + h__1000 + VEGZONS + FWI + f(spatial,model=spde)
)


### this is to tunr the formula using dummy variables for factors
mmatrix<-function(i,dat){
  if(!is.list(i)){
    i<-list(i)  
  }
  lapply(i,function(j){
    mod<-as.character(j)
    mod[3]<-gsub(" \\+ f\\(spatial, model = spde\\)","",mod[3]) # remove spatial effect
    mod[3]<-gsub("-1 \\+ intercept \\+ ","",mod[3]) # remove spatial effect and intercept notation
    mod<-as.formula(paste0(mod[c(1,3)],collapse=" "))
    model.matrix(mod,dat)[,-1]
  })
}
mm<-mmatrix(modell,size)


### this is the model with the dummy variable from the model.matrix as suggested by the pdf E Krainski on using factors
modellmm<-lapply(mm,function(i){
  as.formula(paste("tTotal ~ -1 + intercept +",paste(dimnames(i)[[2]],collapse=" + "),"+ f(spatial, model = spde)"))
})


###################
### model selection

# this runs every model in the model set

registerDoParallel(min(detectCores()-1,length(modell))) 
getDoParWorkers()

q<-0.99 # quantile to estimate
lambda<-10 # pc (like?) prior for xi parameter of gp (see below for a graph of prior)
hyper.gp <- list(theta = list(prior = "loggamma",param = c(1,lambda))) # c(1,15) is supposed to be the default prior. See below

ml<-foreach(i=seq_along(modell),.packages=c("stats","INLA"),.verbose=TRUE) %dopar% {
  # list of variables
  v<-setdiff(all.vars(modell[[i]]),c("tTotal","intercept","spatial","spde"))
  # build the data stack
  spde<-spde # this only to make sure spde is exported to the nodes
  stack.est<-inla.stack(data=list(tTotal=size$tTotal),A=list(A,1),effects=list(c(s.index,list(intercept=1)),data.frame(mm[[i]])),tag="est")
  # run the model with the eb strategy for faster runs (more approximate)
  m<-inla(modellmm[[i]],data=inla.stack.data(stack.est),control.predictor=list(A=inla.stack.A(stack.est)),control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE,return.marginals=FALSE),control.inla=list(strategy='simplified.laplace',int.strategy="eb"),family="gp",control.family=list(control.link=list(quantile=q),hyper=hyper.gp),control.fixed=control.fixed,num.threads=1)
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
#b<-length(modell)
modell[[b]]
bmodel<-modell[[b]]

###################################################################
### build prediction matrices for the map and the prediction graphs
Ap<-inla.spde.make.A(mesh=mesh,loc=coordinates(g))
n<-50 # number of divisions in generated values for the focus variable
Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,n),,drop=FALSE])

################################################
### build newdata with variable values to submit
v<-setdiff(all.vars(bmodel),c("tTotal","Total","intercept","spatial","spde"))
lp<-newdata(x=size[,v,drop=FALSE],v=v,n=n,fun=median,list=FALSE)
lp<-lapply(lp,function(i){mmatrix(bmodel,i)})
lpmed<-mmatrix(bmodel,newdata(x=size[,v,drop=FALSE],v=v,n=1,fun=median,list=FALSE)[[1]][rep(1,length(g)),])[[1]]

########################################################
### bind the data stack for the estimate and for the map
stack.est<-inla.stack(data=list(tTotal=size$tTotal),A=list(A,1),effects=list(c(s.index,list(intercept=1)),data.frame(mm[[b]])),tag="est")
stack.map<-inla.stack(data=list(tTotal=NA),A=list(Ap,1),effects=list(c(s.index,list(intercept=1)),data.frame(lpmed)),tag="map")
full.stack<-inla.stack(stack.est,stack.map)

#######################################
### add a stack for each focus variable
for(i in seq_along(v)){
  le<-nrow(lp[[v[i]]][[1]])
  if(le!=n){
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,le),,drop=FALSE]) # for categorical variables
  }else{
    AA<-Apn # for numerical variables
  }
  stack<-inla.stack(data=list(tTotal=NA),A=list(AA,1),effects=list(c(s.index,list(intercept=1)),data.frame(lp[[v[i]]][[1]])),tag=v[i])     
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

m<-inla(modellmm[[b]],data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.inla=list(strategy='simplified.laplace',int.strategy="eb"),family="gp",control.family=list(list(control.link=list(quantile=q),hyper=hyper.gp)),control.fixed=control.fixed,num.threads=7)


#####################################
### visualize spatial fields

xlim<-range(coordinates(sizesdiv)[,1])
ylim<-range(coordinates(sizesdiv)[,2])

proj<-inla.mesh.projector(mesh,xlim=xlim,ylim=ylim,dims=c(300,300))

mfield<-inla.mesh.project(projector=proj,field=m$summary.random[['spatial']][['mean']])
sdfield<-inla.mesh.project(projector=proj,field=m$summary.random[['spatial']][['sd']])

par(mfrow=c(1,2),mar=c(3,3,2,5))

image.plot(list(x=proj$x,y=proj$y,z=mfield),col=viridis(100),asp=1,main="Spatial field (log scale)") 
axis(1)
axis(2)
plot(swediv,add=TRUE,border=gray(0,0.5))
# not sure how to represent the q quantile form observations...
#plot(sizesdiv,pch=1,cex=0.1*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.13),add=TRUE)
#brks<-c(0.01,0.25,0.50,0.75,0.99)
#legend("topleft",pch=1,pt.cex=3*brks,col=gray(0,0.3),legend=brks,bty="n",title="Probability of location\nbeing an actual fire",inset=c(0.02,0.05))

image.plot(list(x=proj$x,y=proj$y,z=sdfield),col=viridis(100),asp=1,main="sd of spatial field (log scale)") 
axis(1)
axis(2)
plot(swediv,add=TRUE,border=gray(0,0.5))
#plot(sizesdiv,pch=1,cex=0.1*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.13),add=TRUE)



##################################
### build a relative frequency map

### this section is not useful or accurate unless we have values for each pixel

### raw predictions on probability of usage
par(mfrow=c(1,4),oma=c(0,5,0,5))
plot(swediv,border=gray(0,0.25),lwd=0.01,axes=TRUE)
points(sizesdiv,col=alpha("blue",0.2),pch=16,cex=30*size$Total/max(size$Total))
mtext("Fire size",side=4,font=2)

plot(swediv,border=gray(0,0.25),lwd=0.01,axes=TRUE)
points(sizesdiv,col=alpha("blue",0.2),pch=16,cex=sqrt(size$Total)/200)
mtext("sqrt Fire size / 200",side=4,font=2)

p<-transI(m$summary.fitted.values[index[["map"]],"mean"]) # the lambda is to back-transform on the original scale)
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swediv,axes=TRUE)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Predicted fire size', side=4, font=2, line=2.3))
plot(swediv,add=TRUE,border=gray(0,0.25),lwd=0.01)

### sd
p<-m$summary.fitted.values[index[["map"]],"sd"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swediv,axes=TRUE)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='sd of predicted fire size', side=4, font=2, line=2.3))
plot(swediv,add=TRUE,border=gray(0,0.25),lwd=0.01)


####################################################
### graphical predictions with spatial uncertainty

### this section is not that useful because it is a prediction for a given location, hence it includes uncertainty in the spatial field

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
for(i in seq_along(v)){
  p<-m$summary.fitted.values[index[[v[i]]],c("0.025quant","0.5quant","0.975quant")]
  p[]<-lapply(p,transI)
  dat<-data.frame(lp[[v[i]]][[1]])
  if(nrow(p)==n){
    plot(dat[[v[i]]],p[,2],type="l",ylim=c(0,100),xlab=v[i],font=2,ylab="",lty=1,yaxt="n")
    lines(dat[[v[i]]],p[,1],lty=3,lwd=1)
    lines(dat[[v[i]]],p[,3],lty=3,lwd=1)
    points(size[,v[i]],transI(size$tTotal),pch=16,col=gray(0,0.07))
  }else{
    plot(unique(sort(size[,v[i]])),p[,2],type="l",ylim=c(0,100),xlab=v[i],font=2,ylab="",lty=1,yaxt="n")
    segments(x0=as.integer(unique(sort(size[,v[i]]))),x1=as.integer(unique(sort(size[,v[i]]))),y0=p[,1],y1=p[,3],lty=3,lwd=2)
    points(jitter(as.integer(size[,v[i]]),fac=2),transI(size$tTotal),pch=16,col=gray(0,0.07))
  }
  axis(2,las=2)
}
mtext("Fire size in ha",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)


##########################################################
### generate predictions without spatial uncertainty

# page 263 in Zuur

nsims<-500
s<-inla.posterior.sample(nsims,m)
params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  match(i,row.names(s[[1]]$latent))  
}) 
nweights<-grep("spatial",row.names(s[[1]]$latent))

### this is to compare with a quantile model
#par(mfrow=c(round(sqrt(2*length(v)),0),ceiling(sqrt(2*length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
#mq <- rq (log(tTotal) ~ Road1k + Pp_1000 + urbwtr1k + frbreak1k + Ag_1000 + h__1000 + VEGZONS + FWI, data = na.omit(size),tau=q)
#na<-all.vars(mq$formula[[3]])
#grobs<-lapply(na,function(i){
#  if(is.factor(size[,i])){
#    plot(ggpredict(mq,terms=i),limits=c(0,100))  
#  }else{
#    plot(ggpredict(mq,terms=paste(i,"[n=50]")),limits=c(0,100)) 
#  }
#})
#grid.arrange(grobs=grobs,ncol=3)

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,3,2,2),oma=c(0,10,0,0))
for(k in seq_along(v)){
  p<-lapply(1:nsims,function(i){
    betas<-s[[i]]$latent[nparams]
    fixed<-cbind(intercept=1,lp[[v[k]]][[1]]) %*% betas
    ### this if we want a spatial part
    #wk<-s[[i]]$latent[nweights]
    #if(is.factor(size[,v[k]])){
    #  spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    #  spatial<-as.matrix(Apn) %*% wk
    #}
    #p<-exp(fixed+spatial)
    p<-exp(fixed)
    p
  })
  p<-do.call("cbind",p)
  p<-t(apply(p,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  if(nrow(lp[[v[k]]][[1]])==n){
    vals<-lp[[v[k]]][[1]][,v[k]]
    plot(vals,p[,2],type="l",ylim=c(0,100),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    points(size[,v[k]],size$tTotal,pch=1,col=gray(0,0.15))
    lines(vals,p[,2],lwd=2)
    lines(vals,p[,1],lty=3)
    lines(vals,p[,3],lty=3)
  }else{
    plot(unique(sort(size[,v[k]])),p[,2],type="l",ylim=c(0,100),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    points(jitter(as.integer(size[,v[k]])),size$tTotal,pch=1,col=gray(0,0.15))
    segments(x0=as.integer(unique(sort(size[,v[k]]))),x1=as.integer(unique(sort(size[,v[k]]))),y0=p[,1],y1=p[,3],lty=3)
  }
  axis(2,las=2)
}
mtext(paste("Fire size at the",q,"quantile (ha)"),outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)

###############################################
### range and sigma

### this is to show the posteriors of the spatial field

res <- inla.spde.result(m, "spatial", spde)
par(mfrow=c(2,1))
plot(res$marginals.range.nominal[[1]],
     type="l", main="Posterior density for range")
plot(inla.tmarginal(sqrt, res$marginals.variance.nominal[[1]]),
     type="l", main="Posterior density for std.dev.")
par(mfrow=c(1,1))


######################################################################
### check prior and posterior for shape parameter of genPareto

xi<-seq(0,1.2,by=0.001)
f<-function(xi,lambda=10){sqrt(2)*lambda*exp(-sqrt(2)*lambda*xi)}
1-integrate(f,lower=0,upper=0.5,lambda=lambda)$value
plot(xi,f(xi,lambda=lambda),type="l",yaxs="i",xaxs="i",ylim=c(0,max(f(xi,lambda=lambda))),lty=3)
lines(m$marginals.hyperpar[[1]][,1],m$marginals.hyperpar[[1]][,2],lwd=1)
legend("topright",lty=c(3,1),legend=c("Prior","Posterior"),bty="n",cex=1.5)


#################################################################
### generates samples from genPareto distributions

rgp <- function(n, sigma, eta, alpha, xi = 0.001){
  if (missing(sigma)) {
    stopifnot(!missing(eta) && !missing(alpha))
    sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) - 1.0)
  }
  return (sigma / xi * (runif(n)^(-xi) - 1.0))
}
#hist(rgp(n=10000,eta=2.35,alpha=0.98,xi=0.62),breaks=100)


######################################################
### model checking with DHARMa

### from haakon bakka, BTopic112
samples<-inla.posterior.sample(200,m)
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor" # not sure if should use APredictor or Predictor
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-do.call("cbind",samples.effect)
#unique(sapply(strsplit(rownames(m$summary.fitted.values),"\\."),function(i){paste(i[1:min(2,length(i))],collapse=" ")}))

### check with inla model
fitted<-m$summary.fitted.values[index[["est"]],"mean"]
matprob<-apply(s.eff,2,function(i){
  rgp(n=length(i),eta=i,alpha=q,xi=m$summary.hyperpar[1,1]) # this picks a fire size from the genPareto with estimated xi for each obs
})
o<-createDHARMa(simulatedResponse=matprob,observedResponse=size$tTotal,fittedPredictedResponse=fitted,integerResponse=TRUE)
par(mfrow=c(2,2))
plot(o,quantreg=TRUE)
#hist(o$scaledResiduals)


 
################################################################
### this shows the q quantile of observation compared to the simulated quantiles
# not sure if this is the right way to look at it

plot(density(log(size$tTotal)),ylim=range(density(s.eff[,1])$y))
abline(v=quantile(log(size$tTotal),q),col="red")
abline(v=log(quantile(size$tTotal,q)),col="red")
invisible(lapply(1:ncol(s.eff),function(i){
  lines(density(s.eff[,i]),col=gray(0.0,0.02))    
}))


#############################################
### 

for(i in seq_along(fitted)[1:100]){
  hist(as.vector(s.eff[i,]),xlim=range(c(as.vector(s.eff),log(size$tTotal[i]))))
  abline(v=log(size$tTotal[i]),col="red")
}

####################################################################################
### This compares the q quantile of observations to simulated observations from gp

brks<-0:ceiling(max(size$tTotal))
xlim<-c(0,50)
h<-hist(size$tTotal,breaks=brks,xlim=xlim,col="grey70",border="white",ylab="Counts")
abline(v=quantile(size$tTotal,q),lwd=2)
vals<-lapply(1:ncol(s.eff),function(i){
  #abline(v=mean(exp(s.eff[,i])),col=alpha("red",0.1))
  vals<-rgp(n=nrow(s.eff),eta=s.eff[,i],alpha=q,xi=m$summary.hyperpar[1,1])
  hi<-hist(vals,breaks=0:(ceiling(max(vals))),plot=FALSE)
  points(brks[-1]-diff(brks)/2,hi$counts[1:(length(brks)-1)],col=alpha("red",0.2),pch=16)
  abline(v=quantile(vals,q),col=alpha("red",0.2))
})

#### generate samples using VGAM functions
n<-1000
shape<-m$summary.hyperpar[1,1]
eta<-2.5
qs<-q
sigma<-xi*exp(eta)/((1-qs)^(-xi)-1)
xs<-seq(0,50,by=0.01)[-1]
si<-dgpd(xs,location=0,scale=sigma,shape=xi)
w<-which(si<0.005)[1]
plot(xs,si,type="l",yaxs="i",xaxs="i",xlim=c(0,xs[w]))
si<-rgpd(n=n,location=0,scale=sigma,shape=xi)
hist(si,freq=FALSE,add=TRUE,breaks=500)  
  
################################################################################
### plotting predicted pareto curves for a sample of predicted values from obs

#gpcdf<-function(y,sigma,xi){1-(1+xi*(y/sigma))^(-1/xi)} # CDF
gppdf<-function(y,sigma,xi){(1/sigma)*(1+xi*(y/sigma))^(-1*((1/xi)+1))} # CDF

samp<-100
eta<-m$summary.linear.predictor[index[["est"]],"mean"]
hist(eta)
xi<-m$summary.hyperpar[1,1]
sigma<-(xi*exp(median(eta)))/((1-q)^(-xi)-1)
va<-seq(0,10,by=0.01)
plot(va,gppdf(y=va,sigma=sigma,xi=xi),type="l",xaxs="i",yaxs="i",lwd=2,col="red")
quant<-lapply(sample(eta,samp),function(i){
  xi<-m$summary.hyperpar[1,1]
  sigma<-(xi*exp(i))/((1-q)^(-xi)-1)
  lines(va,gppdf(y=va,sigma=sigma,xi=xi),col=gray(0,0.1))
  integrate(gppdf,lower=0,upper=30,sigma=sigma,xi=xi)$value
})
h<-hist(size$tTotal,breaks=c(seq(0,par("usr")[2],length.out=50),max(size$tTotal)),plot=FALSE)
points(h$mids,rescale(h$density,c(0,1.5)),lwd=2,col="darkgreen")
hist(rgp(n=10000,eta=median(eta),alpha=q,xi=xi),breaks=100)
hist(unlist(quant)) ### histograms of % fires below 30 ha (upper in integrate)
#hist(size$tTotal,breaks=seq(0,max(size$tTotal),by=0.05),xlim=c(0,2))

#####

xi<-m$summary.hyperpar[1,1]
hist(size$tTotal,breaks=0:ceiling(max(size$tTotal)),xlim=c(0,50))
abline(v=quantile(size$tTotal,q),col="red")
quant<-lapply(1:ncol(s.eff),function(i){
  eta<-s.eff[,i]
  vals<-rgp(n=length(eta),eta=eta,alpha=q,xi=xi)
  #hist(vals,breaks=0:ceiling(max(vals)),xlim=c(0,50))
  abline(v=quantile(unlist(vals),q),col=gray(0,0.15))
})









####################################################################
############## these variable names might not be relevant anymore
### check same for glm using the two methods
mod1 <- lm(tTotal ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size)
mod1 <- glm(Total ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size, family=Gamma(link="log"))
mod1 <- lm(Total^(0.1) ~ VEGZONSNA + WtrUrb_km + logPop_2017 + Road_dens + trees_age + high_name + ISI + FWI, data = size)
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














