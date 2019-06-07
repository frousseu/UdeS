
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
library(fitdistrplus)
library(actuar)

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
size<-size[sample(1:nrow(size),2000),]

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
lambda<-3.25 # pc (like?) prior for xi parameter of gp (see below for a graph of prior)
hyper.gp <- list(theta = list(prior = "loggamma",param = c(1,lambda))) # c(1,15) is supposed to be the default prior. See below

ml<-foreach(i=seq_along(modell),.packages=c("stats","INLA"),.verbose=TRUE) %dopar% {
  # list of variables
  v<-setdiff(all.vars(modell[[i]]),c("tTotal","intercept","spatial","spde"))
  # build the data stack
  spde<-spde # this only to make sure spde is exported to the nodes
  stack.est<-inla.stack(data=list(tTotal=size$tTotal),A=list(A,1),effects=list(c(s.index,list(intercept=1)),data.frame(mm[[i]])),tag="est")
  # run the model with the eb strategy for faster runs (more approximate)
  m<-inla(modellmm[[i]],data=inla.stack.data(stack.est),control.predictor=list(A=inla.stack.A(stack.est)),control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE,return.marginals=FALSE),control.inla=list(strategy='simplified.laplace',int.strategy="eb"),family="gp",control.family=list(list(control.link=list(quantile=q),hyper=hyper.gp)),control.fixed=control.fixed,num.threads=1)
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


### from haakon bakka, BTopic112
nsims<-200
samples<-inla.posterior.sample(nsims,m)
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor" # not sure if should use APredictor or Predictor
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-do.call("cbind",samples.effect)
xi.eff<-sapply(samples, function(x) x$hyperpar[grep("genPareto",names(x$hyperpar))])
#unique(sapply(strsplit(rownames(m$summary.fitted.values),"\\."),function(i){paste(i[1:min(2,length(i))],collapse=" ")}))

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

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  match(i,row.names(samples[[1]]$latent))  
}) 
nweights<-grep("spatial",row.names(samples[[1]]$latent))

### this is to compare with a quantile model
#par(mfrow=c(round(sqrt(2*length(v)),0),ceiling(sqrt(2*length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
#mq <- rq (tTotal ~ Road1k + Pp_1000 + urbwtr1k + Ag_1000 + h__1000 + FWI, data = na.omit(size),tau=q)
#na<-all.vars(mq$formula[[3]])
#grobs<-lapply(na,function(i){
#  if(is.factor(size[,i])){
#    plot(ggpredict(mq,terms=i),limits=c(0,100),raw=TRUE)  
#  }else{
#    plot(ggpredict(mq,terms=paste(i,"[n=50]")),limits=c(0,100),raw=TRUE) 
#  }
#})
#grid.arrange(grobs=grobs,ncol=3)

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,3,2,2),oma=c(0,10,0,0))
for(k in seq_along(v)){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    fixed<-cbind(intercept=1,lp[[v[k]]][[1]]) %*% betas
    ### this if we want a spatial part
    #wk<-samples[[i]]$latent[nweights]
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

#lambda<-6.5 # lambda 3.25 for a 10% of chance of xi being > 0.5, 6.5 for a 1%
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

### check with inla model
# xi should be sampled as well and not fixed to its mean value
fitted<-m$summary.fitted.values[index[["est"]],"mean"]
matprob<-do.call("cbind",lapply(1:ncol(s.eff),function(i){
  rgp(n=nrow(s.eff),eta=s.eff[,i],alpha=q,xi=xi.eff[i]) #m$summary.hyperpar[1,1]) # this picks a fire size from the genPareto with estimated xi for each obs
}))
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
  #abline(v=mean(exp(s.eff[,i])),col=alpha("blue",0.1))
  vals<-rgp(n=nrow(s.eff),eta=s.eff[,i],alpha=q,xi=xi.eff[i])#m$summary.hyperpar[1,1])
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
gppdf<-function(y,sigma,xi){(1/sigma)*(1+xi*(y/sigma))^(-1*((1/xi)+1))} # PDF

nsamp<-100
samp<-sample(1:ncol(s.eff),nsamp,replace=TRUE)
#eta<-m$summary.linear.predictor[index[["est"]],"mean"][sample(1:nrow(size),nsamp,replace=TRUE)]
eta<-sapply(samp,function(i){sample(s.eff[,i],1)})
#xi<-m$summary.hyperpar[1,1]
xi<-xi.eff[samp]
#hist(eta)
qs<-quantile(size$tTotal,0.99)
qs2<-qs*1.5
va<-seq(0,qs2,by=0.01)
brks<-0:ceiling(max(size$tTotal))
hist(size$tTotal,breaks=brks,xlim=c(0,qs2),freq=FALSE,border="white",col="tomato")
quant<-lapply(1:nsamp,function(i){
  #xi<-m$summary.hyperpar[1,1]
  sigma<-(xi[i]*exp(eta[i]))/((1-q)^(-xi[i])-1)
  lines(va,gppdf(y=va,sigma=sigma,xi=xi[i]),col=gray(0,0.1))
  #integrate(gppdf,lower=0,upper=30,sigma=sigma,xi=xi)$value
  #abline(v=exp(eta[i]),col=gray(0,0.2),lwd=1)
  abline(v=quantile(rgp(n=nrow(size),eta=eta[i],alpha=q,xi=xi[i]),q),col=gray(0,0.2),lwd=1)
})
abline(v=qs,col="tomato",lwd=5)
#f<-fitdist(size$tTotal,distr="pareto")#,start=list(shape=sigma,scale=xi))
#lines(va,dpareto(va,shape=f$estimate[["shape"]],scale=f$estimate[["scale"]]),col="darkgreen",lwd=2)
#h<-hist(size$tTotal,breaks=c(seq(0,par("usr")[2],length.out=50),max(size$tTotal)),plot=FALSE)
#points(h$mids,rescale(h$density,c(0,1.5)),lwd=2,col="darkgreen")
#hist(rgp(n=10000,eta=median(eta),alpha=q,xi=xi),breaks=100)
#hist(unlist(quant)) ### histograms of % fires below 30 ha (upper in integrate)
#hist(size$tTotal,breaks=seq(0,max(size$tTotal),by=0.#),xlim=c(0,2))


#####

xi<-m$summary.hyperpar[1,1]
hist(size$tTotal,breaks=0:ceiling(max(size$tTotal)),xlim=c(0,50))
abline(v=quantile(size$tTotal,q),col="red")
quant<-lapply(1:ncol(s.eff),function(i){
  eta<-s.eff[,i]
  vals<-rgp(n=length(eta),eta=eta,alpha=q,xi=xi.eff[i])
  #hist(vals,breaks=0:ceiling(max(vals)),xlim=c(0,50))
  abline(v=quantile(unlist(vals),q),col=gray(0,0.15))
})



#################################################
### simulations with known GPD from inla.docs
#################################################

### and comparison with quantile regression

rgp <- function(n, sigma, eta, alpha, xi = 0.001)
{
  if (missing(sigma)) {
    stopifnot(!missing(eta) && !missing(alpha))
    sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
  }
  return (sigma / xi * (runif(n)^(-xi) -1.0))
}
n = 1000
n2 = 100
x = runif(n)-0.5
eta = 1+x
alpha = 0.99
xi = 0.6
y = rgp(n, eta = eta, alpha = alpha, xi=xi)
#y = rgpd(n, mu = 
x2 = seq(min(x),max(x),length.out=n2)
x = c(x,x2)
y = c(y,rep(NA,n2))
d<-data.frame(y,x)
r = inla(y ~ 1+x,data = data.frame(y, x),
         family = "gp",
         control.family = list(control.link = list(quantile = alpha)),
         control.predictor = list(compute=TRUE),control.compute=list(config = TRUE),
         verbose=FALSE)


par(mfrow=c(2,2))
### 1
plot(d$x[1:n],d$y[1:n])
lines(x2,exp(r$summary.fitted.values$mean[(n+1):(n+n2)]))
m<-rq(y~x,tau=alpha,data=d[1:n,])
p<-predict(m,data.frame(x=x2))
lines(x2,p,lty=3)
legend("topleft",lty=1:2,legend=c("INLA GPD","Quantile Regression"),title=paste("Quantile =",alpha),bty="n")
### 2
samples<-inla.posterior.sample(100,r)
r$misc$configs$contents
contents<-r$misc$configs$contents
effect<-"Predictor" # not sure if should use APredictor or Predictor
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[1:n]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
s.eff<-do.call("cbind",samples.effect)
plot(density(log(d$y[1:n])),ylim=range(density(s.eff[,1])$y))
abline(v=quantile(log(d$y[1:n]),alpha),lwd=2)
abline(v=log(quantile(d$y[1:n],alpha)),lwd=2)
invisible(lapply(1:ncol(s.eff),function(i){
  lines(density(s.eff[,i]),col=alpha("red",0.1))    
}))
### 3
b<-0.5
brks<-seq(0,ceiling(max(d$y[1:n])),by=b)
xlim<-range(d$y[1:n])
h<-hist(d$y[1:n],breaks=brks,xlim=xlim,col="grey70",border="white",ylab="Counts")
vals<-lapply(1:ncol(s.eff),function(i){
  #abline(v=mean(exp(s.eff[,i])),col=alpha("blue",0.1)) # not sure what this represents
  vals<-rgp(n=nrow(s.eff),eta=s.eff[,i],alpha=alpha,xi=r$summary.hyperpar[1,1])
  hi<-hist(vals,breaks=seq(0,ceiling(max(vals)),by=b),plot=FALSE)
  points(brks[-1]-diff(brks)/2,hi$counts[1:(length(brks)-1)],col=alpha("red",0.05),pch=16)
  abline(v=quantile(vals,alpha),col=alpha("red",0.15))
})
abline(v=quantile(d$y[1:n],alpha),lwd=2)

f<-fitdist(y[1:n],distr="pareto")#,start=list(shape=sigma,scale=xi))
x<-seq(0,max(y[1:n]),length=100)
lines(x,dpareto(x,shape=f$estimate[["shape"]],scale=f$estimate[["scale"]]))


### 4
gppdf<-function(y,sigma,xi){(1/sigma)*(1+xi*(y/sigma))^(-1*((1/xi)+1))} # CDF
samp<-100
sigma<-(xi*exp(mean(eta)))/((1-alpha)^(-xi)-1)
sigma1<-(xi*exp(min(eta)))/((1-alpha)^(-xi)-1)
sigma2<-(xi*exp(max(eta)))/((1-alpha)^(-xi)-1)
va<-seq(0,5,by=0.01)
plot(va,gppdf(y=va,sigma=sigma,xi=xi),type="l",lwd=2,col="black")
lines(va,gppdf(y=va,sigma=sigma1,xi=xi),type="l",lwd=1,col="black")
lines(va,gppdf(y=va,sigma=sigma2,xi=xi),type="l",lwd=1,col="black")
quant<-lapply(sample(r$summary.linear.predictor[1:1000,"mean"],samp),function(i){
  xi2<-r$summary.hyperpar[1,1]
  sigma<-(xi2*exp(i))/((1-alpha)^(-xi2)-1)
  lines(va,gppdf(y=va,sigma=sigma,xi=xi2),col=alpha("red",0.1))
  integrate(gppdf,lower=0,upper=30,sigma=sigma,xi=xi2)$value
})




###########################################
#### fit distributions

sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
vals<-rpareto(n=1000,shape=3,scale=0.5)
#vals<-rgamma(n=1000,shape=3,scale=0.5)
vals<-y[1:n]
f<-fitdist(vals,distr="pareto")#,start=list(shape=sigma,scale=xi))
hist(vals,breaks=50,freq=FALSE)
x<-seq(0,max(vals),length=100)
lines(x,dpareto(x,shape=f$estimate[["shape"]],scale=f$estimate[["scale"]]))

vals<-size$tTotal
f<-fitdist(vals,distr="pareto")#,start=list(shape=sigma,scale=xi))
hist(vals,breaks=1000,freq=FALSE,xlim=c(0,quantile(vals,0.99)))
x<-seq(0,max(vals),length=10000)
lines(x,dpareto(x,shape=f$estimate[["shape"]],scale=f$estimate[["scale"]]))







             