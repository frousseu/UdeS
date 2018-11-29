
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
library(sf)
library(DHARMa)
library(ROCR)
library(doParallel)
library(foreach)

### List of things potentially missing:

# scale variables
# turn coordinates to km or something smaller
# better priors for intercept and factors
# check pc priors for range (sd seems to be done)
# figure out how to correct the shift in the posterior predictive checks

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
newdata<-function(x,v=names(x),n=20,fun=mean,list=FALSE){
  
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
load("~/UdeS/Consultation/GStetcher/Doc/LLF_occur.RData")

############################################################################################
### remove everything with NAs (temporary) and take a random sample to reduce computing time
occ<-llf.occur
occ<-na.omit(occ)
occ<-occ[sample(1:nrow(occ),5000),] # sample location to reduce computing time

######################################################
### check the proportions of 0 and 1 in each VEGZONSNA
tab<-table(occ$PA,occ$VEGZONSNA)
tab[2,]/(tab[1,]+tab[2,])

#####################################################################
### convert high_name to factor and log transform the population size
occ$high_name<-as.factor(occ$high_name)
occ$logPop_2017<-log(occ$Pop_2017+0.2)

#####################################################################
### change vegzone to complete names
occ$VEGZONSNA<-gsub(" ","_",occ$VEGZONSNA)
occ$VEGZONSNA<-as.factor(occ$VEGZONSNA)

#############################################
### build a spatial object with the locations
occs<-occ
coordinates(occs)<-~Longitude+Latitude
proj4string(occs)<-"+init=epsg:4326"
prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
occs<-spTransform(occs,CRS(prj))

#########################################
### build a simple glm to compare results
m1 <- glm (PA ~ VEGZONSNA + trees_age + WtrUrb_km + logPop_2017 + Road_dens + high_name, family = binomial(link = "logit"), data = na.omit(occ))
par(mfrow=c(3,3),mar=c(4,4,3,3))
visreg(m1,scale="response")
par(mfrow=c(1,1))

####################################################################
### look at the variogram with the residuals from the previous model
coords <- coordinates(occs)
v<-variog(coords=coords,data=resid(m1),breaks=seq(0,100000,by=500),max.dist=100000,bin.cloud=TRUE)
plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b") 

##############
### get sweden
swe <- raster::getData("GADM", country = "SWE", level = 0)
swe<-spTransform(swe,proj4string(occs))

####################################################################
### build a mesh and use sweden as a boundary
# smaller values put in here will make a more precise grid, but will take longer to run
mesh<-inla.mesh.2d(loc=coordinates(occs),max.edge=c(10000,40000),offset=c(10000,40000),cutoff=10000,boundary=swe)
plot(mesh,asp=1)

####################################################################
### build spde with pc priors (why use pc priors ?, not sure yet)
spde<-inla.spde2.pcmatern(mesh, # mesh 
                          alpha=2, # smoothness parameter
                          prior.range=c(100000,0.9), # P(practic.range < 100000m) = 0.9
                          prior.sigma=c(3,0.1) # P(sigma > 3) = 0.1
                          )


####################################################################
### sensible priors on factors and intercept

fac<-sapply(occ,is.factor)
fac<-names(fac[fac])
fac<-sapply(fac,function(i){levels(occ[,i])})
fac<-unlist(sapply(seq_along(fac),function(i){paste0(names(fac)[i],fac[[i]])}))
vals<-rep(1/5^2,length(fac))
names(vals)<-fac
vals<-as.list(vals)
vals<-c(vals,list(intercept=1/0.001^2,default=0.001))
control.fixed<-list(prec=vals,mean=list(intercept=0.2774,default=0))


###########################################################
### build the raster/grid that will be used for predictions
g<-makegrid(swe,n=30000)
g<-SpatialPoints(g,proj4string=CRS(proj4string(occs)))
#o<-over(as(g,"SpatialPolygons"),swe) # makes sure pixels touching are included too, does not change much when the grid gets small
o<-over(g,swe)
g<-SpatialPixels(g)
g<-g[apply(o,1,function(i){!all(is.na(i))}),]

############################################
### make spatial index to link with the spde
s.index<-inla.spde.make.index(name="spatial",n.spde=spde$n.spde)

###########################
### make observation matrix
A<-inla.spde.make.A(mesh=mesh,loc=coordinates(occs))

#############
### model set

modell<-list(
  PA ~ -1 + intercept + Road_dens + logPop_2017 + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + high_name + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + high_name + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + high_name + trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + Frbreak_km + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + high_name + WtrUrb_km + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + high_name + Frbreak_km + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + VEGZONSNA + trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + Frbreak_km + trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + Frbreak_km + trees_age + high_name + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + Frbreak_km + trees_age + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + trees_age + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + trees_age + high_name + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + trees_age + high_name + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + Frbreak_km + trees_age + high_name + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + trees_age * VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km * trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + Frbreak_km * trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens * logPop_2017 + Frbreak_km + trees_age + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + Frbreak_km + trees_age + high_name * VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + Road_dens + logPop_2017 + WtrUrb_km + Frbreak_km * trees_age + high_name + VEGZONSNA + f(spatial,model=spde)
)

### this is to do some tests to turn spatial models in non-spatial models
#modell<-lapply(modell,function(i){
#  mod<-as.character(i)
#  mod[3]<-gsub(" \\+ f\\(spatial, model = spde\\)","",mod[3]) # remove spatial effect
#  mod<-as.formula(paste0(mod[c(2,1,3)],collapse=""))
#  mod
#})

###################
### model selection

# this runs every model in the model set

registerDoParallel(detectCores()-1) 
getDoParWorkers()

ml<-vector(mode="list",length=length(modell))

ml<-foreach(i=seq_along(modell),.packages=c("stats","INLA"),.verbose=TRUE) %dopar% {
#for(i in seq_along(modell)){
  ### list of variables
  v<-setdiff(all.vars(modell[[i]]),c("PA","intercept","spatial","spde"))
  ### build the data stack
  spde<-spde # this only to make sure spde is exported to the nodes
  stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(occ[,v,drop=FALSE])),tag="est")
  ### run the model with the eb strategy for faster runs (more approximate)
  inla(modell[[i]],data=inla.stack.data(stack.est),control.predictor=list(A=inla.stack.A(stack.est)),family="binomial",control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE,return.marginals=FALSE),control.inla=list(strategy='gaussian',int.strategy="eb"),control.fixed=control.fixed,num.threads=1)
  # print iterations
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
n<-100 # number of divisions in generated values for the focus variable
Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(312180,6342453),ncol=2)[rep(1,n),,drop=FALSE]) # the graphs are build using a random points in the area

################################################
### build newdata with variable values to submit
v<-setdiff(all.vars(bmodel),c("PA","intercept","spatial","spde"))
lp<-newdata(x=occ[,v,drop=FALSE],v=v,n=n,fun=median,list=TRUE)
lpmed<-lapply(newdata(x=occ[,v,drop=FALSE],v=v,n=1,fun=median,list=TRUE)[[1]],function(i){rep(i,length(g))})

########################################################
### bind the data stack for the estimate and for the map
stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(occ[,v,drop=FALSE])),tag="est")
stack.map<-inla.stack(data=list(PA=NA),A=list(Ap,1),effects=list(c(s.index,list(intercept=1)),lpmed),tag="map")
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
  stack<-inla.stack(data=list(PA=NA),A=list(AA,1),effects=list(c(s.index,list(intercept=1)),lp[[v[i]]]),tag=v[i])     
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
m<-inla(bmodel,Ntrials=1,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),family="binomial",control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.fixed=control.fixed,control.inla=list(strategy='gaussian',int.strategy="eb"))


##################################
### build a relative frequency map

### raw predictions on probability of usage
par(mfrow=c(1,4),oma=c(0,5,0,5))
plot(swe,border=gray(0,0.25),lwd=0.01)
points(occs,col=alpha(ifelse(occ$PA==1,"red","blue"),0.4),pch=ifelse(occ$PA==1,16,16),cex=ifelse(occ$PA==1,0.15,0.15))
p<-m$summary.fitted.values[index[["map"]],"mean"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Probability of being an actual fire', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)

### prediction on relative risk
p<-m$summary.fitted.values[index[["map"]],"mean"]
p<-p/(sum(occ$PA)/nrow(occ))
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(brks-1,rev(brewer.pal(11,"RdBu")))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='Relative risk', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)

### sd
p<-m$summary.fitted.values[index[["map"]],"sd"]
gp<-SpatialPixelsDataFrame(g,data=data.frame(p=p))
rgp<-raster(gp)
brks <- seq(min(p),max(p),by=0.01)
cols<-colo.scale(300,rev(brewer.pal(11,"RdYlGn")))
plot(swe)
plot(rgp,col=cols,axes=FALSE,box="n",legend.shrink=1, legend.width=4,add=TRUE,axis.args=list(at=pretty(brks,n=10), labels=pretty(brks,n=10)),
     legend.args=list(text='sd of the probability of being an actual fire', side=4, font=2, line=2.3))
plot(swe,add=TRUE,border=gray(0,0.25),lwd=0.01)


#####################################
### graphical predictions

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
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


######################################################
### model checking and posterior predictive checks

# from bayesian regression modeling with INLA, Faraway frpom google books
# the same thing is shown on Blangardio et Cameletti on p. 168
# something does not work cause response is 1 or 0 and there are NaN in predictions
post.predicted.pval<-vector(mode="numeric",length=nrow(occ))
for(i in 1:nrow(occ)){
  post.predicted.pval[i]<-inla.pmarginal(q=occ$PA[i],marginal=m$marginals.fitted.values[[i]])
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
  rbinom(length(i),size=1,prob=inla.link.invlogit(i))  
})
o<-createDHARMa(simulatedResponse=matprob,observedResponse=occ$PA,fittedPredictedResponse=prob,integerResponse=TRUE)
par(mfrow=c(2,2))
plot(o,quantreg=TRUE)
hist(o$scaledResiduals)

### check same for glm using the two methods

bm<-as.character(modell[[b]])
bm[3]<-gsub(" \\+ f\\(spatial, model = spde\\)","",bm[3]) # remove spatial effect and intercept notation
bm[3]<-gsub("\\-1 \\+ intercept \\+ ","",bm[3]) # remove spatial effect and intercept notation
bm<-as.formula(paste0(bm[c(2,1,3)],collapse=""))
bm

mod <- glm(bm, data = occ, family=binomial)
simulationOutput <- simulateResiduals(fittedModel = mod)
plot(simulationOutput,quantreg=TRUE)
s<-simulateResiduals(mod)
o<-createDHARMa(simulatedResponse = s[["simulatedResponse"]],observedResponse = occ$PA, fittedPredictedResponse = fitted(mod),integerResponse = T)
plot(o,quantreg=TRUE)
hist(o$scaledResiduals)

hist(m$summary.fitted.values[1:nrow(occ),"mean"])
hist(inla.link.invlogit(predict(mod)))

range(m$summary.fitted.values[1:nrow(occ),"mean"])
range(inla.link.invlogit(predict(mod)))

par(ask=TRUE)
plot(m,plot.prior=TRUE)
par(ask=FALSE)


int<-m$marginals.fixed$intercept
plot(int[,1],int[,2],type="l")
hist(inla.link.invlogit(int[,1]))
hist(inla.link.invlogit(seq(-3,3,by=0.01)))


mu<-3
alpha<-0.5
pc.prec<-function(tau,mu,alpha){
  lambda<--log(alpha)/mu
  dtau<-(lambda/2)*(tau^(-3/2))*exp(-lambda*tau^(-1/2))
  sigma<-1/sqrt(dtau)
  sigma
}
v<-seq(0.001,1000,by=0.001)
p<-pc.prec(tau=1/v^2,mu,alpha)
plot(v,p,type="l",ylim=c(0,50),xlim=c(0,5),pos=c(0,0))
p<-1/sqrt(inla.pc.dprec(prec=1/v^2,mu,alpha))
plot(v,p,type="l",ylim=c(0,50),xlim=c(0,5),pos=c(0,0))




res = inla.spde.result(m, "spatial", spde)
par(mfrow=c(2,1))
plot(res$marginals.range.nominal[[1]],
     type="l", main="Posterior density for range")
plot(inla.tmarginal(sqrt, res$marginals.variance.nominal[[1]]),
     type="l", main="Posterior density for std.dev.")
par(mfrow=c(1,1))



### From Baaka BTopic122, spatial field
local.plot.field = function(field, mesh,xlim=c(200000,1000000),ylim=c(6100000,7700000), ...){
  stopifnot(length(field) == mesh$n)
  proj = inla.mesh.projector(mesh, xlim = xlim,ylim = ylim, dims=c(500, 500))
  field.proj = inla.mesh.project(proj, field)
  n.col = 20
  image.plot(list(x = proj$x, y=proj$y, z = field.proj),xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1,asp=1, ...)
}

local.plot.field(m$summary.random[['spatial']][['sd']],mesh)
plot(swe,add=TRUE)
coo<-coordinates(swe)
lwr<-exp(res$summary.log.range.nominal[["0.025quant"]])
mea<-exp(res$summary.log.range.nominal[["mean"]])
upr<-exp(res$summary.log.range.nominal[["0.975quant"]])
arrows(coo[1,1]-0.5*mea, coo[1,2], coo[1,1]+0.5*mea, coo[1,2], length=0.05, angle=90, code=3, lwd=3)
#points(occs,cex=0.1)


##########################################
### confusion matrix with training data

predprob<-m$summary.fitted.values[index[["est"]],"0.5quant"]
cm<-table(as.logical(occ$PA),rbinom(length(predprob),size=1,prob=predprob))
cm

sensitivity<-cm[2,2]/sum(cm[2,]) # true positive rate
specificity<-cm[1,1]/sum(cm[1,]) # true negative rate

##########################################
### confusion matrix on hold-out data


##########################################
### confusion matrix on hold-out data

### code adapted from Myer et al. 2017 (spatiotemporal mosquitoes, https://doi.org/10.1002/ecs2.1854) supplementary material

### but see Boyce et al. 2002 pour le use-availability design

roc.pred<-prediction(predprob,occ$PA)
roc.perf<-performance(roc.pred, measure="tpr", x.measure="fpr") #"tpr" means true positive rate, "fpr" is false positive rate
#tiff(filename="testfig5.tiff",width=3,height=3,units="in",res=300,pointsize=8,compression="lzw")
plot(roc.perf)
abline(a=0,b=1)
#dev.off()
#This code snippet finds the optimal cutoff point
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, roc.pred))

#Obtain the AUC value
auc.perf = performance(roc.pred, measure = "auc")
auc.perf@y.values


##########################################
### cross-validation/cpo/pit measures ?

plot(m) #? PIT and CPO valid with bernoulli response?


#################################
###

#image(inla.mesh.project(mesh,field=m$summary.fitted.values[inla.stack.index(full.stack,tag="latent")$data,"mean"]),dims=c(10,10))
projgrid <- inla.mesh.projector(mesh, dims=c(500,500))
xmean <- inla.mesh.project(projgrid, m$summary.random$spatial$mean)
xsd <- inla.mesh.project(projgrid, m$summary.random$spatial$sd)
image(xmean,asp=2,col=heat.colors(100))
res<-inla.spde2.result(m,"spatial",spde)
plot(res[["marginals.range.nominal"]][[1]], type = "l",main = "Nominal range, posterior density")




################################
### lgcp

# it would be better to do a lgcp for a point process

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

locs<-occs[occs$PA==1,]
locs$x<-coordinates(locs)[,1]
locs$y<-coordinates(locs)[,2]

plot(sweden)
plot(locs,add=TRUE)

o<-owin(xrange=bbox(locs)[1,],yrange=bbox(locs)[2,])
X<-ppp(coordinates(locs)[,1],coordinates(locs)[,2],window=o)

plot(Kest(X))
plot(Kinhom(X))
plot(Ginhom(X))


# inhomogeneous pattern of maples
#X <- unmark(split(lansing)$maple)

Q<-quadscheme(X)

# (1) intensity function estimated by model-fitting
# Fit spatial trend: polynomial in x and y coordinates
fit <- ppm(X, ~ polynom(x,y,8), Poisson(),data=locs@data,covariates=NULL)
fit <- ppm(Q ~ Road_dens + logPop_2017 + WtrUrb_km + Frbreak_km + trees_age + high_name + VEGZONSNA, Poisson(),data=locs@data)
fit <- ppm(Q ~ Road_dens + logPop_2017 + WtrUrb_km + Frbreak_km + trees_age + high_name + VEGZONSNA, Poisson(),data=locs@data)

# (a) predict intensity values at points themselves,
#     obtaining a vector of lambda values
lambda <- predict(fit, locations=X, type="trend")
# inhomogeneous K function
Ki <- Kinhom(X, lambda, nlarge=10000)
Ki <- Ginhom(X, lambda, nlarge=2000)
#Ki <- Kest(X, lambda)
plot(Ki,xlim=c(0,5000))
  
n<-200  
x<-runif(n,0,50)
y<-runif(n,0,50)

x<-sapply(x,function(i){rnorm(10,i,1)})
y<-sapply(y,function(i){rnorm(10,i,1)})

plot(x,y)

o<-owin(xrange=range(x),yrange=range(y))
X<-ppp(x,y,window=o)

plot(envelope(X,Kest))
  
  



