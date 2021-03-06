
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
library(fields)
library(viridisLite)
library(raster)

### List of things potentially missing:

# scale variables
# turn coordinates to km or something smaller
# better priors for intercept and factors
# check pc priors for range (sd seems to be done)
# figure out how to correct the shift in the posterior predictive checks

cat("\014")

#############
### load data
#load("~/UdeS/Consultation/GStetcher/Doc/LLF_occur.RData")
load("~/UdeS/Consultation/GStetcher/Doc/MSB_occ_nodup.RData")

#####################################################################
### source newdata and toseq

source("https://raw.githubusercontent.com/frousseu/UdeS/master/GStecher/newdata.R")

############################################################################################
### remove everything with NAs (temporary) and take a random sample to reduce computing time
occ<-MSB.occ
#occ<-llf.occur
#occ<-na.omit(occ)
#temp1<-occ[occ$high_name_100m=="Others (deciduous and others)" & occ$PA==1,][1,]
temp1<-occ[occ$VE=="Alpin zon" & occ$PA==1,][1,]
occ<-rbind(occ[sample(1:nrow(occ),1000),],temp1) # sample location to reduce computing time and keep a 1 for Alpin

######################################################
### check the proportions of 0 and 1 in each VEGZONSNA
tab<-table(occ$PA,occ$VEGZONSNA)
tab[2,]/(tab[1,]+tab[2,])

#####################################################################
### convert high_name_100m to factor and log transform the population size
occ$high_name_100m<-as.factor(gsub("\\(|\\)| ","",occ$high_name_100m))
occ$logpop_raster_100m<-log(occ$pop_raster_100m+0.2)

#####################################################################
### change vegzone to complete names
occ$VEGZONSNA<-gsub(" ","_",occ$VEGZONSNA)
occ$VEGZONSNA<-as.factor(occ$VEGZONSNA)

#############################################
### build a spatial object with the locations
occs<-occ
coordinates(occs)<-~Long+Lat
proj4string(occs)<-"+init=epsg:4326"
occsll<-occs
prj<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
occs<-spTransform(occs,CRS(prj))

plot(occs,pch=16,col=alpha(ifelse(occs$PA==0,"blue","red"),0.35))

#########################################
### build a simple glm to compare results
#m1 <- glm (PA ~ VEGZONSNA + Trees_age_100m + WtrUrb_100m + logpop_raster_100m + NSkog_100m + high_name_100m, family = binomial(link = "logit"), data = na.omit(occ))

m1 <- glm (PA ~ NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m + Trees_age_100m + high_name_100m + VEGZONSNA, family = binomial(link = "logit"), data = na.omit(occ))

par(mfrow=c(3,3),mar=c(4,4,3,3))
visreg(m1,scale="response",ylim=0:1)
par(mfrow=c(1,1))

####################################################################
### look at the variogram with the residuals from the previous model
coords <- coordinates(occs)
#v<-variog(coords=coords,data=resid(m1),breaks=seq(0,100000,by=500),max.dist=100000,bin.cloud=TRUE)
#plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b") 

##############
### get sweden
swe <- raster::getData("GADM", country = "SWE", level = 0)
swe<-spTransform(swe,proj4string(occs))

####################################################################
### build a mesh and use sweden as a boundary
# smaller values put in here will make a more precise grid, but will take longer to run

#bound<-inla.nonconvex.hull(coordinates(occs)[sample(1:nrow(occs),100),], concave=-0.5, resolution=20)
mesh<-inla.mesh.2d(loc=coordinates(occs),max.edge=c(10000,30000),offset=c(10000,30000),cutoff=10000,boundary=swe)
plot(mesh,asp=1)

####################################################################
### build spde with pc priors (why use pc priors ?, not sure yet)
spde<-inla.spde2.pcmatern(mesh, # mesh 
                          alpha=2, # smoothness parameter
                          prior.range=c(100000,0.9), # P(practic.range < 100000m) = 0.9
                          prior.sigma=c(3,0.1), # P(sigma > 3) = 0.1,
                          constr=T # forces the spatial field to sum to zero
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
vals<-c(vals,list(intercept=1/5^2,default=0.001)) # not sure if the intercept should be fixed or not in relation to the know prop of 0 or 1s
control.fixed<-list(prec=vals,mean=list(intercept=0,default=0),expand.factor.strategy = "model.matrix")


###########################################################
### build the raster/grid that will be used for predictions
g<-makegrid(swe,n=2000)
g<-SpatialPoints(g,proj4string=CRS(proj4string(occs)))
#o<-over(as(g,"SpatialPolygons"),swe) # makes sure pixels touching are included too, does not change much when the grid gets small
#test1<-st_as_sf(g)
#test2<-st_as_sf(swe)
#test<-st_intersects(test1,test2)
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
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + VEGZONSNA + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + high_name_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + high_name_100m + VEGZONSNA + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + high_name_100m + Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + Firebrk_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + high_name_100m + WtrUrb_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + high_name_100m + Firebrk_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + VEGZONSNA + Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + Firebrk_100m + Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + Firebrk_100m + Trees_age_100m + high_name_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + Firebrk_100m + Trees_age_100m + VEGZONSNA + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Trees_age_100m + VEGZONSNA + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Trees_age_100m + high_name_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Trees_age_100m + high_name_100m + VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m + Trees_age_100m + high_name_100m + VEGZONSNA + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Trees_age_100m * VEGZONSNA + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m * Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + Firebrk_100m * Trees_age_100m + f(spatial,model=spde),
  PA ~ -1 + intercept + NSkog_100m * logpop_raster_100m + Firebrk_100m + Trees_age_100m + f(spatial,model=spde),
  #PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m + Trees_age_100m + high_name_100m * VEGZONSNA + f(spatial,model=spde),
  PA ~ -1 + intercept + NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m * Trees_age_100m + high_name_100m + VEGZONSNA + f(spatial,model=spde)
)

### this is to do some tests to turn spatial models in non-spatial models
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
mm<-mmatrix(modell,occ)


#mm<-lapply(modell,function(i){
#  mod<-as.character(i)
#  mod[3]<-gsub(" \\+ f\\(spatial, model = spde\\)","",mod[3]) # remove spatial effect
#  mod[3]<-gsub("-1 \\+ intercept \\+ ","",mod[3]) # remove spatial effect and intercept notation
#  mod<-as.formula(paste0(mod[c(1,3)],collapse=" "))
#  model.matrix(mod,occ)[,-1]
#})

### this is the model with the dummy variable from the model.matrix as suggested by the pdf E Krainski on using factors
modellmm<-lapply(mm,function(i){
  as.formula(paste("PA ~ -1 + intercept +",paste(dimnames(i)[[2]],collapse=" + "),"+ f(spatial, model = spde)"))
})


###################
### model selection

# this runs every model in the model set

registerDoParallel(min(length(modell),detectCores()-2)) 
getDoParWorkers()

#ml<-vector(mode="list",length=length(modell))

ml<-foreach(i=seq_along(modell),.packages=c("stats","INLA"),.verbose=TRUE) %dopar% {
#for(i in seq_along(modell)){
  ### list of variables
  v<-setdiff(all.vars(modell[[i]]),c("PA","intercept","spatial","spde"))
  ### build the data stack
  spde<-spde # this only to make sure spde is exported to the nodes
  #stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(occ[,v,drop=FALSE])),tag="est") # without model.matrix
  stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),data.frame(mm[[i]])),tag="est") # with model.matrix see avoiding problems with factors in inla.stack by E. Krainski
  #https://groups.google.com/forum/#!searchin/r-inla-discussion-group/factors$20avoiding%7Csort:date/r-inla-discussion-group/iimQ1t1onE0/1fLyJbJcBAAJ
  ### run the model with the eb strategy for faster runs (more approximate)
  m<-inla(modellmm[[i]],data=inla.stack.data(stack.est),control.predictor=list(A=inla.stack.A(stack.est)),family="binomial",control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE,return.marginals=FALSE),control.inla=list(strategy='simplified.laplace',int.strategy="eb"),control.fixed=control.fixed,num.threads=1)
  m
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
#Apn<-inla.spde.make.A(mesh=mesh,loc=matrix(c(NA,NA),ncol=2)[rep(1,n),,drop=FALSE]) # the graphs are build using a random points in the area

################################################
### build newdata with variable values to submit
v<-setdiff(all.vars(bmodel),c("PA","intercept","spatial","spde"))
lp<-newdata(x=occ[,v,drop=FALSE],v=v,n=n,fun=median,list=FALSE)
lp<-lapply(lp,function(i){mmatrix(bmodel,i)})
lpmed<-mmatrix(bmodel,newdata(x=occ[,v,drop=FALSE],v=v,n=1,fun=median,list=FALSE)[[1]][rep(1,length(g)),])[[1]]

########################################################
### bind the data stack for the estimate and for the map
#stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),as.list(occ[,v,drop=FALSE])),tag="est") ### old stack not based on explicit model.matrix
stack.est<-inla.stack(data=list(PA=occ$PA),A=list(A,1),effects=list(c(s.index,list(intercept=1)),data.frame(mm[[b]])),tag="est")
stack.map<-inla.stack(data=list(PA=NA),A=list(Ap,1),effects=list(c(s.index,list(intercept=1)),data.frame(lpmed)),tag="map")
full.stack<-inla.stack(stack.est,stack.map)

#######################################
### add a stack for each focus variable
for(i in seq_along(v)){
  le<-nrow(lp[[v[i]]][[1]])
  if(le!=n){
    AA<-inla.spde.make.A(mesh=mesh,loc=matrix(c(819006,6545844),ncol=2)[rep(1,le),,drop=FALSE]) # for categorical variables
  }else{
    AA<-Apn # for numerical variables
  }
  stack<-inla.stack(data=list(PA=NA),A=list(AA,1),effects=list(c(s.index,list(intercept=1)),data.frame(lp[[v[i]]][[1]])),tag=v[i])     
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
m<-inla(modellmm[[b]],Ntrials=1,data=inla.stack.data(full.stack),control.predictor=list(A=inla.stack.A(full.stack),compute=TRUE,link=1),family="binomial",control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),control.fixed=control.fixed,control.inla=list(strategy='simplified.laplace',int.strategy="eb"))



#####################################
### visualize spatial fields

xlim<-range(coordinates(occs)[,1])
ylim<-range(coordinates(occs)[,2])

proj<-inla.mesh.projector(mesh,xlim=xlim,ylim=ylim,dims=c(300,300))

mfield<-inla.mesh.project(projector=proj,field=m$summary.random[['spatial']][['mean']])
sdfield<-inla.mesh.project(projector=proj,field=m$summary.random[['spatial']][['sd']])

par(mfrow=c(1,2),mar=c(3,3,2,5))

image.plot(list(x=proj$x,y=proj$y,z=mfield),col=viridis(100),asp=1,main="Spatial field (logit scale)") 
axis(1)
axis(2)
plot(swe,add=TRUE,border=gray(0,0.5))
plot(occs,pch=1,cex=3*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.3),add=TRUE)
brks<-c(0.01,0.25,0.50,0.75,0.99)
legend("topleft",pch=1,pt.cex=3*brks,col=gray(0,0.3),legend=brks,bty="n",title="Probability of location\nbeing an actual fire",inset=c(0.02,0.05))

image.plot(list(x=proj$x,y=proj$y,z=sdfield),col=viridis(100),asp=1,main="sd of spatial field (logit scale)") 
axis(1)
axis(2)
plot(swe,add=TRUE,border=gray(0,0.5))
plot(occs,pch=1,cex=3*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.3),add=TRUE)




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
  dat<-data.frame(lp[[v[i]]][[1]])
  if(nrow(p)==n){
    plot(dat[[v[i]]],p[,2],type="l",ylim=c(0,1),xlab=v[i],font=2,ylab="",lty=1,yaxt="n")
    lines(dat[[v[i]]],p[,1],lty=3)
    lines(dat[[v[i]]],p[,3],lty=3)
  }else{
    plot(unique(sort(occ[,v[i]])),p[,2],type="l",ylim=c(0,1),xlab=v[i],font=2,ylab="",lty=1,yaxt="n")
    segments(x0=as.integer(unique(sort(occ[,v[i]]))),x1=as.integer(unique(sort(occ[,v[i]]))),y0=p[,1],y1=p[,3],lty=3)
  }
  axis(2,las=2)
  points(occ[,v[i]],jitter(occ$PA,amount=0.035),pch=16,col=gray(0,0.15))
}
mtext("Probability of being an actual fire",outer=TRUE,cex=1.2,side=2,xpd=TRUE,line=2)


##########################################################
### generate predictions without spatial effet

# page 263 in Zuur

nsims<-100
s<-inla.posterior.sample(nsims,m)
params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  match(i,row.names(s[[1]]$latent))  
}) 
nweights<-grep("spatial",row.names(s[[1]]$latent))

### this is to compare with a glm model
#par(mfrow=c(round(sqrt(2*length(v)),0),ceiling(sqrt(2*length(v)))),mar=c(4,4,3,3),oma=c(0,10,0,0))
#m1 <- glm (PA ~ NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m + Trees_age_100m + high_name_100m + VEGZONSNA, family = binomial(link = "logit"), data = na.omit(occ))
#visreg(m1,scale="response",ylim=0:1)

par(mfrow=c(round(sqrt(length(v)),0),ceiling(sqrt(length(v)))),mar=c(4,3,2,2),oma=c(0,10,0,0))
for(k in seq_along(v)){
  p<-lapply(1:nsims,function(i){
    betas<-s[[i]]$latent[nparams]
    fixed<-cbind(intercept=1,lp[[v[k]]][[1]]) %*% betas
    ### this if we want a spatial part
    #wk<-s[[i]]$latent[nweights]
    #if(is.factor(occ[,v[k]])){
    #  spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(819006,6545844),ncol=2)[rep(1,nlevels(occ[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    #  spatial<-as.matrix(Apn) %*% wk
    #}
    #p<-inla.link.invlogit(fixed+spatial)
    p<-inla.link.invlogit(fixed)
    p
  })
  p<-do.call("cbind",p)
  p<-t(apply(p,1,function(i){c(quantile(i,0.0275),mean(i),quantile(i,0.975))}))
  if(nrow(lp[[v[k]]][[1]])==n){
    vals<-lp[[v[k]]][[1]][,v[k]]
    plot(vals,p[,2],type="l",ylim=c(0,1),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    abline(sum(occ$PA)/nrow(occ),0,col=gray(0,0.2))
    points(occ[,v[k]],jitter(occ$PA,amount=0.035),pch=1,col=gray(0,0.15))
    lines(vals,p[,2],lwd=2)
    lines(vals,p[,1],lty=3)
    lines(vals,p[,3],lty=3)
  }else{
    plot(unique(sort(occ[,v[k]])),p[,2],type="l",ylim=c(0,1),xlab=v[k],font=2,ylab="",lty=1,yaxt="n")
    abline(sum(occ$PA)/nrow(occ),0,col=gray(0,0.2))
    points(jitter(as.integer(occ[,v[k]])),jitter(occ$PA,amount=0.035),pch=1,col=gray(0,0.15))
    segments(x0=as.integer(unique(sort(occ[,v[k]]))),x1=as.integer(unique(sort(occ[,v[k]]))),y0=p[,1],y1=p[,3],lty=3)
  }
  axis(2,las=2)
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
samples<-inla.posterior.sample(100,m)
m$misc$configs$contents
contents<-m$misc$configs$contents
effect<-"APredictor"
id.effect<-which(contents$tag==effect)
ind.effect<-contents$start[id.effect]-1+(1:contents$length[id.effect])[index[["est"]]]
samples.effect<-lapply(samples, function(x) x$latent[ind.effect])
#s.eff<-t(matrix(unlist(samples.effect),byrow=T,nrow=length(samples.effect)))
s.eff<-do.call("cbind",samples.effect)

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
bm[3]<-gsub("-1 \\+ intercept \\+ ","",bm[3]) # remove spatial effect and intercept notation
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
  image.plot(list(x = proj$x, y=proj$y, z = field.proj),xlim = xlim, ylim = ylim, col = viridis(n.col), nlevel=n.col+1,asp=1, ...)
}

local.plot.field(m$summary.random[['spatial']][['mean']],mesh)
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
image(xmean,asp=2,col=viridis(100))
res<-inla.spde2.result(m,"spatial",spde)
plot(res[["marginals.range.nominal"]][[1]], type = "l",main = "Nominal range, posterior density")




################################
### lgcp

# it would be better to do a lgcp for a point process

ds2<-occs[occs$PA==1,]

covList<-with(ds2@data,list(Populati_2=log(Populati_2),NSkog_100m=log(NSkog_100m),MDC=log(MDC)))

ds2$NSkog_100m2<-log(ds2$NSkog_100m)

r <- raster(ncol = 100, nrow = 200, ext = extent(ds2))
r <- rasterize(ds2, r, field = 1, fun = "count", background = 0)
plot(r)

roads<-rasterize(ds2[,"NSkog_100m2"],r, field = "NSkog_100m2", fun = mean, background = 0)
plot(roads)

fit<-lgcp(formula=~NSkog_100m2,
          data=ds2[1:50,],
          grid=20,
          covariates=list(NSkog_100m2=roads),
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
fit <- ppm(Q ~ NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m + Trees_age_100m + high_name_100m + VEGZONSNA, Poisson(),data=locs@data)
fit <- ppm(Q ~ NSkog_100m + logpop_raster_100m + WtrUrb_100m + Firebrk_100m + Trees_age_100m + high_name_100m + VEGZONSNA, Poisson(),data=locs@data)

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
  
  



