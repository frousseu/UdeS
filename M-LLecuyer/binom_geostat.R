
library(PrevMap)
library(splancs)
library(geoR)
library(gstat)
library(sp)
library(AICcmodavg)
library(GSIF)
library(scales)
library(ncf)
library(car)
library(data.table)
library(readxl)
library(raster)
library(visreg)
library(PrevMap)
library(bbmle)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(velox)



################################################
### load data
################################################


d<-as.data.frame(read_excel("C:/Users/User/Documents/Lou/LandEcoCorrected_changed.xlsx"),stringsAsFactors=FALSE)
d<-head(d,-1)


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


#########################################
### import raster and build pred grid
#########################################

r1<-"C:/Users/User/Documents/Lou/Mature_forest_2015_include_bajos_secondary.tif"
r2<-"C:/Users/User/Documents/Lou/Mature_forest_2015_not_include_bajos_secondary.tif"
r3<-"C:/Users/User/Documents/Lou/Landcover_2015_extended.tif"
code<-read.table("C:/Users/User/Documents/Lou/Legend.txt",header=TRUE,sep=",",stringsAsFactors=FALSE)
r <- stack(r3)

### build pred grid
ras <- raster(xmn=bbox(r)[1,1],xmx=bbox(r)[1,2],ymn=bbox(r)[2,1],ymx=bbox(r)[2,2],ncols=100,nrows=100)
ras[] <- runif(ncell(ras))
plot(ras)
g<-as(ras,"SpatialPixelsDataFrame")
proj4string(g)<-proj4string(r)
w<-10000
p<-gBuffer(SpatialPoints(coordinates(g)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
text(coordinates(p)[,1],coordinates(p)[,2],1:length(p))

### extract values from grid
v<-velox(r)
e<-v$extract(p,fun=function(x){paste(x,collapse="_")})

### compute proportions
ans<-lapply(e,function(i){
	x<-table(unlist(strsplit(i,"_")))
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


### get proportions at obs
p<-gBuffer(SpatialPoints(coordinates(ds)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
eobs<-v$extract(p,fun=function(x){paste(x,collapse="_")})
ans<-lapply(eobs,function(i){
	x<-table(unlist(strsplit(i,"_")))
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


### newdata and model
model<-Attack~Secondary
newdata<-g@data[,"Secondary",drop=FALSE]


#####################################
### simple glm without RSA
#####################################

glm1<-glm(model,data=d,family=binomial)
pred_glm<-predict(glm1,newdata,type="response")


###################################################################
### explore variogram
###################################################################

coords<-as.matrix(d[,c("X","Y")])
v<-variog(coords=coords,data=resid(glm1),breaks=seq(0,50000,by=500))
fitv<-variofit(v,ini.cov.pars=c(2,5000),cov.model="matern",fix.nugget=FALSE,nugget=0,fix.kappa=TRUE,kappa=0.5)
plot(v)
lines(fitv)


#################################################
### GSIF package
#################################################

dat<-as.data.frame(d)[,-1]
names(dat)[1:2]<-c("x","y")

### model
gsif<-fit.regModel(model,rmatrix=dat,predictionDomain=g,method="GLM",fit.family=binomial(link="logit"),stepwise=FALSE,vgmFun="Exp")

### predictions
pred_gsif<-predict(gsif,g)



###########################################
### use PRevMap for a binomial MCML model
###########################################

# Giorgi et Diggle 2016

# MLE based on MCMC

### assume the nb of draw is one for a bernoulli process ?
d$nbevent<-1

### controls and starting values
control.mcmc<-control.mcmc.MCML(n.sim=10000,burnin=1000,thin=8,h=NULL,c1.h = 0.01,c2.h = 1e-04)
par0<-c(coef(glm1),c(0.99,5700,0.15)) # inputs from the variogram

### first evaluation
mcml1<-binomial.logistic.MCML(model,units.m=~nbevent,par0=par0,coords=~X+Y,data=d,control.mcmc=control.mcmc,kappa=0.5,start.cov.pars=c(par0[4],par0[5]/par0[3]))

### need a second estimation with starting values from the first evaluation (don't understand why)
par1<-coef(mcml1)
mcml1<-binomial.logistic.MCML(model,units.m=~nbevent,par0=par1,coords=~X+Y,data=d,control.mcmc=control.mcmc,kappa=0.5,start.cov.pars=c(par1[4],par1[5]/par1[3]))

### prediction
pred_mcml <- spatial.pred.binomial.MCML(mcml1,coordinates(g),predictors=newdata,control.mcmc = control.mcmc, type = "marginal",scale.predictions = "prevalence",standard.errors = TRUE, thresholds = 0.2,scale.thresholds = "prevalence")

### extra codes
#poly <- coords[chull(coords),]
#grid.pred <- gridpts(poly, xs = 10000, ys = 10000)

#####################################
### use glgm from Bonat & Ribeiro
#####################################

# Bonat & Ribeiro 2015

# MLE based on the Laplace approximation

### source code is their supplementary material
source("C:/Users/User/Documents/GitHub/UdeS/M-LLecuyer/functionssglmm.r")

### start values
init <- start.values.glgm(model, family="binomial", data=d[,c("x","y",colnames(model.frame(model,d)))], coords=d[,c("x","y")],nugget=TRUE, ntrial=1)

### Fitting Binomial SGLMM
glgm1 <- glgm(model, cov.model = "matern", kappa = log(0.5), inits = init, data=d, coords = d[,c("x","y")],nugget=TRUE, family="binomial", ntrial = 1,method.optim = "BFGS", method.integrate = "NR")

### prediction
pred_glgm<-prediction(glgm1,as.data.frame(coordinates(g)))

# I think the prediction function does not allow covariates



#####################################
### compare predictions on a grid
#####################################

g$glm<-pred_glm
g$mcml<-pred_mcml$prevalence$predictions
g$glgm<-inv.logit(pred_glgm)
g$gsif<-pred_gsif@predicted$Attack

gg<-as(g,"SpatialPolygonsDataFrame")

par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(0,0,0,3))
plot(gg,col=gray(1-g$glm),border=gray(1-g$glm));text(par("usr")[1],par("usr")[4]-5000,"glm",xpd=TRUE,adj=c(-1,1),cex=3)
plot(gg,col=gray(1-g$mcml),border=gray(1-g$mcml));text(par("usr")[1],par("usr")[4]-5000,"mcml",xpd=TRUE,adj=c(-1,1),cex=3)
plot(gg,col=gray(1-g$glgm),border=gray(1-g$glgm));text(par("usr")[1],par("usr")[4]-5000,"glgm",xpd=TRUE,adj=c(-1,1),cex=3)
plot(gg,col=gray(1-g$gsif),border=gray(1-g$gsif));text(par("usr")[1],par("usr")[4]-5000,"gsif",xpd=TRUE,adj=c(-1,1),cex=3)
legend("right",col=gray(1-(0:10/10)),legend=round((0:10/10),1),border=NA,pt.cex=2.8,pch=15,bty="n")


#######################################
### Correct away to compute the AIC?
#######################################

# Hoeting et al. 2006

# mle =
# n =
# p =
# k =


sAICc<-function(mle,n,p,k){
  -2*mle+2*n*(p+k+1)/(n-p-k-2)
}

sAICc(as.numeric(fit$log.lik),n=101,p=3,k=2)
aictab(list(m))



#######################################################################
### check convergence of mcml
#######################################################################

par(mfrow=c(3,3))
S.mean <- apply(pred.MCML$samples, 2, mean)
acf(S.mean,main = "")
plot(S.mean,type = "l")
plot(ecdf(S.mean[1:5000]), main = "")
lines(ecdf(S.mean[5001:10000]), col = 2, lty = "dashed")

ind.S <- sample(1:nrow(g), 2)
acf(pred.MCML$samples[ind.S[1],], main = "")
plot(pred.MCML$samples[ind.S[1], ],ylab = paste("Component n.", ind.S[1]), type = "l")
plot(ecdf(pred.MCML$samples[ind.S[1], 1:5000]), main = "")
lines(ecdf(pred.MCML$samples[ind.S[1], 5001:10000]),col = 2, lty = "dashed")















