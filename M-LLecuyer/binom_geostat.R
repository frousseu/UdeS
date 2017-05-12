
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



################################################
### load data
################################################


d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/LandEcoCorrected_changed.xlsx"),stringsAsFactors=FALSE)

ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=16 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
ds<-d
coordinates(ds)<-~X+Y
proj4string(ds) <- ll
ds<-spTransform(ds,CRS(prj))

d$X<-coordinates(ds)[,1]
d$Y<-coordinates(ds)[,2]

d$x<-d$X
d$y<-d$Y


#####################################
### data and grid for predictions
#####################################

sca<-c("_p","_m","_g","tg")
base<-c("Attack","Cat_Typ","Dens_Liv")

### prediction grid
ds$Cat_Typ<-as.factor(ds$Cat_Typ)
e <- extent(coordinates(ds))
r <- raster(e, ncol=50, nrow=50)
vars<-c(names(d)[grep(sca[4],names(d))],"Dens_Liv")
x<-rasterize(ds,r,vars,fun=mean)
plot(x)
g<-as(x,"SpatialPixelsDataFrame")
proj4string(g)<-CRS(proj4string(ds))

### newdata
newdata<-g@data[,"P_indFor_tg",drop=FALSE]


#####################################
### simple glm without RSA
#####################################

glm1<-glm(Attack~P_indFor_tg,data=d,family=binomial)
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
gsif<-fit.regModel(Attack ~ P_indFor_tg,rmatrix=dat,predictionDomain=g,method="GLM",fit.family=binomial(link="logit"),stepwise=FALSE,vgmFun="Exp")

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
mcml1<-binomial.logistic.MCML(formula=Attack~P_indFor_tg,units.m=~nbevent,par0=par0,coords=~X+Y,data=d,control.mcmc=control.mcmc,kappa=0.5,start.cov.pars=c(par0[4],par0[5]/par0[3]))

### need a second estimation with starting values from the first evaluation (don't understand why)
par1<-coef(mcml1)
mcml1<-binomial.logistic.MCML(formula=Attack~P_indFor_tg,units.m=~nbevent,par0=par1,coords=~X+Y,data=d,control.mcmc=control.mcmc,kappa=0.5,start.cov.pars=c(par1[4],par1[5]/par1[3]))

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
source("C:/Users/rouf1703/Documents/UdeS/GitHub/UdeS/M-LLecuyer/functionssglmm.r")

### start values
init <- start.values.glgm(Attack ~ P_indFor_tg, family="binomial", data=d[,c("x","y","Attack","P_indFor_tg")], coords=d[,c("x","y")],nugget=TRUE, ntrial=1)

### Fitting Binomial SGLMM
glgm1 <- glgm(Attack ~ P_indFor_tg, cov.model = "matern", kappa = log(0.5), inits = init, data=d, coords = d[,c("x","y")],nugget=TRUE, family="binomial", ntrial = 1,method.optim = "BFGS", method.integrate = "NR")

### prediction
pred_glgm<-prediction(glgm1,as.data.frame(coordinates(g)))



#####################################
### compare predictions on a grid
#####################################

g$glm<-pred_glm
g$mcml<-pred_mcml$prevalence$predictions
g$glgm<-inv.logit(pred_glgm)
g$gsif<-pred_gsif@predicted$Attack

gg<-as(g,"SpatialPolygonsDataFrame")

par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(0,0,0,3))
plot(gg,col=gray(g$glm),border=gray(g$glm));text(par("usr")[1],par("usr")[4]-5000,"glm",xpd=TRUE,adj=c(0,1))
plot(gg,col=gray(g$mcml),border=gray(g$mcml));text(par("usr")[1],par("usr")[4]-5000,"mcml",xpd=TRUE,adj=c(0,1))
plot(gg,col=gray(g$glgm),border=gray(g$glgm));text(par("usr")[1],par("usr")[4]-5000,"glgm",xpd=TRUE,adj=c(0,1))
plot(gg,col=gray(g$gsif),border=gray(g$gsif));text(par("usr")[1],par("usr")[4]-5000,"gsif",xpd=TRUE,adj=c(0,1))
legend("right",col=gray(0:10/10),legend=round(0:10/10,1),border=NA,pt.cex=2.8,pch=15,bty="n")


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
