
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
library(spaMM)
library(geoRglm)
library(INLA)



################################################
### load data
################################################


d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/LandEcoCorrected_changed.xlsx"),stringsAsFactors=FALSE)
d<-head(d,-1)
d


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

d$Cat_Typ<-as.factor(d$Cat_Typ)


#########################################
### import raster and build pred grid
#########################################

r1<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Mature_forest_2015_include_bajos_secondary.tif"
r2<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Mature_forest_2015_not_include_bajos_secondary.tif"
r3<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Landcover_2015_extended.tif"
code<-read.table("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Legend.txt",header=TRUE,sep=",",stringsAsFactors=FALSE)
code<-code[code$Code!=0,]# ignore background values in raster
r <- stack(r3)

### build pred grid
ras <- raster(xmn=bbox(r)[1,1],xmx=bbox(r)[1,2],ymn=bbox(r)[2,1],ymx=bbox(r)[2,2],ncols=100,nrows=100)
ras[] <- runif(ncell(ras))
plot(ras)
g<-as(ras,"SpatialPixelsDataFrame")
proj4string(g)<-proj4string(r)
<<<<<<< HEAD
w<-10000
=======
w<-5000
>>>>>>> a8aa53a0963524745df21b98ff3f9f2b04ea264a
p<-gBuffer(SpatialPoints(coordinates(g)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
text(coordinates(p)[,1],coordinates(p)[,2],1:length(p))

### extract values from grid
v<-velox(r)
e<-v$extract(p,fun=function(x){paste(x,collapse="_")})

### compute proportions
ans<-lapply(e,function(i){
	x<-table(unlist(strsplit(i,"_")))
	x<-x[names(x)!=0] # ignore background values in raster
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
g<-g[!is.na(g@data[,1]),] #remove NA values for grid without only background data


### get proportions at obs
p<-gBuffer(SpatialPoints(coordinates(ds)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
eobs<-v$extract(p,fun=function(x){paste(x,collapse="_")})
ans<-lapply(eobs,function(i){
	x<-table(unlist(strsplit(i,"_")))
	x<-x[names(x)!=0] # ignore background values in raster
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
names(d)[which(names(d)=="Agriculture\t")]<-"Agriculture"
names(g)[which(names(g)=="Agriculture\t")]<-"Agriculture"
g$Cat_Typ<-"C"

### newdata and model

g$Forest<-g$Selva_baja+g$Selva_mediana+g$Selva_alta_mediana+g$Subcaducifolia
d$Forest<-d$Selva_baja+d$Selva_mediana+d$Selva_alta_mediana+d$Subcaducifolia
model<-Attack~Cat_Typ+Secondary+Forest+Pasture+Milpa+Bajos
#model<-Attack~Cat_Typ
newdata<-g@data[,attributes(terms(model))$term.labels,drop=FALSE]
#newdata$Cat_TypC<-1
#newdata$Cat_TypS<-0
newdata$Cat_Typ<-factor("C",levels=c("B","C","S"))
#newdata<-newdata[,-match("Cat_Typ",names(newdata))]

#####################################
### simple glm without RSA
#####################################

glm1<-glm(model,data=d,family=binomial)
pred_glm<-predict(glm1,newdata,type="response")
visreg(glm1,scale="response")

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
gsif1<-fit.regModel(model,rmatrix=dat,predictionDomain=g,method="GLM",fit.family=binomial(link="logit"),stepwise=FALSE,vgmFun="Mat")

### predictions
pred_gsif<-predict(gsif1,g)



###########################################
### use PRevMap for a binomial MCML model
###########################################

# Giorgi et Diggle 2016

# MLE based on MCMC

### assume the nb of draw is one for a bernoulli process ?
d$nbevent<-1

### controls and starting values
control.mcmc<-control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8,h=NULL,c1.h = 0.01,c2.h = 1e-04)
par0<-c(coef(glm1),c(0.99,5700,0.15)) # inputs from the variogram
sp<-length(coef(glm1))

### first evaluation
mcml1<-binomial.logistic.MCML(model,units.m=~nbevent,par0=par0,coords=~X+Y,data=d,control.mcmc=control.mcmc,kappa=0.5,start.cov.pars=c(par0[sp+2],par0[sp+3]/par0[sp+1]))

### need a second estimation with starting values from the first evaluation (don't understand why)
par1<-coef(mcml1)
mcml1<-binomial.logistic.MCML(model,units.m=~nbevent,par0=par1,coords=~X+Y,data=d,control.mcmc=control.mcmc,kappa=0.5,start.cov.pars=c(par1[sp+2],par1[sp+3]/par1[sp+1]))

### prediction
pred_mcml <- spatial.pred.binomial.MCML(mcml1,coordinates(g),predictors=newdata,control.mcmc = control.mcmc, type = "marginal",scale.predictions = "prevalence",standard.errors = FALSE, thresholds = NULL,scale.thresholds = NULL)

### extra codes
#poly <- coords[chull(coords),]
#grid.pred <- gridpts(poly, xs = 10000, ys = 10000)

#####################################
### use glgm from Bonat & Ribeiro
#####################################

# Bonat & Ribeiro 2015

# MLE based on the Laplace approximation

### source code is their supplementary material
source("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/functionssglmm.r")

### start values
init <- start.values.glgm(model, family="binomial", data=d[,c("x","y",colnames(model.frame(model,d)))], coords=d[,c("x","y")],nugget=TRUE, ntrial=1)

### Fitting Binomial SGLMM
glgm1 <- glgm(model, cov.model = "matern", kappa = log(0.5), inits = init, data=d[,c(colnames(model.frame(model,d)))], coords = d[,c("x","y")],nugget=TRUE, family="binomial", ntrial = 1,method.optim = "BFGS", method.integrate = "NR")

### prediction
pred_glgm<-prediction(glgm1,as.data.frame(coordinates(g)))

# I think the prediction function does not allow covariates

#####################################
### geoRglm with MCMC
#####################################


### can't run it because of trend error

#geodat<-as.geodata(d,coords.col = 2:3,data.col="Attack",covar.col="Secondary",units.m=rep(1,nrow(d)))
geodat<-as.geodata(d,coords.col = 2:3,data.col="Attack",units.m=rep(1,nrow(d)))
mccontrol<-mcmc.control(S.scale=1, S.start="random", burn.in=2000, thin=3, n.iter=10000, phi.start=1,phi.scale=1)
kcontrol<-krige.glm.control(type.krige = "sk",obj.model=NULL,cov.model="exponential",cov.pars=c(3,5000),beta=c(1,1),kappa=0.5,nugget=0.15)
trend <- trend.spatial(~1)
ocontrol<-output.glm.control(sim.posterior=FALSE, sim.predict=FALSE, keep.mcmc.sim=FALSE, quantile=FALSE,inference=TRUE,messages=TRUE)

geoRglm1<-binom.krige(geodat,units.m = "default", locations = coordinates(g)[1:500,],mcmc.input=mccontrol, krige=kcontrol, output=ocontrol)

geomodel<-list(family="binomial",cov.pars=c(3,5000),beta=c(1,0.5,0.5,1,1,1),trend="2nd",cov.model="matern",nugget=0.3)
m<-glsm.mcmc(geodat,coords=geodat$coords,data=geodat$data,units.m="default",model=geomodel,mcmc.input=mccontrol,messages=TRUE)

p<-glsm.krige(m,locations=coordinates(g)[1:10,],micro.scale=NULL,trend.l="2nd")

#test<- create.mcmc.coda(m, mcmc.input = list(thin = 1))
#autocorr.plot(test)

spa<-corrHLfit(Attack~Forest+Secondary+Pasture+Cat_Typ+Matern(1|x+y),data=d,family=binomial)
#spa<-corrHLfit(Attack~Cat_Typ+Matern(1|x+y),data=d,family=binomial)
pspa<-predict(spa,newdata=cbind(newdata,g@data,as.data.frame(coordinates(g))))
#filled.mapMM(spa)

form <- Attack ~ Secondary
i<-inla(model,data=d,family="binomial",Ntrials=1)


#####################################
### compare predictions on a grid
#####################################

g$glm<-pred_glm
g$mcml<-pred_mcml$prevalence$predictions
g$glgm<-inv.logit(pred_glgm)
g$gsif<-pred_gsif@predicted$Attack
g$spa<-pspa
  
gg<-as(g,"SpatialPolygonsDataFrame")

par(mfrow=c(2,3),mar=c(0,0,0,0),oma=c(0,0,0,3))
plot(gg,col=gray(1-g$glm),border=gray(1-g$glm));text(par("usr")[1],par("usr")[4]-5000,"glm",xpd=TRUE,adj=c(-1,1),cex=3)
points(d$X,d$Y,pch=1,cex=1,col=ifelse(d$Attack==0,"white","black"))
plot(gg,col=gray(1-g$mcml),border=gray(1-g$mcml));text(par("usr")[1],par("usr")[4]-5000,"mcml",xpd=TRUE,adj=c(-1,1),cex=3)
plot(gg,col=gray(1-g$glgm),border=gray(1-g$glgm));text(par("usr")[1],par("usr")[4]-5000,"glgm",xpd=TRUE,adj=c(-1,1),cex=3)
plot(gg,col=gray(1-g$gsif),border=gray(1-g$gsif));text(par("usr")[1],par("usr")[4]-5000,"gsif",xpd=TRUE,adj=c(-1,1),cex=3)
plot(gg,col=gray(1-g$spa),border=gray(1-g$spa));text(par("usr")[1],par("usr")[4]-5000,"spaMM",xpd=TRUE,adj=c(-1,1),cex=3)
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

sAICc(as.numeric(mcml1$log.lik),n=101,p=3,k=6)
aictab(list(glm1))



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















