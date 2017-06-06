



library("geoR")
library("PrevMap")
data("loaloa")

############################
### PrevMap
############################

loaloa$logit <- log((loaloa$NO_INF + 0.5)/(loaloa$NO_EXAM - loaloa$NO_INF + 0.5))
profile.kappa <- shape.matern(formula = logit ~ 1,coords = ~ LONGITUDE + LATITUDE,data = loaloa, set.kappa = seq(0.2,1.5, length = 15),start.par = c(0.2,0.05), coverage = 0.95)
c(profile.kappa$lower, profile.kappa$upper)
profile.kappa$kappa.hat

coords <- as.matrix(loaloa[, c("LONGITUDE", "LATITUDE")])
vari <- variog(coords = coords, data = loaloa$logit,uvec = c(0, 0.1, 0.15, 0.2, 0.4, 0.8, 1.4, 1.8, 2, 2.5, 3))
vari.fit <- variofit(vari, ini.cov.pars = c(2, 0.2),cov.model = "matern",fix.nugget = FALSE, nugget = 0 ,fix.kappa = TRUE, kappa = 0.5)
par(mfrow = c(1,2))
plot(coords, pch = 20, asp = 1, cex = 0.5, main = "(a)")
plot(vari, main = "(b)")
lines(vari.fit)
vari.fit

fit.glm <- glm(cbind(NO_INF, NO_EXAM - NO_INF) ~ 1, data = loaloa,family = binomial)
par0 <- c(coef(fit.glm), vari.fit$cov.pars, vari.fit$nugget)
c.mcmc <- control.mcmc.MCML(n.sim = 10000, burnin = 2000,thin = 8, h = (1.65)/(nrow(loaloa) ^ (1/6)))
fit.MCML1 <- binomial.logistic.MCML(formula = NO_INF ~ 1,units.m = ~ NO_EXAM, par0 = par0,coords = ~ LONGITUDE + LATITUDE, data = loaloa,control.mcmc = c.mcmc,kappa = 0.5, start.cov.pars = c(par0[3], par0[4]/par0[2]))
fit.MCML1$log.lik

############################
### geostatsp
############################

library("geostatsp")
data("swissRain")
swissRain$lograin = log(swissRain$rain)
fit =  geostatsp:::glgm(formula="lograin", 
                data=swissRain,
                grid=30,
                covariates=swissAltitude, 
                family="gaussian", 
                buffer=20000,
                priorCI=list(sd=c(0.01, 5), 
                range=c(50000,500000),
                sdNugget = c(0.01, 5)), 
                control.mode=list(theta=c(1.6,-0.25,2.9),
                restart=TRUE)
)


cols<-alpha(rev(terrain.colors(100)),1)
cols<-rev(gray(seq(0,1,by=0.01)))
covr<-raster(g[,c("Secondary","Pasture")])
covr<-stack(g[,c("Secondary","Pasture")])
predr<-aggregate(covr,fac=3)
glm1<-glm(Attack~Cat_Typ+Secondary+Pasture,data=d,family="binomial")
#summary(glm1)
#visreg(glm1,scale="response")
p<-predict(glm1,data.frame(as.matrix(predr),Cat_Typ="C"),type="response")
tempr<-predr[[1]]
tempr[]<-p

#plot(as(predr,"SpatialPixels"),col=colo.scale(p,terrain.colors(100)))


m2<-glgm(formula=Attack~Cat_Typ+Secondary+Pasture, 
                        data=ds,
                        grid=predr,
                        covariates=predr, 
                        family="binomial", 
                        buffer=5000,
                        shape=2,
                        priorCI=list(sd=c(u=0.2, alpha=0.05),range=c(2000,50000)),
                        control.compute=list(dic=TRUE,waic=TRUE,mlik=TRUE))

summary(m2$inla)
m$parameters$summary

### check predictions
par(mfrow=c(1,2))
plot(tempr,col=cols)
plot(ds,add=TRUE,pch=1)
plot(m2$raster[["predict.invlogit"]],col=cols)
plot(ds,add=TRUE,pch=1)

### check priors and posteriors
plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density', lwd=2,lty=2)
lines(m$parameters$sd$posterior, lty=1, lwd=2)
legend("topright", lty=2:1, legend=c("prior","posterior"))

plot(m$parameters$range$posterior,type="l",xlim = c(0,50*1000),xlab='range (m)', ylab='density', lwd=2,lty=1)
lines(m$parameters$range$prior, lty=2, lwd=2)
legend("topright", lty=2:1, legend=c("prior","posterior"))

### possible model averaging
ml<-list(m2$inla,m1$inla)
item<-"marginals.fixed"
INLABMA:::fitmargBMA2(ml,item="marginals.fixed",ws=c(1,1))




param=c(shape=2,range=9900,variance=3)
u=seq(0,4,len=20)
uscale = sqrt(8*param['shape'])* u / param['range']
theMaterns = cbind(
  dist=u, 
  manual=	param['variance']*
    ( 1/(gamma(param['shape']) * 2^(param['shape']-1)  ) ) * 
    uscale^param['shape'] * besselK( uscale , param['shape']),
  geostatsp=matern(u, param=param)
)
head(theMaterns)
matplot(theMaterns[,'dist'], 
        theMaterns[,c('manual','geostatsp')],
        col=c('red','blue'), type='l')
legend('topright', fill=c('red','blue'),
       legend=c('manual','geostatsp'))



mymatern = function(u, phi, kappa) {
  uscale = sqrt(8 * kappa) * u/phi
  res = (1/(gamma(kappa) * 2^(kappa - 1))) * uscale^kappa *
    besselK(uscale, kappa)
  res[u == 0] = 1
  res
}

u<-seq(0,15000,by=1)
plot(u,mymatern(u,9000,2),type="l",ylim=c(0,1))







