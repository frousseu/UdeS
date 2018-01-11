
#ds2<-ds[,c("Attack","Cat_Typ","Bridge_p","Branch_p")]
#names(ds2)<-c("y","f","x1","x2")

#writeOGR(ds2,dsn="C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="ds",driver="ESRI Shapefile",overwrite=TRUE)


#

library(geostatsp)
library(INLA)
library(sp)
library(rgdal)

ds<-readOGR(dsn="C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="ds")



fit<-glgm(y~f+x1+x2,
          data=ds,
          grid=20,
          covariates=NULL, 
          family="binomial", 
          buffer=10000,
          shape=1,
          priorCI=list(sd=c(0.5,3),range=c(2000,40000)),
          control.compute=list(waic=TRUE,dic=TRUE,mlik=TRUE),
          control.fixed=list(mean=list(fC=0,fS=0,x1=0,x2=0),prec=list(fC=1/((5)^2),fS=1/((5)^2),x1=1/((0.5)^2),x2=1/((0.5)^2))),
          num.threads=1
)

#summary(fit$inla)

#fit$inla$all.hyper$fixed

#autoplot(fit$inla,which=1,priors=TRUE)
plot_fixed_marginals(fit$inla, priors = TRUE, CI = TRUE)


x<-seq(-5,5,by=0.01)
plot(x,dnorm(x,0,0.5),type="l")
plot(x,dnorm(x,0,0.5),type="l")
lines(fit$inla$marginals.fixed$x2[,1],fit$inla$marginals.fixed$x2[,2])


for(i in 1:3){
  m<-inla.rerun(fit$inla)
  print(c(m$waic$waic,m$dic$dic))
}

summary(glm(ml[[16]][[3]],data=d,family="binomial"))

m<-fit

par(mfrow=c(1,2))

ma<-max(c(m$parameters$sd$prior[,2],m$parameters$sd$posterior[,2]))
plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density',lty=2,xlim=c(0,5),ylim=c(0,ma))
lines(m$parameters$sd$posterior,lty=1)
#lines(seq(0,10,by=0.1),dlgamma(seq(0,10,by=0.1),m$parameters$sd$params.intern$param[1],m$parameters$sd$params.intern$param[2]),col="blue")
legend("topright", lty=2:1, legend=c("prior","posterior"))

ma<-max(c(m$parameters$range$prior[,2],m$parameters$range$posterior[,2]))
plot(m$parameters$range$posterior,type="l",xlim = c(0,50*1000),xlab='range (m)', ylab='density',lty=1,ylim=c(0,ma))
lines(m$parameters$range$prior,lty=2)
#lines(dgamma(seq(0,50000,by=10),m$parameters$range$params.intern[1],m$parameters$range$params.intern[2]),col="red")
legend("topright", lty=2:1, legend=c("prior","posterior"))
