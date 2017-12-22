
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
          priorCI=list(sd=c(0.7,2.5),range=c(5000,20000)),
          control.compute=list(waic=TRUE,dic=TRUE,mlik=TRUE)
)

summary(fit$inla)

for(i in 1:10){
  m<-inla.rerun(fit$inla)
  print(c(m$waic$waic,m$dic$dic))
}

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
