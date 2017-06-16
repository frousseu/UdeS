
library(MODIStsp)
library(raster)
library(velox)
library(mgcv)
library(rgeos)
library(sp)
library(rgdal)
library(tmap)
library(data.table)
library(scales)
library(quantreg)
library(FlexParamCurve)
#library(phenex)

#MODIStsp()

## load RData
path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/MODIS/VI_16Days_250m_v6/Time_Series/RData/"
l<-list.files(path)
for(i in seq_along(l)){
  load(paste0(path,l[i]))
  assign(paste0("r",i),raster_ts)
}
rm(raster_ts)

rv<-r5 #r5 is the initial raster used in plot ndvi3
rd<-r3
#v<-velox(r[[1:100]])

ram<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram")
plot(ram)
ram<-gBuffer(ram,width=-0.01)
plot(ram,add=TRUE)

erv_raw<-extract(rv,ram)
erd_raw<-extract(rd,ram)
#ee<-do.call("cbind",lapply(strsplit(v$extract(ram,fun=function(i){paste(i,collapse="_")}),"_"),as.integer))
erv<-erv_raw[[1]]
erd<-erd_raw[[1]]

day<-lapply(1980:2017,function(i){
  seq.Date(as.Date(paste0(i,"-01-01")),as.Date(paste0(i,"-12-31")),by=1)
})
names(day)<-1980:2017
julp<-sapply(strsplit(dimnames(erv)[[2]],"_"),function(i){
  as.integer(i[4]) 
})
datep<-as.Date(sapply(strsplit(dimnames(erv)[[2]],"_"),function(i){
  day[[as.character(i[3])]][as.integer(i[4])] 
}),origin="1970-01-01")
ind<-sapply(strsplit(dimnames(erv)[[2]],"_"),function(i){
  i[2]
})
sat<-substr(sapply(strsplit(dimnames(erv)[[2]],"_"),function(i){
  i[1]
}),1,3)
#erv<-t(apply(erv,1,function(i){identity(i)}))


########################################
### raw images and models

d<-as.data.table(data.frame(x=rep(1:ncol(erv),nrow(erv)),datep=datep,julp=julp,ind=ind,sat=sat,year=as.integer(substr(datep,1,4)),y=as.vector(t(erv)),jul=as.vector(t(erd))),date)
d$y<-d$y/10000
d<-d[!is.na(d$jul) & !is.na(d$y),]
d<-d[,.(y,ind,sat,year,julp,jul,np=.N),by=.(datep)]
d<-d[,.(y,ind,sat,year,julp,np,median=quantile(y,0.5,na.rm=TRUE),mean=mean(y,na.rm=TRUE),n=.N),by=.(datep,jul)]
d<-as.data.frame(d)
d$date<-as.Date(sapply(1:nrow(d),function(i){
  k<-ifelse(d$jul[i]<d$julp[i],1,0)
  day[[as.character(d$year[i]+k)]][d$jul[i]]
}),origin="1970-01-01")
d<-d[d$year>=2003,]
d$datex<-as.numeric(d$date)
d<-d[order(d$datep,d$jul),]



png("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc/ndvi4.png",width=22,height=10,units="in",res=300)
par(mar=c(7,4,4,4))
plot(d$date,d$y,col=gray(0,0.1),xaxt="n",xlab="Date",ylab="NDVI",type="n",xlim=c(min(d$date)-300,max(d$date)))

years<-seq_along(unique(d$year))
invisible(lapply(years,function(i){
  val<-as.numeric(day[[as.character(unique(d$year)[i])]])
  rect(xleft=val-0.5,ybottom=-2000,xright=val+0.5,ytop=12000,col=alpha("grey20",seq_along(val)/length(val)/4),border=NA,xpd=FALSE)
}))
axis.Date(1, at=seq(min(d$date,na.rm=TRUE),max(d$date,na.rm=TRUE),by="1 month"), format="%Y-%m-%d",las=2,cex.axis=0.5)

points(d$date,d$y,col=gray(0,0.025))
d2<-unique(d[,c("date","mean","sat","n","np")])
points(d2$date,d2$mean,col=alpha("green4",d2$n/d2$np),pch=ifelse(d2$sat=="MOD",16,17),cex=1)

sa<-d$sat%in%c("MYD","MOD")
m1<-gam(y~s(datex,k=100),data=d[sa,])
m2<-loess(y~datex,data=d[sa,],span=0.03)
sx<-seq(min(d$datex,na.rm=TRUE),max(d$datex,na.rm=TRUE),by=1)
p1<-predict(m1,data.frame(datex=sx),type="response")
p2<-predict(m2,data.frame(datex=sx),type="response")
lines(sx,p1,col=alpha("blue",0.35),lwd=4)
lines(sx,p2,col=alpha("red",0.35),lwd=4)

peak<-NULL

invisible(peak<-lapply(years[-length(years)],function(i){
  year<-unique(d$year)[i]
  ### up
  dd<-d[d$datep>=paste0(year-1,"-12-01") & d$datep<=paste0(year,"-09-30"),]
  lo1<-list(Asym=0,xmid=12000,scal=0.5,c=0.1)
  up1<-list(Asym=1,xmid=18000,scal=50,c=0.4)
  m1<-nls(y~Asym/(1+exp((xmid-datex)/scal))+c,data=dd,start=list(Asym=0.5,xmid=quantile(dd$datex,0.5,na.rm=TRUE),scal=10,c=0.2),control=list(minFactor=1e-12,maxiter=50),lower=lo1,upper=up1,algorithm="port")
  se<-seq(min(dd$datex),max(dd$datex),by=1) 
  se<-se[-c(1:20,(length(se)-19):length(se))]
  lines(se,predict(m1,data.frame(datex=se)),col=alpha("green4",0.85),lwd=4)
  # bounds
  if(i==1){
    se2<-seq(min(dd$datex)-70,max(dd$datex)+70,by=1)-400
    c2<-mean(d$y[d$jul%in%c(330:365,1:60)])
    Asym2<-mean(d$y[d$jul%in%170:270])-c2
    with(lo1,lines(se2,Asym2/(1+exp((coef(m1)["xmid"]-400-se2)/scal))+c2,col=alpha("green4",0.85),lwd=4))
    with(up1,lines(se2,Asym2/(1+exp((coef(m1)["xmid"]-400-se2)/scal))+c2,col=alpha("green4",0.85),lwd=4))
  }
  
  co<-as.list(coef(m1))
  peak<-c(peak,co$xmid)

  ### optimums
  l<-logistic_optimum(alpha=co$Asym,beta=-co$xmid/co$scal,gamma=1/co$scal)
  l<-unique(unlist(l))
  invisible(lapply(l,function(j){
    lines(rep(j,2),c(-0.3,-0.2),col="green4",lwd=1)
  }))
  
  ### steepness makes it difficult for convergence...
  dd<-d[d$datep>=paste0(year,"-07-01") & d$datep<=paste0(year+1,"-03-01"),]
  lo2<-c(Asym=-1,xmid=12000,scal=0.1,c=co$Asym+co$c)
  up2<-c(Asym=0,xmid=18000,scal=50,c=co$Asym+co$c)
  m2<-tryCatch(nls(y~Asym/(1+exp((xmid-datex)/scal))+c,data=dd,start=list(Asym=-co$Asym,xmid=quantile(dd$datex,0.5,na.rm=TRUE),scal=0.1,c=co$Asym+co$c),control=list(minFactor=1e-12,maxiter=50),lower=lo2,upper=up2,algorithm="port"),error=function(j){TRUE})
  if(!isTRUE(m2)){
    se<-seq(min(dd$datex),max(dd$datex),by=1) 
    se<-se[-c(1:20,(length(se)-19):length(se))]
    lines(se,predict(m2,data.frame(datex=se)),col=alpha("green4",0.85),lwd=4)
  }
  
  

return(peak)}))

legend("topright",title="NDVI",pch=c(1,16,17,NA,NA,NA),lwd=c(NA,NA,NA,4,4,4),col=c(gray(0,0.3),alpha("green4",0.5),alpha("green4",0.5),alpha("blue",0.35),alpha("red",0.35),"green4"),legend=c("Value in a 250m pixel","Moy. Aqua sat.","Moy. Terra sat.","GAM","LOESS","Double logistic"),bty="n",inset=c(0.05,0))

dev.off()


gu<-as.integer(format(as.Date(unlist(peak),origin="1970-01-01"),"%j"))
gu<-gu-mean(gu)
hist(gu)

###
fun<-function(){
  plot(ram,add=TRUE)  
}
plot(r[[1:10]],addfun=fun)


### visualisation prediction (dynamic)
tmap_mode("view")
tm_shape(rv[["MYD13Q1_NDVI_2009_233"]])+tm_raster(alpha=0.9,n=10,palette=rev(terrain.colors(10)))+tm_shape(ram)+tm_borders(lwd=5)+tm_layout(basemaps = c("Esri.WorldImagery","HERE.hybridDay"))




################################################
### Derivatives logistic
# Asym


### alpha beta gamma
logistic<-function(x,alpha=1,beta=1,gamma=1,c=0){
  
  d0<-function(alpha,beta,gamma,c){
    alpha/(1+exp(-beta-gamma*x))+c
  }
  
  d1<-function(alpha,beta,gamma,c){
    alpha*gamma*exp(-beta-gamma*x)*(1+exp(-beta-gamma*x))^(-2)
  }
  
  d2<-function(alpha,beta,gamma,c){
    alpha*gamma^2*exp(-beta-gamma*x)*(exp(-beta-gamma*x)-1)*(1+exp(-beta-gamma*x))^(-3)
  }
  
  d3<-function(alpha,beta,gamma,c){
    alpha*gamma^3*exp(-beta-gamma*x)*(1-4*exp(-beta-gamma*x)+exp(-beta-gamma*x)^2)*(1+exp(-beta-gamma*x))^(-4)
  }
  
  d4<-function(alpha,beta,gamma,c){
    alpha*gamma^4*exp(-beta-gamma*x)*(-1+(11*exp(-beta-gamma*x))-(11*(exp(-beta-gamma*x)^2))+exp(-beta-gamma*x)^3)*(1+exp(-beta-gamma*x))^(-5)
  }
  
  y0<-d0(alpha,beta,gamma,c)
  y1<-d1(alpha,beta,gamma,c)
  y2<-d2(alpha,beta,gamma,c)
  y3<-d3(alpha,beta,gamma,c)
  y4<-d4(alpha,beta,gamma,c)
  
  col<-gray((0:4)/5)
  plot(x,y0,ylim=range(c(y0,y1,y2,y3,y4)),type="n")
  lines(x,y0,lwd=4,col=col[1])
  lines(x,y1,lwd=2,col=col[2])
  lines(x,y2,lwd=2,col=col[3])
  lines(x,y3,lwd=2,col=col[4])
  lines(x,y4,lwd=2,col=col[5])
  
  legend("right",lwd=c(4,2,2,2,2),col=col,legend=c("logistic",paste("derivative",1:4)))
  
}

logistic(seq(-10,20,by=0.01),alpha=2,beta=1/5,gamma=0.8)

logistic_optimum<-function(alpha=1,beta=1,gamma=1,c=0){
  #logisitic function derivative's optimums
  l<-list()
  l[[1]]<--beta/gamma
  l[[2]]<-c(-(log(2+sqrt(3))+beta)/gamma,-(log(2-sqrt(3))+beta)/gamma)  
  l[[3]]<-c(-(log(5+2*sqrt(6))+beta)/gamma,-beta/gamma,-(log(5-2*sqrt(6))+beta)/gamma) 
  l
  
}

l<-logistic_optimum(seq(-10,20,by=0.01),alpha=2,beta=1/5,gamma=0.8)
l<-unique(unlist(l))
lapply(l,function(i){
  lines(rep(i,2),c(-1000,1000),lty=2)
})







