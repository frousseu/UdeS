
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
library(phenex)

#MODIStsp()

## load RData
path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/MODIS/VI_16Days_250m_v6/Time_Series/RData/"
l<-list.files(path)
for(i in seq_along(l)){
  load(paste0(path,l[i]))
  assign(paste0("r",i),raster_ts)
}
rm(raster_ts)

r<-r5
#v<-velox(r[[1:100]])

ram<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram")
plot(ram)
ram<-gBuffer(ram,width=-0.01)
plot(ram,add=TRUE)

ee<-extract(r,ram)
#ee<-do.call("cbind",lapply(strsplit(v$extract(ram,fun=function(i){paste(i,collapse="_")}),"_"),as.integer))
ee<-ee[[1]]
day<-substr(seq.Date(as.Date("2008-01-01"),as.Date("2008-12-31"),by=1),6,10)
xd<-as.Date(sapply(strsplit(dimnames(ee)[[2]],"_"),function(i){
  paste(i[3],day[as.integer(i[4])],sep="-")  
}))
ind<-sapply(strsplit(dimnames(ee)[[2]],"_"),function(i){
  i[2]
})
sat<-substr(sapply(strsplit(dimnames(ee)[[2]],"_"),function(i){
  i[1]
}),1,3)
ee<-t(apply(ee,1,function(i){identity(i)}))


#######################################
### double logistic

d<-as.data.table(data.frame(x=jitter(rep(1:ncol(ee),nrow(ee)),factor=0),xd=xd,ind=ind,sat=sat,year=as.integer(substr(xd,1,4)),y=as.vector(t(ee))))
d<-d[,.(xd,y,ind,sat,year,median=quantile(y,0.5,na.rm=TRUE),mean=mean(y,na.rm=TRUE)),by=.(x)]
d<-d[d$year>=2003,]
d$x<-as.numeric(factor(d$xd))
d$date<-as.numeric(as.factor(substr(d$xd,6,10)))
d$y2<-ifelse(d$y<=0,0,d$y)
d<-as.data.frame(d)

dd<-d[d$xd>="2004-01-01" & d$xd<"2004-12-31",]
plot(dd$x,dd$y,col=gray(0,0.1))
#nls(y2~(L/1+exp(-k*(date-x0))),data=as.data.frame(d)[samp,],start=list(L=6000,k=1,x0=6))
m1<-nls(y2~SSlogis(x,Asym,xmid,scal),data=as.data.frame(dd))
m2<-nlrq(y2~SSlogis(x,Asym,xmid,scal),data=as.data.frame(dd),tau=0.5)
m3<-loess(y2~x,data=as.data.frame(dd),span=0.3)
se<-unique(dd$x)
lines(se,predict(m1,data.frame(x=se)))
lines(se,predict(m2,data.frame(x=se)))
lines(se,predict(m3,data.frame(x=se)))

#ddd<-ddply(dd,.(
#m<-modelNDVI(dd$y,year=2004,method="DLogistic")
#for (i in m){ plot(i) } 


########################################
### raw images and models

d<-as.data.table(data.frame(x=jitter(rep(1:ncol(ee),nrow(ee)),factor=0),xd=xd,ind=ind,sat=sat,year=as.integer(substr(xd,1,4)),y=as.vector(t(ee))))
d$y<-d$y/10000
d<-d[,.(xd,y,ind,sat,year,median=quantile(y,0.5,na.rm=TRUE),mean=mean(y,na.rm=TRUE)),by=.(x)]
d<-d[d$year>=2003,]
d$x<-as.numeric(factor(d$xd))
d<-as.data.frame(d)



png("C:/Users/rouf1703/Documents/ndvi2.png",width=22,height=10,units="in",res=300)
par(mar=c(7,4,4,4))
plot(d$x,d$y,col=gray(0,0.1),xaxt="n",type="n")

years<-seq_along(unique(d$year))
periods<-seq_along((unique(substr(d$xd,6,10))))
invisible(lapply(years,function(i){
   rect(xleft=periods-0.5+(i-1)*max(periods),ybottom=-2000,xright=periods+0.5+(i-1)*max(periods),ytop=12000,col=alpha("grey20",periods/max(periods)/4),border=NA,xpd=FALSE)
}))

axis(1,at=d$x[!duplicated(d$x)][seq(1,length(d$x),by=4)],labels=d$xd[!duplicated(d$x)][seq(1,length(d$x),by=4)],las=2,cex.axis=0.5)
#axis.Date(1, at=unique(d$xd), format="%Y-%m-%d",las=2,cex.axis=0.5)

points(d$x,d$y,col=gray(0,0.03))
d2<-unique(d[,c("x","mean","sat")])
points(d2$x,d2$mean,col=alpha("green4",0.85),pch=ifelse(d2$sat=="MOD",16,17),cex=1)
#lines(d$x,d$median,col="darkgreen",lwd=2)

sa<-d$sat%in%c("MYD","MOD")
m1<-gam(y~s(x,k=100),data=d[sa,])
m2<-loess(y~x,data=d[sa,],span=0.03)
sx<-seq(0,max(d$x),by=1)
p1<-predict(m1,data.frame(x=sx),type="response")
p2<-predict(m2,data.frame(x=sx),type="response")
lines(sx,p1,col=alpha("blue",0.35),lwd=4)
lines(sx,p2,col=alpha("red",0.35),lwd=4)

invisible(lapply(years,function(i){
  year<-unique(d$year)[i]
  ### up
  dd<-d[d$xd>=paste0(year-1,"-12-01") & d$xd<=paste0(year,"-09-30"),]
  m1<-nls(y~Asym/(1+exp((xmid-x)/scal))+c,data=dd,start=list(Asym=0.5,xmid=quantile(dd$x,0.5),scal=1.2,c=0.2),control=list(minFactor=1e-12,maxiter=50),lower=c(0.3,0,0.8,0.1),upper=c(0.9,10000,3,0.3),algorithm="port")
  se<-unique(dd$x)
  lines(se,predict(m1,data.frame(x=se)),col="green4",lwd=4)
  ### down
  #dd<-d[d$xd>=paste0(year,"-07-01") & d$xd<=paste0(year+1,"-02-01"),]
  #m1<-nlrq(y~Asym/(1+exp((xmid-x)/scal))+c,data=dd,start=list(Asym=-0.4,xmid=quantile(dd$x,0.5),scal=1.3,c=0.58),control=list(minFactor=1e-12,maxiter=50))
  #se<-unique(dd$x)
  #lines(se,predict(m1,data.frame(x=se)),col="green4",lwd=4)
  #lines(se,0.35/(1+exp((292-se)/0.8))+0.16,col="black",lwd=8)

}))

legend("topright",title="NDVI",pch=c(16,17,NA,NA,NA),lwd=c(NA,NA,4,4,4),col=c(alpha("green4",0.85),alpha("green4",0.85),alpha("blue",0.35),alpha("red",0.35),"green4"),legend=c("Moy Aqua","Moy Terra","GAM","LOESS","Logisitc"),bty="n",inset=c(0.05,0))

dev.off()


#dev.off()


fun<-function(){
  plot(ram,add=TRUE)  
}
plot(r[[1:10]],addfun=fun)


r2<-projectRa

### visualisation prediction (dynamic)
tmap_mode("view")
tm_shape(r[[24]])+tm_raster(alpha=0.9,n=10,palette=rev(terrain.colors(10)))+tm_shape(ram)+tm_borders(lwd=5)+tm_layout(basemaps = c("Esri.WorldImagery","HERE.hybridDay"))






data(avhrr)
data(modis)


ex<-d$year==2008

ndvi.list1 <- modelNDVI(ndvi.values=d$y[ex]/10000, 
                        year.int=2008, multipleSeasons=FALSE, correction="bise", 
                        method="Growth", MARGIN=2, doParallel=TRUE, slidingperiod=40)

for (ndvi.ob in ndvi.list1){ plot(ndvi.ob) } 







