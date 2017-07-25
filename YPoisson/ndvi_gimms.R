

library(foreach)
library(doParallel)
library(signal)
library(robustbase)
library(MODIS)
library(plyr)
library(gimms)
library(ncdf4)
library(rgdal)
library(rasterVis)
library(rgeos)
library(FRutils)
library(signal)
library(scales)
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
library(signal)
library(zoo)
library(rasterVis)
library(FRutils)
library(phenex)

# y = a year as an integer
# j = the julian day to convert
# j2 = the julian day that marks the beginning of a block

###logisitic function derivative's optimums
logistic_optimum<-function(alpha=1,beta=1,gamma=1,offset=0){
  l<-list()
  l[[1]]<-as.list(-beta/gamma)
  l[[2]]<-as.list(data.frame(t(cbind(-(log(2+sqrt(3))+beta)/gamma,-(log(2-sqrt(3))+beta)/gamma))))  
  l[[3]]<-as.list(data.frame(t(cbind(-(log(5+2*sqrt(6))+beta)/gamma,-beta/gamma,-(log(5-2*sqrt(6))+beta)/gamma))))
  l
}

### alpha beta gamma
logistic<-function(x,alpha=1,beta=1,gamma=1,offset=0){
  d0<-function(alpha,beta,gamma,offset){
    alpha/(1+exp(-beta-gamma*x))+offset
  }
  d1<-function(alpha,beta,gamma,offset){
    alpha*gamma*exp(-beta-gamma*x)*(1+exp(-beta-gamma*x))^(-2)
  }
  d2<-function(alpha,beta,gamma,offset){
    alpha*gamma^2*exp(-beta-gamma*x)*(exp(-beta-gamma*x)-1)*(1+exp(-beta-gamma*x))^(-3)
  }
  #d3<-function(alpha,beta,gamma,c){
  #  alpha*gamma^3*exp(-beta-gamma*x)*(1-4*exp(-beta-gamma*x)+exp(-beta-gamma*x)^2)*(1+exp(-beta-gamma*x))^(-4)
  #}
  #d4<-function(alpha,beta,gamma,c){
  #  alpha*gamma^4*exp(-beta-gamma*x)*(-1+(11*exp(-beta-gamma*x))-(11*(exp(-beta-gamma*x)^2))+exp(-beta-gamma*x)^3)*(1+exp(-beta-gamma*x))^(-5)
  #}
  y0<-d0(alpha,beta,gamma,offset)
  y1<-d1(alpha,beta,gamma,offset)
  y2<-d2(alpha,beta,gamma,offset)
  #y3<-d3(alpha,beta,gamma,c)
  #y4<-d4(alpha,beta,gamma,c)
  list(unname(y0),unname(y1),unname(y2))
}


jul2nd<-function(y,j,j2=NULL){
  inc<-0
  if(!is.null(j2)){
    inc<-ifelse(j<j2,1,0) 
  }
  as.integer(as.Date(paste0(y+inc,"-01-01")))+j-1
}

findminmax<-function(x,n=1,beg="06-01",end="11-01",max=TRUE){
  stopifnot(!is.null(names(x)))
  d<-substr(names(x),6,10)
  bloc<-d>=beg & d<=end
  run<-rle(bloc)
  l<-Map(":",c(1,head(cumsum(run[[1]]),-1))[run[[2]]],cumsum(run[[1]])[run[[2]]])
  res<-lapply(l,function(i){
    r<-base:::rank(ifelse(max,-1,1)*x[i])
    val<-sort(r)[1:n]
    index<-i[match(val,r)]
    index   
  })
  res
}

fitLog<-function(x,mmdate=c("12-01","09-15")){
  years<-as.integer(unique(substr(names(x),1,4)))
  years<-years[years<=2016]
  l<-lapply(years,function(i){
    paste(c(i-1,i),mmdate,sep="-")  
  })
  peak<-lapply(l,function(i){
    sx<-x[which(names(x)>=i[1] & names(x)<=i[2])]
    d<-data.frame(y=sx,x=as.integer(as.Date(names(sx))))
    min_xmid<-as.integer(as.Date(i[1])+120)
    max_xmid<-as.integer(as.Date(i[2])-75)
    lo1<-list(Asym=0,xmid=min_xmid,scal=10,offset=0.0)
    up1<-list(Asym=1,xmid=max_xmid,scal=40,offset=0.8)
    m1<-tryCatch(nls(y~Asym/(1+exp((xmid-x)/scal))+offset,data=d,start=list(Asym=0.5,xmid=min_xmid+(max_xmid-min_xmid)/2,scal=20,offset=0.2),control=list(minFactor=1e-12,maxiter=500),lower=lo1,upper=up1,algorithm="port"),error=function(j){TRUE})
    if(!isTRUE(m1)){
      se<-seq(min(d$x),max(d$x),by=1) 
      coef(m1)
    }else{
      NA  
    }
  })
  peak
}

#############################################################
##### Get regions of interest ###############################

# A list of SpatialPolygonsDataFrame for which a metric is wanted

### Refuges Alberta BC
z<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/YPoisson/Doc",layer="largezone_BC_Alberta")
buff<-gBuffer(bbox2pol(z),width=50000)
pol<-spTransform(z,CRS("+init=epsg:4326")) # en latlon
buff<-spTransform(buff,CRS("+init=epsg:4326")) # en latlon

#############################################################
##### Get GIMMS raster ######################################

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/GIMMS/"

x<-list.files(path)
ts<-monthlyIndices(x, version = 1, timestamp = T)
doy<-as.Date(as.character(ts+round(c(diff(ts),17)/2,0)))
r<-rasterizeGimms(x=paste0(path,x),ext=buff,cores=6L)
#names(r)<-as.character(doy)
rd<-stack(setValues(r,as.integer(format(rep(doy,each=ncell(r)),"%j"))))

#tmap_mode("view")
#tm_shape(pol)+tm_borders(lwd=5,alpha=0.3,col="red")+tm_shape(r[[1]])+tm_raster(palette=rev(terrain.colors(100)))+tm_layout(basemaps=c("Esri.WorldImagery","Esri.WorldShadedRelief","Esri.NatGeoWorldMap"))

#############################################################
##### Extract raster values for each region #################

#registerDoParallel(6) 
#getDoParWorkers()

# build a list of matrices of each pixel ndvi value for each polygon
#v<-foreach(i=1:length(pol),.packages=c("raster","sp")) %dopar% {
#r<-aggregate(r,fac=2) # only to speed things up
v<-extract(r)

#vd<-foreach(i=1:length(pol),.packages=c("raster","sp")) %dopar% {
#rd<-aggregate(rd,fac=2) # only to speed things up
vd<-extract(rd)

#############################################################
##### Compute metrics #######################################

### les variables seront:
  
#NDVImax_log: maximum ndvi according to logistic curve
#NDVImaxtime_log: timing of the maximum ndvi according to logistic curve (does not make much sense, unless a double-l is used
#NDVIgu_log: slope at maximum green-up
#NDVIgutime_log: timing of the maximum green-up
#NDVIgutimebeg_log: timing of the beginning of the green-up, defined by the maximum of the second derivative
#NDVIgutimeend_log: timing of the end of the green-up, defined by the minimum of the second derivative

#NDVImax_sg: same thing, but with the Savitsky-Golay filters
#NDVImaxtime_sg
#NDVIgu_sg
#NDVIgutime_sg

cell<-lapply(1:nrow(v[1:nrow(v),]),function(i){
    pl<-FALSE
    val<-v[i,]
    if(all(is.na(val))){
      val[seq_along(val)]<-1  
    }
    jul<-vd[i,]
    names(val)<-doy
    names(jul)<-doy
    
    ### Savitsky-Golay Filter
    s0<-sgolayfilt(na.approx(val),n=7,p=3,m=0) ### hacky take out of first year that begins in 07-08, thus only max available, not gu
    s1<-sgolayfilt(na.approx(val),n=15,p=3,m=1)
    names(s0)<-names(val)
    names(s1)<-names(val)
    #pos<-findminmax(s1,n=5,beg="03-01",end="07-01")
    #sg<-as.Date(sapply(pos,function(k){mean(as.Date(names(s1)[k]))}))
    pos0<-unlist(findminmax(s0,n=1,beg="03-01",end="10-01"))
    pos1<-unlist(findminmax(s1,n=1,beg="03-01",end="07-01"))
    
    NDVImax_sg<-s0[pos0]
    NDVImaxtime_sg<-as.Date(names(s0)[pos0])
    NDVIgu_sg<-s1[pos1]
    NDVIgutime_sg<-as.Date(names(s1)[pos1])
    
    ### Logistic curve
    mLog<-fitLog(val[!is.na(val)])[-1] #take out firt year for gimms
    
    scal<-sapply(mLog,function(k){k["scal"]})
    xmid<-sapply(mLog,function(k){k["xmid"]})
    Asym<-sapply(mLog,function(k){k["Asym"]})
    offset<-sapply(mLog,function(k){k["offset"]})
    
    NDVIgu_log<-logistic(xmid,alpha=Asym,beta=-xmid/scal,gamma=1/scal,offset=offset)[[2]] # do not take the scale, but the max slope
    NDVIgutime_log<-as.Date(xmid)
    NDVImax_log<-Asym+offset
    #NDVImaxtime_log<-as.Date(names(s1)[pos1]) # not really defined unless using double-logistic
    
    l<-logistic_optimum(alpha=Asym,beta=-xmid/scal,gamma=1/scal)
    NDVIgutimebeg_log<-sapply(l[[2]],"[",1)
    NDVIgutimeend_log<-sapply(l[[2]],"[",2)
    
    if(pl){ 
      plot(as.Date(names(val)),val,ylim=c(-0.2,1),xaxt="n")
      axis.Date(1,at=as.Date(paste0(substr(names(val),1,4),"-01-01")),las=2)
      lines(as.Date(names(val)),s0)
      points(as.Date(names(val)),s1,col="red",cex=0.5)
      lines(as.Date(names(val)),s1,col="red")
      abline(0,0)
      h<-sapply(mLog,function(k){k["offset"]})+sapply(mLog,function(k){k["Asym"]})/2
      axis.Date(1,at=NDVIgutime_log,las=2,cex.axis=0.7,col.axis=alpha("green4",0.5),format="%m-%d")
      points(NDVIgutime_sg,h,col="red",pch=16)
      points(NDVIgutime_log,h,col="green4",pch=16)
      points(NDVIgutime_sg,NDVIgu_sg,col="red",pch=16)
      points(NDVImaxtime_sg,NDVImax_sg,col="black",pch=16)
      points(NDVIgutimebeg_log,h,col="green4",pch=16)
      points(NDVIgutimeend_log,h,col="green4",pch=16)
    }
    ans<-rbind(NDVIgutime_sg,NDVIgu_sg,NDVImaxtime_sg=NDVImaxtime_sg[-1],NDVImax_sg=NDVImax_sg[-1],NDVIgu_log,NDVIgutime_log,NDVImax_log,NDVIgutimebeg_log,NDVIgutimeend_log)
    #names(ans)<-substr(doy,1,4)
    if(!i%%50){
      print(paste(i))
    }
    if(all(val==1)){
      ans[]<-NA
      ans
    }else{
      ans
    }
  })

### lr is a list that stores all rasters
lr<-vector(mode="list",length=nrow(cell[[1]]))
names(lr)<-dimnames(cell[[1]])[[1]]
for(i in 1:nrow(cell[[1]])){
  peak<-do.call("rbind",lapply(cell,function(j){
    j[i,]#-i[8,]
  }))
  rr<-subset(r,1:ncol(peak))
  if(any(grep("time",dimnames(cell[[1]])[[1]][i]))){
    rr<-setValues(rr,matrix(as.integer(format(as.Date(peak),"%j")),ncol=ncol(peak)))
  }else{#rr[]<-peak
    rr<-setValues(rr,peak)
  }
  names(rr)<-unique(substr(doy,1,4))[-1]
  lr[[i]]<-rr
}
  
### show all years for a given value
levelplot(lr$NDVIgutime_sg,col.regions=rasterTheme()$regions$col,cuts=99)
dif<-lr$NDVIgutime_log-mean(lr$NDVIgutime_log,na.rm=TRUE)
mm<-range(extract(dif),na.rm=TRUE)
levelplot(dif,at=seq(mm[1],mm[2],length.out=100),col.regions=colo.scale(1:99,c("blue3","white","tomato")))

#tmap_mode("view")
#tm_shape(pol)+tm_borders(lwd=1,alpha=0.5,col="black")+tm_shape(rr[[19]])+tm_raster(palette=rasterTheme()$regions$col,alpha=0.9,n=20)+tm_layout(basemaps=c("Esri.WorldImagery","Esri.WorldShadedRelief","Esri.NatGeoWorldMap"))

#writeRaster(lr[[6]],"C:/Users/rouf1703/Documents/rast.tif",format="GTiff",overwrite=TRUE)
#rast<-stack("C:/Users/rouf1703/Documents/r2014.tif")

#e<-extract(lr[[1]],pol)
#test<-unstack(lr[[2]])

### arrange data as a stack of different values for each year
lr2<-lapply(lr,unstack)
lr2<-lapply(seq_along(lr2[[1]]),function(i){
  res<-stack(lapply(lr2,function(j){
    j[[i]]
  }))  
  names(res)<-names(lr)
  res
})
names(lr2)<-1982:2015
#plot(lr2[[21]])

### write yearly rasters
for(i in seq_along(lr2)){
  writeRaster(lr2[[i]],paste0("C:/Users/rouf1703/Documents/","r",paste0(names(lr2)[i],".tif")),format="GTiff",overwrite=TRUE)
}
#rast<-stack("C:/Users/rouf1703/Documents/r1982.tif")

### suppose now yout want summarized values for each polygon
#plot(lr2[[21]][[1]])
#plot(pol,add=TRUE)
#e<-extract(lr2[[1]],pol[1,],fun=mean)


e<-extract(rast,pol[1,])

v<-velox(stack(lr[[6]]))
#v<-velox(rast)
test<-v$extract(pol,fun=function(i){mean(i,na.rm=TRUE)})
plot(pol,col=colo.scale(test[,6],terrain.colors(100)),border="white")


#pol2<-pol
#pol2@data<-as.data.frame(test)
#pol2<-st_as_sf(pol2)
#plot(pol2,max.plot=34,sf.colors=terrain.colors(100),border="white")



