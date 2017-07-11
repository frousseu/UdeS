
library(foreach)
library(doParallel)
library(signal)
library(robustbase)
library(MODIS)
library(plyr)


### This script is for extracting ndvi/evi metrics from RasterStack objects using a set of regions defined by polygons

# RasterStacks may be of the MODIS or GIMMS type
# A raster of julian dates is assumed to accompany each value raster

# y = a year as an integer
# j = the julian day to convert
# j2 = the julian day that marks the beginning of a block


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
  r<-rank(ifelse(max,-1,1)*x)
  res<-lapply(l,function(i){
    val<-sort(r[i])[1:n]
    index<-match(val,r)
    index   
  })
  res
}


#### test logistic



fitLog<-function(x,mmdate=c("12-01","09-15"),plot=FALSE){
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
    lo1<-list(Asym=0,xmid=min_xmid,scal=15,c=0.0)
    up1<-list(Asym=1,xmid=max_xmid,scal=40,c=0.8)
    m1<-tryCatch(nls(y~Asym/(1+exp((xmid-x)/scal))+c,data=d,start=list(Asym=0.5,xmid=min_xmid+(max_xmid-min_xmid)/2,scal=20,c=0.2),control=list(minFactor=1e-12,maxiter=500),lower=lo1,upper=up1,algorithm="port"),error=function(j){TRUE})
    if(!isTRUE(m1)){
      se<-seq(min(d$x),max(d$x),by=1) 
      if(plot){  
        #plot(as.Date(d$x),d$y)
        lines(as.Date(se),predict(m1,data.frame(x=se)),col=alpha("green4",0.5),lwd=4)
      }
      coef(m1)
    }else{
      NA  
    }
  })
  peak
}

fitGau<-function(x,mmdate=c("12-01","10-15"),plot=FALSE){
  years<-as.integer(unique(substr(names(x),1,4)))
  years<-years[years<=2016]
  l<-lapply(years,function(i){
    paste(c(i-1,i),mmdate,sep="-")  
  })
  peak<-lapply(l,function(i){
    sx<-x[which(names(x)>=i[1] & names(x)<=i[2])]
    d<-data.frame(y=sx,x=as.integer(as.Date(names(sx))))
    min_mu<-as.integer(as.Date(i[1])+120)
    max_mu<-as.integer(as.Date(i[2])-0)
    #lo1<-list(c=0.0)
    #up1<-list(c=0.7)
    m1<-tryCatch(nls(y~(a/(sigma*sqrt(2*pi)))*exp(-((x-mu)^2)/(2*sigma^2))+c,data=d,start=list(a=100,mu=max_mu-30,sigma=5,c=0.3),control=list(minFactor=1e-12,maxiter=500),algorithm="port"),error=function(j){TRUE})
    if(!isTRUE(m1)){
      se<-seq(min(d$x),max(d$x),by=1) 
      if(plot){  
        #plot(as.Date(d$x),d$y)
        lines(as.Date(se),predict(m1,data.frame(x=se)),col=alpha("blue4",0.5),lwd=4)
      }
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

### Ram mountain
ram<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc",layer="ram")
pol<-gBuffer(ram,width=0.0)
#pol<-bbox2pol(c(-72.26,-72.15,45.29,45.40))
pol<-SpatialPolygonsDataFrame(pol,data.frame(id=1),match.ID=FALSE)

#tmap_mode("view")
#tm_shape(pol)+tm_polygons(alpha=0.3)+tm_borders(lwd=5)+tm_layout(basemaps=c("Esri.WorldImagery"))

### Refuges Alberta BC
z<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/YPoisson/Doc",layer="largezone_BC_Alberta")
pol<-spTransform(z,CRS(proj4string(r)))


#############################################################
##### Get MODIS raster ######################################

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/MODIS/VI_16Days_250m_v6/Time_Series/RData/"

load(paste0(path,"MOD13Q1_MYD13Q1_NDVI_49_2000_1_2017_RData.RData"))
modis<-raster_ts #r5 is the initial raster used in plot ndvi3
load(paste0(path,"MOD13Q1_MYD13Q1_DOY_49_2000_1_2017_RData.RData"))
modis_jul<-raster_ts
load(paste0(path,"MOD13Q1_MYD13Q1_Rely_49_2000_1_2017_RData.RData"))
modis_rely<-raster_ts
load(paste0(path,"MOD13Q1_MYD13Q1_VI_QA_49_2000_1_2017_RData.RData"))
modis_QA<-raster_ts
rm(raster_ts)
l<-strsplit(names(modis),"_")
modis_doy<-as.Date(paste0(sapply(l,"[",3),"-01-01"))+as.integer(sapply(l,"[",4))-1

m<-match("MYD",substr(names(modis),1,3))
#keep<-m:length(names(modis))
keep<-1:length(names(modis))
#keep<-grep("MOD",names(r))

### subset modis to only get data from when Aqua was also used
modis<-subset(modis,keep)
modis_jul<-subset(modis_jul,keep)
modis_doy<-modis_doy[keep]

# Alpha and Terra data are offset by a certain amount in their blocks, but the specific date on which the ndvi is measured for a pixel may be earlier in a later block. Thus, although blocks are ordered, it does not mean that the specific dates are ordered. Ideally, order dates for each pixel one the observations are out of the complete matrix

#for(i in seq_along(dim(modis_jul))){
#  year<-unlist(strsplit(names(modis)[i],"_"))[[3]]
#  dateseq<-as.integer(seq.Date(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),by=1))
#  jul2n<-function(x){
#    dateseq[x]  
#  }
#  res<-calc(modis_jul[,,i],fun=jul2n)
#  modis_jul<-setValues(modis_jul,layer=i,  
#}


#day<-lapply(1980:2017,function(i){
#  seq.Date(as.Date(paste0(i,"-01-01")),as.Date(paste0(i,"-12-31")),by=1)
#})



#rc<-SpatialPointsDataFrame(coordinates(modis),proj4string=CRS(proj4string(modis)),data.frame(id=seq_len(dim(modis)[1]*dim(modis)[2])))


#############################################################
##### Get GIMMS raster ######################################

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/GIMMS/"

x<-list.files(path)
ts<-monthlyIndices(x, version = 1, timestamp = T)
gimms_doy<-as.Date(as.character(ts+round(c(diff(ts),17)/2,0)))
gimms<-rasterizeGimms(x=paste0(path,x),ext=pol,cores=6L) # clipping
names(gimms)<-as.character(gimms_doy)
gimms_jul<-stack(setValues(gimms,as.integer(format(rep(gimms_doy,each=ncell(gimms)),"%j"))))







#############################################################
##### Extract raster values for each region #################

r<-gimms # determine series to use
rd<-gimms_jul
doy<-gimms_doy


registerDoParallel(6) 
getDoParWorkers()

# build a list of matrices of each pixel ndvi value for each polygon
v<-foreach(i=1:length(pol),.packages=c("raster","sp")) %dopar% {
  extract(r,pol[i,])[[1]]  
}
### remove all NA values for some pixels
#v<-lapply(v,function(i){
#  k<-!apply(i,1,function(j){
#    all(is.na(j))  
#  })
#  i[k,]
#})

vd<-foreach(i=1:length(pol),.packages=c("raster","sp")) %dopar% {
  extract(rd,pol[i,])[[1]]  
}
#vd<-lapply(vd,function(i){
#  k<-!apply(i,1,function(j){
#    all(is.na(j))  
#  })
#  i[k,]
#})


#############################################################
##### Compute metrics #######################################


peak_cell<-lapply(seq_along(v),function(i){
  lapply(1:nrow(v[[i]]),function(j){
    pl<-FALSE
    val<-v[[i]][j,]/1
    if(all(is.na(val))){
      val[seq_along(val)]<-1  
    }
    jul<-vd[[i]][j,]
    names(jul)<-doy
    #jul2<-as.integer(sapply(strsplit(names(jul),"_"),"[",4))
    jul2<-as.integer(format(doy,"%j"))
    k<-is.na(jul)
    if(any(k)){
      jul[k]<-jul2[k] # missing values in precise julian days are given the beginning of the block  
    }
    name<-unname(as.Date(jul2nd(as.integer(substr(doy,1,4)),jul,jul2),origin="1970-01-01")) # blocks are ordered, but not necessarily precise dates
    names(val)<-name
    o<-order(name)
    val<-val[o]
    #val2<-val
    #sup<-c(diff(val,lag=2)>0.4,FALSE)
    #val[sup]<-NA
    s1<-sgolayfilt(na.spline(val),n=41,p=3,m=1)
    names(s1)<-names(val)
    #pos<-findminmax(s1,n=5,beg="03-01",end="07-01")
    #sg<-as.Date(sapply(pos,function(k){mean(as.Date(names(s1)[k]))}))
    pos<-unlist(findminmax(s1,n=1,beg="03-01",end="07-01"))
    sg<-as.Date(names(s1)[pos])
    s0<-sgolayfilt(na.spline(val),n=21,p=3,m=0)
    #data.frame(na.spline(v[[i]][j,]),doy,names(s1),max=as.numeric(seq_along(doy)%in%pos)) # verif
    if(!j%%20)
      print(paste(i,j))
    if(pl){  
      plot(as.Date(names(val)),val,ylim=c(-0.2,1),xaxt="n")
      #points(as.Date(names(val[sup])),na.spline(val)[sup],pch=16)
      #points(as.Date(names(val2[sup])),val2[sup],pch=8)
      axis.Date(1,at=as.Date(paste0(substr(names(val),1,4),"-01-01")),las=2)
      lines(as.Date(names(val)),s0)
      points(as.Date(names(val)),s1*7,col="red",cex=0.5)
      abline(0,0)
    }
    mLog<-fitLog(val[!is.na(val)],plot=pl)[-1] #take out firt year for gimms
    logi<-as.Date(sapply(mLog,function(k){k["xmid"]}))
    #xx<<-val[!is.na(val)]
    #gaus<-as.Date(fitGau(val[!is.na(val)],plot=pl))[-1]
    if(pl){  
      h<-sapply(mLog,function(k){k["c"]})+sapply(mLog,function(k){k["Asym"]})/2
      axis.Date(1,at=logi,las=2,cex.axis=0.7,col.axis=alpha("green4",0.5),format="%m-%d")
      points(sg,h,col="red",pch=16)
      points(logi,h,col="green4",pch=16)
    }
    ans<-rbind(sg,logi)#,gaus)
    if(all(val==1)){
      ans[]<-NA
      ans
    }else{
      ans
    }
  })
})

peak1<-lapply(peak_cell,function(i){
  ii<-do.call("rbind",lapply(i,function(j){j[1,]}))
  as.Date(colMeans(ii,na.rm=TRUE),origin="1970-01-01")  
})
peak2<-lapply(peak_cell,function(i){
  ii<-do.call("rbind",lapply(i,function(j){j[2,]}))
  as.Date(colMeans(ii,na.rm=TRUE),origin="1970-01-01")  
})

peak1g<-peak1
peak2g<-peak2

#peak1m<-peak1
#peak2m<-peak2

peak22<-sapply(peak2,function(i){median(as.integer(format(i,"%j")),na.rm=TRUE)})
plot(pol,col=colo.scale(peak22))
as.Date(range(peak22),"1970-01-01")
#tmap_mode("view")
#tm_shape(pol)+tm_borders(lwd=2,col="red")+tm_layout(basemaps=c("Esri.WorldImagery"))


### gimms and modis compare (peak1 peak2)

j1g<-as.Date(as.integer(format(peak1g[[1]],"%j")))
j2g<-as.Date(as.integer(format(peak2g[[1]],"%j")))
j1m<-as.Date(as.integer(format(peak1m[[1]],"%j")))
j2m<-as.Date(as.integer(format(peak2m[[1]],"%j")))
ylim<-range(as.Date(c("1970-04-01","1970-07-01")))


#j1g<-scale(j1g)[,1]
#j2g<-scale(j2g)[,1]
#j1m<-scale(j1m)[,1]
#j2m<-scale(j2m)[,1]
#ylim<-range(c(j1g,j2g))

plot(as.integer(substr(peak1g[[1]],1,4)),j1g,pch=16,col=alpha("red",0.3),ylim=ylim,type="l",las=2,lwd=4,xlim=c(1980,2016))
lines(as.integer(substr(peak2g[[1]],1,4)),j2g,pch=16,col=alpha("red",0.3),lwd=4,lty=2)
#lines(as.integer(substr(peak1m[[1]],1,4)),j1m,pch=16,col=alpha("blue",0.3),ylim=range(c(j1m,j2m)),type="l",las=2,lwd=4)
#lines(as.integer(substr(peak2m[[1]],1,4)),j2m,pch=16,col=alpha("blue",0.3),lwd=4,lty=2)
abline(0,0)


gimmsSG<-peak1g[[1]][match(1981:2016,substr(peak1g[[1]],1,4))]
gimmsLO<-peak2g[[1]][match(1981:2016,substr(peak2g[[1]],1,4))]
modisSG<-peak1m[[1]][match(1981:2016,substr(peak1m[[1]],1,4))]
modisLO<-peak2m[[1]][match(1981:2016,substr(peak2m[[1]],1,4))]

res<-data.frame(years=1981:2016,gimmsSG,gimmsLO,modisSG,modisLO)
res$modisSG[res$years==2000]<-NA #

res2<-ddply(res,.(years),function(i){format(i[-1],"%j")})
names(res2)[2:ncol(res2)]<-paste0(names(res2)[2:ncol(res2)],"jul")

res<-merge(res,res2)
fwrite(res,"C:/Users/rouf1703/Documents/UdeS/Consultation/L-ARenaud/Doc/greenup_ts_all.csv",row.names=FALSE,sep=";")










