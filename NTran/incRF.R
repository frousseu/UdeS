

library(data.table)
library(plyr)
library(signal)
library(randomForest)
library(scales)
library(FRutils)
library(weatherData)
library(lubridate)

### BOF data
Data <- read.table("C:/Users/User/Documents/incFR/DayNightData.txt", header = FALSE) # This is a .txt file that is not yet ready to be processed with 
Data2 <- txt.converter(table = Data, date.col = 1, time.col = 2, temp.col = 3, date.format = "%m/%d/%y", time.format = "%H:%M:%S" )


path<-"C:/Users/User/Documents/incFR"
f<-read.table("C:/Users/User/Documents/incFR/focales.txt",header=TRUE,stringsAsFactors=FALSE)
f$date<-paste(paste0("20",substr(f$date,7,8)),substr(f$date,4,5),substr(f$date,1,2),sep="-")
f$time<-ifelse(nchar(f$heure)==7,paste0("0",f$heure),f$heure)
f$DateTime<-paste(f$date,f$time)
f$dt<-as.POSIXct(f$DateTime)
f$on<-ifelse(f$Fsur.oeuf=="on",1,0)
#f$Fsur.oeuf<-ifelse(f$Fsur.oeuf=="debut","on",ifelse(f$Fsur.oeuf=="fin","off",f$Fsur.oeuf))
f<-dlply(f,.(nichoir,date))

l<-list.files(path,pattern=".csv")
x<-lapply(l,function(i){
 	x<-fread(paste(path,i,sep="/"))
 	x<-as.data.frame(x)
 	x$nichoir<-as.integer(substr(i,1,4))
 	x$dt<-as.POSIXct(x$DateTime)
 	x
	
})

x<-do.call("rbind",x)
x<-rbind(data.frame(Data2,nichoir=9999,dt=as.POSIXct(Data2$DateTime)),x)
x$date<-substr(x$DateTime,1,10)
x<-x[order(x$nichoir,x$dt),]
x<-dlply(x,.(nichoir,date))

x2<-lapply(names(x),function(i){
	xx<-x[[i]]
	m<-match(i,names(f))
	if(is.na(m)){
		on<-rep(NA,nrow(xx))
	}else{
		ff<-f[[m]]
		on<-sapply(xx$dt,function(j){
			di<-j-ff$dt
			if(all(di>0) || all(di<0)){
				NA
			}else{
				ff$on[tail(which(di>=0),1)]
			}
		})
	}
	xx$on<-on
	xx
})

names(x2)<-names(x)
x<-x2


plotInc<-function(x,mm=c(0,1)){
	d<-as.numeric(diff(x$dt,units="secs"))
	invisible(sapply(seq_along(d),function(i){
	  lines(c(x$dt[i],x$dt[i]+d[i]),mm[rep(x$on[i],2)+1],lwd=10,col=alpha("red",0.3),lend=1)  	
	}))
}



########################
### classify with RF
########################

x<-lapply(x,function(i){
  
  i$inc<-i$on
  i$temp<-i$Reading
  i$jul<-format(as.Date(substr(i$dt,1,10)),"%j")
  
  smoothn<-ceiling(nrow(i)/4)
  smoothn<-ifelse(smoothn%%2==0,smoothn+1,smoothn)
  
  s0<-sgolayfilt(i$temp,m=0,p=3,n=smoothn)
  s1<-sgolayfilt(i$temp,m=1,p=3,n=5)
  s2<-sgolayfilt(i$temp,m=2,p=3,n=5)
  s3<-sgolayfilt(i$temp,m=3,p=3,n=5)
  
  i$temps<-i$temp-s0
  i$tempss<-scale(i$temps)
  i$s0<-s0
  i$s1<-s1
  i$s2<-s2
  i$s3<-s3
  
  s11<-sgolayfilt(i$temps,m=1,p=3,n=5)
  s22<-sgolayfilt(i$temps,m=2,p=3,n=5)
  s33<-sgolayfilt(i$temps,m=3,p=3,n=5)
  
  s111<-sgolayfilt(i$temps,m=1,p=3,n=5)
  s222<-sgolayfilt(i$temps,m=2,p=3,n=5)
  s333<-sgolayfilt(i$temps,m=3,p=3,n=5)
  
  i$s11<-s11
  i$s22<-s22
  i$s33<-s33
  
  i$s111<-s111
  i$s222<-s222
  i$s333<-s333
  
  i$sign<-as.factor(ifelse(i$s1>=0,1,0))
  i$time<-sapply(strsplit(substr(i$dt,12,19),":"),function(j){j<-as.numeric(j);j[1]*3600+j[2]*60+j[3]})
  
  i$runif<-runif(nrow(i))
  
  #t<-getDetailedWeather("CYSC",i$date[1])
  #t$dt<-t$Time
  #ii<-cbind(i,id=1:nrow(i))
  #m<-merge(ii[,c("id","dt")],t,all=TRUE)
  #i$temp_ext<-na.spline(m$TemperatureC)[!is.na(m$id)]
  
  i
  
})

####
## FR

d<-do.call("rbind",x)
d$inc<-as.factor(d$inc)
d2<-d[!is.na(d$inc),]
#d2<-d2[sample(1:nrow(d2),300),]
d2<-tail(d2,900)

rf<-randomForest(inc~temp+temps+tempss+s1+s2+s3+s11+s22+s33+s111+s222+s333+sign+jul+time,data=d2,ntree=10000)
#rf<-randomForest(inc~runif,data=d[!is.na(d$inc),][1:100,],ntree=5)
varImpPlot(rf)
rf$confusion
reprtree:::plot.getTree(rf)

#p<-predict(rf,d,type="prob")[,2]
#lines(p-0.5,col=alpha("green4",0.5),lwd=4)
#p<-predict(rf,d)
#points((as.integer(p)-1.5)*2,pch=16,col="red",cex=0.4)



########

#layout(matrix(1:length(x),ncol=7,byrow=TRUE))

par(mar=c(2,2,0,0),oma=c(3,3,3,3),ask=TRUE)
for(i in seq_along(x)){
	m<-match(names(x)[i],names(f))
	
	bof
	bofdat<-x[[i]]
	bofdat$DateTime<-bofdat$dt
	bofdat<-bofdat[,1:2]
	bofdat$temp<-bofdat$Reading
	Boffed <- bof(object = bofdat, pT = 0.1, nT = -0.2, obs.lag = 2, plag = 1, nlag = 1, inctemp = 25, nwindow = 5, pwindow = 5, nullwindow = 0)
	#plot.bof.correct(Boffed, tdif = 2, date.format = "%Y-%m-%d", pch = 19, cex = 0.3)
	
	if(!is.na(m)){
		plot(x[[i]]$dt,x[[i]]$Reading,type="l",ylim=c(9,42),col=alpha("black",0.8),xlim=range(f[[m]]$dt)+c(-7200,7200),lwd=2,main=names(x)[i],xaxt="n")
		plotInc(f[[m]],mm=c(40.5,42.5))    
		plotInc(f[[m]],mm=c(40.1,42.1))  
	}else{
		plot(x[[i]]$dt,x[[i]]$Reading,type="l",ylim=c(9,42),col=alpha("black",0.8),lwd=2,main=names(x)[i],xaxt="n")	
	}
	axis.POSIXct(1,at=pretty(x[[i]]$dt,15))
	lines(x[[i]]$dt,(x[[i]]$on)+39,col="blue")
	p<-predict(rf,x[[i]],type="prob")[,2]
	lines(x[[i]]$dt,(p*7.5)+17.5,col=alpha("green4",0.3),lwd=2)
	points(x[[i]]$dt,(p*7.5)+17.5,col=alpha("green4",0.3),pch=16,cex=1)
	lines(x[[i]]$dt,x[[i]]$s0,col=alpha("black",0.5),lwd=2)
	lines(x[[i]]$dt,(x[[i]]$s1*3)+10,col=alpha("darkred",0.5),lwd=2)
	lines(x[[i]]$dt,(x[[i]]$s2*2)+10,col=alpha("red",0.5),lwd=2)
	lines(x[[i]]$dt,(x[[i]]$s3*2)+10,col=alpha("orange",0.5),lwd=2)

	abline(17.5,0,col=alpha("green4",0.3),lty=2,lwd=2)
	abline(21.25,0,col=alpha("green4",0.3),lty=2,lwd=2)
	abline(25,0,col=alpha("green4",0.3),lty=2,lwd=2)

	points(x[[i]]$dt,2*(ifelse(p>=0.5,1,0))+40.5,col=alpha("green4",0.3),pch=16,cex=1)
	#points(x[[i]]$dt,p+42.5,col=alpha(colo.scale(p-0.5,c("blue","white","red"),center=TRUE),0.9),pch=16,cex=0.7)
	points(Boffed$DateTime,2*(as.numeric(factor(Boffed$typ,levels=c("out","in")))-1)+40.1,col=alpha("blue",0.3),pch=16,cex=1)
	
	t<-getDetailedWeather("CYSC",x[[i]]$date[1])
	lines(t$Time,t$TemperatureC,col=alpha("red",0.1),lwd=10)
	
	#lines(x[[i]]$dt,x[[i]]$temp_ext,col=alpha("darkred",0.1),lwd=10)
	
	abline(10,0)
	
}




#########################
### The Bof
#########################

Data <- read.table("C:/Users/User/Documents/incFR/DayNightData.txt", header = FALSE) # This is a .txt file that is not yet ready to be processed with the function
head(Data)
str(Data)

Data2 <- txt.converter(table = Data, date.col = 1, time.col = 2, temp.col = 3, date.format = "%m/%d/%y", time.format = "%H:%M:%S" )
str(Data2)

Boffed <- bof(object = Data2, pT = 0.1, nT = -0.1, obs.lag = 2, plag = 1, nlag = 1, inctemp = 25, nwindow = 5, pwindow = 5, nullwindow = 0)
plot.bof.correct(Boffed, tdif = 2, date = "2012-05-25", date.format = "%Y-%m-%d", pch = 19, cex = 0.3)

BoutTable <- table.bof(Boffed, obs.lag = 2, dur.units = "mins", ID = "NestBox4-2010") 
#The argument ID is useful if you want to add an ID column (e.g. if you bind multiple nestbox)

hist(BoutTable$duration, breaks = 100)
LongBouts <- BoutTable[BoutTable$duration > 100, ]

DayTable <- day.sum(Boffed, obs.lag = 2, ID = "NestBox4-2010") 
DayTable2 <- day.sum(object = Boffed, obs.lag = 2, sunrise = "06:00:00", sunset = "20:00:00", ID = "NestBox4-2010", time.format = "%H:%M:%S")




#########################
#########################
#########################
#########################
							
y<-bofdat$temp[1:100]
TS<-1
n<-21
plot(y)
lines(y)

points(sgolayfilt(y,m=0,p=3,n=n,ts=TS),col="red");lines(sgolayfilt(y,m=0,p=3,n=n,ts=TS),col="red")
points(sgolayfilt(y,m=1,p=3,n=n,ts=TS)+23,col="blue");lines(sgolayfilt(y,m=1,p=3,n=n,ts=TS)+23,col="blue");abline(23,0,col="blue")








