
library(readxl)
library(plyr)
library(scales)



d<-read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/RBradley/Doc/PLFA data.xls",sheet=2,skip=6)

### format data
n<-names(d)
n<-n[!n%in%c("Peak name",NA,"NA")]
n<-gsub(" ","",gsub(" - ","_",n))
n<-sapply(strsplit(n,"_"),function(i){paste(i[2],i[1],sep="_")})
n<-rep(n,each=5)
names(d)<-c("peak",n)
d<-d[,-((ncol(d)-3):(ncol(d)))]

d<-as.data.frame(t(d),stringsAsFactors=FALSE)
names(d)<-unname(d[1,])
d<-d[-1,]
names(d)[1]<-"plot"
d$sample<-substr(d$plot,nchar(d$plot),nchar(d$plot))
d$plot<-gsub("i","I",substr(d$plot,1,2))
d$season<-sapply(strsplit(rownames(d),"_"),function(i){i[1]})
d$year<-as.integer(substr(sapply(strsplit(rownames(d),"_"),function(i){i[2]}),1,4))
d$pair<-ifelse(d$plot%in%c("SI","BI"),1,2)
keep<-c("plot","pair","sample","season","year")
d<-d[,c(keep,setdiff(names(d),keep))]
d[(length(keep)+1):ncol(d)]<-lapply(d[(length(keep)+1):ncol(d)],as.numeric)


# ne pas oublier d'enlever le control
ag<-names(rev(sort(colSums(d[,-seq_along(keep)],na.rm=TRUE))))[1:15] # prendre les 15 ag les plus abondants
d<-d[,c(keep,ag)]
#d[(length(keep)+1):ncol(d)]<-lapply(d[(length(keep)+1):ncol(d)],scale)



x<-ddply(d,.(pair,season,year),function(i){
  w1<-which(i$plot%in%c("BI","SR"))
  w2<-which(i$plot%in%c("SI","CO")) 
  n<-1:15
  x<-as.matrix(i[w1,-seq_along(keep)])[,n]
  y<-as.matrix(i[w2,-seq_along(keep)])[,n]
  center<-colMeans(y)
  SigmaInv<-solve(var(d[,-seq_along(keep)][,n]))
  sqrt(apply(x,1,function(i){(i-center) %*% SigmaInv %*% (i-center)}))
  #sqrt(mahalanobis(x,center=unname(center),cov=SigmaInv,tol=1e-80,inverted=TRUE))
})
x<-x[order(x$pair,x$year,ifelse(x$season=="spring",1,2)),]
x$time<-rep(1:5,2)
x<-x[,c("pair","season","year","time","6","7","8","9","10")]


#windows()
plot(1:5,1:5,ylim=c(0,max(as.numeric(unlist(x[,-(1:4)])),na.rm=TRUE)),type="n",xaxt="n",yaxt="n",xlab="Time",ylab="Mahalanobis Distance")
axis(1,at=1:5,label=paste(x$season[1:5],x$year[1:5],sep="_"))
axis(2,las=2)
for(i in 1:5){points(rep(x$time[i],5),unlist(x[i,5:9]),col=alpha("blue",0.5),pch=1,lwd=2,cex=1.5)}
for(i in 6:10){points(rep(x$time[i],5),unlist(x[i,5:9]),col=alpha("red",0.5),pch=1,lwd=2,cex=1.5)}
x2<-ddply(x,.(pair,time),function(i){
  cbind(i[rep(1,5),1:4],dm=unlist(i[1,5:9]))
})
m1<-lm(dm~time,data=x2[x2$pair==1,])
m2<-lm(dm~time,data=x2[x2$pair==2,])
abline(m1,col=alpha("blue",0.5),lwd=3)
abline(m2,col=alpha("red",0.5),lwd=3)
legend("topright",pch=1,lwd=3,pt.cex=1.5,col=alpha(c("red","blue",0.5)),legend=c("BI vs. SI","SR vs. CO"))
visreg(m1)
visreg(m2)
