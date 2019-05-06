
library(readxl)
library(multifunc) # https://doi.org/10.1111/2041-210X.12143
library(FRutils)
library(scales)

d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/CPotvin/Multifunction and SME data.xlsx",sheet=1,range="A3:AC25",col_names=TRUE))
y<-unlist(as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/CPotvin/Multifunction and SME data.xlsx",sheet=1,range="A3:AC3",col_names=FALSE)))

vars<-c("SOC3","CMIC","RC","SR","HAW","LL","CWD","LD","AGB","CO")
v<-rep(vars,times=c(3,2,2,1,3,3,3,3,3,3))
y[grep("2001|2005|2006|2005-2006",y)]<-"2001-2006"
y[grep("2011|2012|2013",y)]<-"2011-2013"
y[grep("2015|2016|2017",y)]<-"2015-2017"

names(d)<-c("id","size","diversity",v)


l<-lapply(vars,function(i){
  w<-which(names(d)%in%i)
  x<-do.call("rbind",lapply(1:length(w),function(k){d[,c("id","size","diversity")]}))
  data.frame(x,year=rep(y[w],each=nrow(d)),metric=i,value=unname(unlist(d[,w])),stringsAsFactors=FALSE)
})
d<-do.call("rbind",l)

d<-reshape(d, idvar = setdiff(names(d),c("metric","value")), timevar = "metric", direction = "wide")
names(d)<-gsub("value.","",names(d))


##################
### type of var

stock<-c("AGB","RB","SOC3")
flux<-setdiff(vars,stock)

air<-c("AGB","CWD","CO","LL","LD") # are the two litter air?
soil<-setdiff(vars,air)

####################################
### reflecting

re<-c("RC","HAW","LD","CO")
d[,re]<--1*d[,re]


####################################
### subset
d<-d[,-match(c("SR","CMIC","RC","LL","LD"),names(d))]


########################################
### subset and remove all cols with NAs

ld<-split(d,d$year)
#ld<-list(d)

ld<-lapply(ld,function(j){
  j[,!apply(j,2,function(i){any(is.na(i))})]
})

#ld<-lapply(ld,function(j){xv<-names(j)[names(j)%in%vars];for(i in xv){j[,i]<-1};j})

ld<-lapply(ld,function(j){
  xv<-names(j)[names(j)%in%vars]
  j[xv]<-lapply(j[xv],rescale,to=0:1)
  j
})

coly<-c("brown3","cornflowerblue","chartreuse4")
mgp<-c(3,0.25,0)
tcl<--0.2
ths<-seq(0.01,0.99,by=0.01)
  
par(mfrow=c(2,2),mar=c(4,4,1,3))
invisible(ans<-lapply(seq_along(ld),function(k){
  i<-ld[[k]]
  th<-getFuncsMaxed(i,names(i)[names(i)%in%vars], threshmin=min(ths), threshmax=max(ths), prepend=c("plot","diversity","year"), maxN=3)
  #slopes<-getCoefTab(funcMaxed ~ diversity, data=th, coefVar=c("diversity"), family=quasipoisson(link="identity"))
    
  plot(0,0,ylim=c(0-0.25,sum(names(i)%in%vars)+0.5),xlim=range(th$diversity),type="n",xlab="",ylab="",font=2,font.lab=2,xaxt="n",yaxt="n")
  title(ylab="Nb of functions over threshold",line=1,cex.lab=1.2)
  title(xlab="Diversity",line=1.2,cex.lab=1.2)
  axis(1,at=sort(unique(i$diversity)),label=sort(unique(i$diversity)),tcl=tcl,mgp=mgp)
  axis(2,at=0:(sum(names(i)%in%vars)),label=0:(sum(names(i)%in%vars)),las=2,tcl=tcl,mgp=mgp)
  mtext(i$year[1],side=3,line=-1.5,font=2,adj=c(0.95))
    
  thr<-sort(unique(th$thresholds))
  l<-split(th,th$thresholds)
  ramp<-c("grey85",coly[k],"grey1")
  cols<-alpha(colo.scale(1:length(l),ramp),0.75)
  co<-lapply(seq_along(l),function(j){
    k<-l[[j]]
    points(jitter(k$diversity,fac=0.5),jitter(k$funcMaxed,fac=1.25),col=cols[j],pch=1,cex=0.75,lwd=0.2) 
    k$year<-as.factor(k$year)
    print(j)
    m<-glm(funcMaxed~diversity,data=k,family=poisson(link="log"))
  
    div<-seq(1,5,by=0.1)
    yb<-levels(k$year)
    e<-expand.grid(diversity=div,year=yb)
    pred<-predict(m,newdata=e,type="response")
    e$pred<-pred
    lp<-split(e,e$year)
    lapply(lp,function(p){
      lines(div,p$pred,col=cols[j],lwd=2)
    })
    res<-as.data.frame(summary(m)$coef)
    res<-cbind(th=k$thresholds[1],year=i$year[1],var=row.names(res),res,as.data.frame(confint(m)))
    res
  })
  ## legend
  legend_image <- as.raster(matrix(alpha(colo.scale(1:length(l),ramp),0.75), ncol=1))
  rec<-c(par("usr")[2]+0.04,par("usr")[3],par("usr")[2]+0.13,par("usr")[4])
  rasterImage(legend_image,rec[1],rec[2],rec[3],rec[4],xpd=TRUE,lwd=1)
  #at<-seq(min(ths),max(ths),length.out=10)
  at<-c(min(ths),seq(0.1,0.9,by=.1),max(ths))
  labs<-rescale(c(at,range(ths)),to=c(1,0))
  labs<-at[1:length(at)]
  text(x=rec[3],y=rev(seq(rec[2],rec[4],l=length(at))),cex=0.8,labels=round(labs,2),adj=c(-0.2,0.5),font=0,xpd=TRUE)
  box(col="grey70")
  do.call("rbind",co)
}))
  
  

### plot of diversity coef vs. th with CIs
ans2<-lapply(ans,function(i){i[i$var=="diversity",]})
ylim<-range(unlist(do.call("rbind",ans2)[,c("Estimate","2.5 %","97.5 %")]))
ylim<-c(-0.2,0.4)
plot(0,0,type="n",xlim=range(do.call("rbind",ans2)$th),ylim=ylim,ylab="",xlab="",yaxt="n",xaxt="n")
title(ylab="Slope of diversity",line=1.7,cex.lab=1.2)
title(xlab="Threshold",line=1.2,cex.lab=1.2)
axis(1,tcl=tcl,mgp=mgp)
axis(2,las=2,tcl=tcl,mgp=mgp)
abline(0,0,lty=2)
lapply(seq_along(ans2),function(i){
  x<-ans2[[i]]
  polygon(c(x$th,rev(x$th),x$th[1]),c(x$"2.5 %",rev(x$"97.5 %"),x$"2.5 %"[1]),col=alpha(coly[i],0.15),border=NA)
  lines(x$th,x$"2.5 %",col=alpha(coly[i],0.7),lwd=1)
  lines(x$th,x$"97.5 %",col=alpha(coly[i],0.7),lwd=1)
  points(x$th,x$Estimate,col=alpha(coly[i],0.5),pch=ifelse(x$"Pr(>|z|)"<=0.05,16,16),cex=1.25,lwd=1)
  m<-loess(Estimate~th,data=x,span=0.7)
  #m<-gam(Estimate~s(th,k=6),data=x)
  p<-predict(m,data.frame(th=x$th))
  lines(x$th,p,lwd=2,col=alpha(coly[i],0.5))
  box(col="grey70")
})
legend("topleft",title="Coefficient and CIs\nof slope and smoother",legend=sapply(ans,function(i){i$year[1]}),col=alpha(coly,0.15),cex=1,pch=15,bty="n",pt.cex=3,y.intersp=1.5,inset=c(0.03,0.04))
legend("topleft",title="Coefficient and CIs\nof slope and smoother",legend=sapply(ans,function(i){i$year[1]}),col=alpha(coly,0.5),cex=1,pch=16,bty="n",pt.cex=1.25,y.intersp=1.5,inset=c(0.03,0.04))



### histograms of p values
ans2<-do.call("rbind",ans)
ans2<-split(ans2,ans2$var)
par(mfrow=c(2,3))
lapply(ans2,function(i){
  hist(i$"Pr(>|z|)",breaks=seq(0,1,by=0.05),main=i$var[1],border=NA,col="darkgreen")  
})


###########################################
###########################################
###########################################

############################################################
### for having the diversity:year interaction in the model

#ld<-split(d,d$year)
ld<-list(d)

ld<-lapply(ld,function(j){
  j[,!apply(j,2,function(i){any(is.na(i))})]
})

#ld<-lapply(ld,function(j){xv<-names(j)[names(j)%in%vars];for(i in xv){j[,i]<-1};j})

ld<-lapply(ld,function(j){
  xv<-names(j)[names(j)%in%vars]
  j[xv]<-lapply(j[xv],rescale,to=0:1)
  j
})

ld[[1]]$year<-factor(ld[[1]]$year)

coly<-c("brown3","cornflowerblue","chartreuse4")
mgp<-c(3,0.25,0)
tcl<--0.2
ths<-seq(0.05,0.95,by=0.01)

par(mfrow=c(2,2),mar=c(4,4,1,3))

  i<-ld[[1]]
  th<-getFuncsMaxed(i,names(i)[names(i)%in%vars], threshmin=min(ths), threshmax=max(ths), prepend=c("plot","diversity","year"), maxN=3)
  #slopes<-getCoefTab(funcMaxed ~ diversity, data=th, coefVar=c("diversity"), family=quasipoisson(link="identity"))
  
  plot(0,0,ylim=c(0-0.25,sum(names(i)%in%vars)+0.5),xlim=range(th$diversity),type="n",xlab="",ylab="",font=2,font.lab=2,xaxt="n",yaxt="n")
  title(ylab="Nb of functions over threshold",line=1,cex.lab=1.2)
  title(xlab="Diversity",line=1.2,cex.lab=1.2)
  axis(1,at=sort(unique(i$diversity)),label=sort(unique(i$diversity)),tcl=tcl,mgp=mgp)
  axis(2,at=0:(sum(names(i)%in%vars)),label=0:(sum(names(i)%in%vars)),las=2,tcl=tcl,mgp=mgp)
  #mtext(i$year[1],side=3,line=-1.5,font=2,adj=c(0.95))
  
  thr<-sort(unique(th$thresholds))
  l<-split(th,th$thresholds)
  
  #co<-lapply(seq_along(l),function(j){
  for(j in seq_along(l)){  
    k<-l[[j]]
    #points(jitter(k$diversity,fac=0.5),jitter(k$funcMaxed,fac=1.25),col=cols[j],pch=1,cex=0.75,lwd=0.2) 
    k$year<-as.factor(k$year)
    print(j)
    m<-glm(funcMaxed~diversity*year,data=k,family=poisson)
    print(summary(m))
    div<-seq(1,5,by=0.1)
    yb<-levels(k$year)
    e<-expand.grid(diversity=div,year=yb)
    pred<-predict(m,newdata=e,type="response")
    e$pred<-pred
    lp<-split(e,e$year)
    lapply(seq_along(lp),function(y){
      ramp<-c("grey80",coly[y])
      cols<-alpha(colo.scale(1:length(l),ramp),0.75)
      p<-lp[[y]]
      lines(div,p$pred,col=cols[j],lwd=2)#col=cols[j],lwd=2)
    })
    res<-as.data.frame(summary(m)$coef)
    #res<-cbind(th=k$thresholds[1],year=i$year[1],var=row.names(res),res,as.data.frame(confint(m)))
    res
    visreg2d(m,"diversity","year",scale="response")
  }
  ## legend
  legend_image <- as.raster(matrix(alpha(colo.scale(1:length(l),ramp),0.75), ncol=1))
  rec<-c(par("usr")[2]+0.04,par("usr")[3],par("usr")[2]+0.13,par("usr")[4])
  #rasterImage(legend_image,rec[1],rec[2],rec[3],rec[4],xpd=TRUE,lwd=1)
  #at<-seq(min(ths),max(ths),length.out=10)
  at<-c(min(ths),seq(0.1,0.9,by=.1),max(ths))
  labs<-rescale(c(at,range(ths)),to=c(1,0))
  labs<-at[1:length(at)]
  text(x=rec[3],y=rev(seq(rec[2],rec[4],l=length(at))),cex=0.8,labels=round(labs,2),adj=c(-0.2,0.5),font=0,xpd=TRUE)
  box(col="grey70")
  do.call("rbind",co)









