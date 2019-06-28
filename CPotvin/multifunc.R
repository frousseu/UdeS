
library(readxl)
library(multifunc) # https://doi.org/10.1111/2041-210X.12143
library(FRutils)
library(scales)
library(viridisLite)
library(DHARMa)
library(colorspace)
library(ggplot2)
library(data.table)

#################################################################
### For the three blocks with common functions
#################################################################

### version for Multifunction and SME data.xlsx
#d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/CPotvin/Multifunction and SME data.xlsx",sheet=1,range="A3:AC25",col_names=TRUE))
#y<-unlist(as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/CPotvin/Multifunction and SME data.xlsx",sheet=1,range="A3:AC3",col_names=FALSE)))
#vars<-c("SOC3","CMIC","SRQ","SR","HB","LL","CWD","LD","AGB","CO")
#v<-rep(vars,times=c(3,2,2,1,3,3,3,3,3,3))


### version for Multifunction and SME data June13-2019.xlsx
d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/CPotvin/Multifunction and SME data June13-2019.xlsx",sheet=1,range="A3:AK24",col_names=FALSE))
y<-unlist(as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/CPotvin/Multifunction and SME data June13-2019.xlsx",sheet=1,range="A2:AK2",col_names=FALSE)))
vars<-c("LLB","SOC3M","SOC3","MB","CMIC","SRQ","SR","HB","CWD","LLD","AGB","CO","CRB")
varnames<-c("Leaf Litter Biomass","Soil Organic Carbon Mass","Soil Organic Carbon %","Mean Bas","CMIC","Soil Respiratory Quotient","Soil Respiration","Herbaceous Biomass","Coarse Woody Debris","Leaf Litter Decomposition","Aboveground Biomass","Canopy Openness","Coarse Root Biomass")
cat(paste(paste(varnames,paste0("(",vars,")")),collapse=", "))
v<-rep(vars,times=c(3,3,3,2,2,2,1,3,3,3,3,3,3))

y[grep("2001|2005|2006|2005-2006",y)]<-"Early"
y[grep("2011|2012|2013",y)]<-"Mid"
y[grep("2015|2016|2017",y)]<-"Late"

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

stock<-c("AGB","SRB","SOC3M","SOC3")
flux<-setdiff(vars,stock)

air<-c("AGB","CWD","CO","LLB","LLD") # are the two litter air?
soil<-setdiff(vars,air)

####################################
### reflecting

re<-c("SRQ","HB","LLD","CO")
d[,re]<--1*d[,re]


####################################
### subset
d2<-d # for full data in 2017
d<-d[,-match(c("SR","CMIC","SRQ","LLB","LLD","SOC3M","MB"),names(d))]

d2<-d2[d2$year=="Late",]
d2$year<-paste(d2$year,"full",sep=" - ")
d2<-d2[,-match(c("SOC3M","MB","CMIC","LLD"),names(d2))] # full data in 2017

########################################
### subset and remove all cols with NAs

ld<-split(d,d$year)
ld<-ld[c("Early","Mid","Late")]
ld<-c(ld,"Late - full"=list(d2))

ld<-lapply(ld,function(j){
  j[,!apply(j,2,function(i){any(is.na(i))})]
})

ldun<-ld
ld<-lapply(ld,function(j){
  xv<-names(j)[names(j)%in%vars]
  j[xv]<-lapply(j[xv],rescale,to=0:1)
  j
})


coly<-c("brown3","cornflowerblue","chartreuse4")
mgp<-c(3,0.25,0)
tcl<--0.2
ths<-seq(0.01,0.99,by=0.01)
ramp<-rev(viridis(100))

par(mfrow=c(2,3),mar=c(3,2,1,0),oma=c(2,4,0,2))
n1<-20
n2<-3
#mat<-matrix(1:6,nrow=2,byrow=FALSE)
#mat<-matrix(rep(c(1:6),each=n1),nrow=2,byrow=T)
mat<-cbind(
do.call("cbind",lapply(1:n1,function(i){1:2})),
do.call("cbind",lapply(1:n1,function(i){3:4})),
do.call("cbind",lapply(1:n1,function(i){5:6})),
do.call("cbind",lapply(1:n1,function(i){7:8})),
do.call("cbind",lapply(1:n2,function(i){rep(9,2)})))
layout(mat)
layout.show(9)
invisible(ans<-lapply(seq_along(ld),function(k){
  i<-ld[[k]]
  maxN<-min(table(i$diversity))+1 # Byrnes uses the min+1 of the minimum sample for a given diversity level (use to be maxN=3)
  print(maxN)
  th<-getFuncsMaxed(i,names(i)[names(i)%in%vars], threshmin=min(ths), threshmax=max(ths), prepend=c("plot","diversity","year"), maxN=maxN)
  #slopes<-getCoefTab(funcMaxed ~ diversity, data=th, coefVar=c("diversity"), family=quasipoisson(link="identity"))
    
  plot(0,0,ylim=c(0-0.25,sum(names(i)%in%vars)+0.5),xlim=range(th$diversity),type="n",xlab="",ylab="",font=2,font.lab=2,xaxt="n",yaxt="n")
  if(k%in%c(1:4)){
    axis(2,at=0:(sum(names(i)%in%vars)),label=0:(sum(names(i)%in%vars)),las=2,tcl=tcl,mgp=mgp)
  } 
  if(k==2){
    mtext("Diversity",1,line=2,cex=1.2,adj=1.15)
  }  
  axis(1,at=sort(unique(i$diversity)),label=sort(unique(i$diversity)),tcl=tcl,mgp=mgp)
  mtext(i$year[1],side=3,line=-1.5,font=2,adj=c(0.05))
  
  
  thr<-sort(unique(th$thresholds))
  l<-split(th,th$thresholds)
  cols<-alpha(colo.scale(1:length(l),ramp),1)
  
  
  ### fit models
  com<-lapply(seq_along(l),function(j){
    k<-l[[j]]
    k$year<-as.factor(k$year)
    m<<-glm(funcMaxed~diversity,data=k,family=quasipoisson(link="identity"),start=c(2,0.5))
    
    #v<-seq(1,5,by=0.1)
    #p<-predict(m,data.frame(diversity=v),type="response")
    #plot(jitter(k$diversity),jitter(k$funcMaxed))
    #lines(v,p)
    #Sys.sleep(0.00)
    #crap<<-k
    
    #sims<-simulateResiduals(m)
    #plot(sims)
    #print(m$deviance/m$df.residual)
    m
  })
  
  
  ### extract coefficients
  coe<-lapply(seq_along(com),function(h){
    m<-com[[h]]
    res<-as.data.frame(summary(m)$coef)
    suppressMessages(res<-cbind(th=l[[h]]$thresholds[1],year=i$year[1],var=row.names(res),res,as.data.frame(confint(m))))
    res
  })
  coe<-do.call("rbind",coe)
  coe<-coe[coe$var=="diversity",]
  
  
  ### plot and build CIs
  div<-seq(min(d$diversity),max(d$diversity),by=0.1)
  e<-data.frame(diversity=div)
  co<-lapply(seq_along(l),function(j){
    k<-l[[j]]
    #points(jitter(k$diversity,fac=0.5),jitter(k$funcMaxed,fac=1.25),col=cols[j],pch=1,cex=0.75,lwd=0.2) 
    k$year<-as.factor(k$year)
    m<-com[[j]]

    pred<-predict(m,newdata=e,type="response")
    e$pred<-pred
    res<-coe[j,]
    cis<-as.vector(res[,c("2.5 %","97.5 %")])
    sw<-if(all(cis>0) | all(cis<0)){sw<-TRUE}else{sw<-FALSE}
    #sw<-TRUE
    lines(div,e$pred,col=cols[j],lwd=5)
    #lines(div,e$pred,col=ifelse(sw,"red",alpha(cols[j],0.0)),lwd=1)
    points(0.92,head(e$pred,1),col=ifelse(sw,"black",alpha(cols[j],0.0)),pch=15)
    points(5.08,tail(e$pred,1),col=ifelse(sw,"black",alpha(cols[j],0.0)),pch=15)
    res
  })
  box(col="grey70")
  
  #browser()
  
  ### find and plot max slope
  mc<-c(which(coe[,"2.5 %"]>0)[1],which.max(coe[,"Estimate"]),rev(which(coe[,"2.5 %"]>0))[1])
  lapply(seq_along(mc),function(j){
    if(!is.na(mc[j])){
      m<-com[[mc[j]]]
      pred<-predict(m,newdata=e,type="response")
      e$pred<-pred
      if(j%in%c(1,3)){
        lines(div,e$pred,col="red",lwd=2,lty=2)
      }else{
        lines(div,e$pred,col="red",lwd=2)
      }
    }
  })
  

  ### plot of diversity coef vs. th with CIs
  #ylim<-range(unlist(coe[,c("Estimate","2.5 %","97.5 %")]))
  #ylim<-c(-0.2,0.4)
  ylim<-c(-0.5,1.5)
  plot(0,0,type="n",xlim=range(coe$th),ylim=ylim,ylab="",xlab="",yaxt="n",xaxt="n")
  if(k%in%1:4){
    axis(2,las=2,tcl=tcl,mgp=mgp)
  }
  if(k==2){
    mtext("Threshold",1,line=2,cex=1.2,adj=1.2)
  }
  axis(1,tcl=tcl,mgp=mgp)

  polygon(c(coe$th,rev(coe$th),coe$th[1]),c(coe$"2.5 %",rev(coe$"97.5 %"),coe$"2.5 %"[1]),col=alpha("black",0.15),border=NA)
  #lines(coe$th,coe$"2.5 %",col=alpha("black",0.3),lwd=1)
  #lines(coe$th,coe$"97.5 %",col=alpha("black",0.3),lwd=1)
  if(k%in%1:4){
    lapply(1:nrow(coe),function(h){
      lines(rep(coe$th[h],2),c(coe$"97.5 %"[h],coe$"2.5 %"[h]),col=cols[h],lwd=5,lend=2)  
      if(h%in%mc){
        if(match(h,mc)%in%c(1,3)){
          lines(rep(coe$th[h],2),c(coe$"97.5 %"[h],coe$"2.5 %"[h]),col="red",lwd=2,lend=2,lty=2)
        }else{
          lines(rep(coe$th[h],2),c(coe$"97.5 %"[h],coe$"2.5 %"[h]),col="red",lwd=2,lend=2)
        }
      }
    })
  }
  # white polygons to cut color gradient lines
  polygon(c(coe$th,rev(coe$th),coe$th[1]),c(coe$"2.5 %",rep(-10,length(coe[[1]])),-10),col="white",border=NA,xpd=FALSE)
  polygon(c(coe$th,rev(coe$th),coe$th[1]),c(coe$"97.5 %",rep(10,length(coe[[1]])),10),col="white",border=NA,xpd=FALSE)
  points(coe$th,coe$Estimate,col=alpha("black",0.75),pch=ifelse(coe$"Pr(>|z|)"<=0.05,1,1),cex=2,lwd=2)
  abline(0,0,lty=2)
  #m<-loess(Estimate~th,data=coe,span=0.7)
  #m<-gam(Estimate~s(th,k=6),data=coe)
  #p<-predict(m,data.frame(th=coe$th))
  #lines(coe$th,p,lwd=2,col=alpha(coly[i],0.5))
  box(col="grey70")

  
  mtext(bquote("   T"[min] == ~ .(format(coe$th[mc[1]],nsmall=2)) ~ "   T"[mde] == ~ .(format(coe$th[mc[2]],nsmall=2)) ~ "   T"[max] == ~ .(format(coe$th[mc[3]]),nsmall=2)),1,line=-2,adj=0.2)

  #legend("topleft",title="Coefficient and CIs\nof slope and smoother",legend=sapply(ans,function(i){i$year[1]}),col=alpha(coly,0.15),cex=1,pch=15,bty="n",pt.cex=3,y.intersp=1.5,inset=c(0.03,0.04))
  #legend("topleft",title="Coefficient and CIs\nof slope and smoother",legend=sapply(ans,function(i){i$year[1]}),col=alpha(coly,0.5),cex=1,pch=16,bty="n",pt.cex=1.25,y.intersp=1.5,inset=c(0.03,0.04))
}))


legend("topleft",legend=c("Slope estimate for a given threshold","Significant slope",expression("T"[min]*" & "*"T"[max]),expression("T"[mde])),lwd=c(NA,NA,2,2),pt.lwd=c(2,1,NA,NA),lty=c(NA,NA,2,1),pch=c(1,15,NA,NA),col=c("black","black","red","red"),cex=1.7,bty="n",pt.cex=c(2,1,NA,NA),inset=c(0.03,0))


## legend
plot(0,0,xlim=0:1,ylim=0:1,xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n",type="n")
rec<-c(-0.4,0,0.2,1)
legend_image <- as.raster(matrix(alpha(colo.scale(1:length(l),ramp),1), ncol=1))
rasterImage(legend_image,rec[1],rec[2],rec[3],rec[4],xpd=TRUE,lwd=1)
#at<-seq(min(ths),max(ths),length.out=10)
at<-c(min(ths),seq(0.1,0.9,by=.1),max(ths))
labs<-rescale(c(at,range(ths)),to=c(1,0))
labs<-at[1:length(at)]
text(x=rec[3],y=rev(seq(rec[2],rec[4],l=length(at))),cex=1,labels=round(labs,2),adj=c(-0.2,0.5),font=1,xpd=TRUE)

mtext("Nb of functions over threshold",2,line=2,cex=1.2,outer=TRUE,adj=0.85)  
mtext("Slope of diversity with CI",2,line=2,cex=1.2,outer=TRUE,adj=0.20)
mtext("Threshold",4,cex=1.2,outer=TRUE)



o1<-match(intersect(names(ld[[1]]),vars),vars)
o2<-order(varnames[o1])
cat(paste(paste(varnames[o1][o2],paste0("(",vars[o1][o2],")")),collapse=", "))


o1<-match(setdiff(names(ld[[4]]),names(ld[[1]])),vars)
o2<-order(varnames[o1])
cat(paste(paste(varnames[o1][o2],paste0("(",vars[o1][o2],")")),collapse=", "))


### histograms of p values
#ans<-coe
#ans<-split(ans,ans$var)
#par(mfrow=c(2,3))
#lapply(ans,function(i){
#  hist(i$"Pr(>|z|)",breaks=seq(0,1,by=0.05),main=i$var[1],border=NA,col="darkgreen")  
#})


#####################
### diversity vs. function

### with 0:1 variables used in multifunctionality
dats<-rbindlist(lapply(ld,function(i){
  melt(i,measure.vars = setdiff(names(i),c("id","size","diversity","year")),variable.name = "func", value.name = "y")  
}))
dats[,funcs:=paste0(varnames[match(func,vars)],"\n",func)]
dats[,periods:=factor(year,levels=unique(year))]

ggplot(aes(x=diversity,y=y),data=dats)+geom_point(size=3)+facet_grid(~periods~funcs)+theme_bw(base_size=15)+stat_smooth(method="lm",colour="black",size=2)+xlab("\nDiversity")+ylab("Standardized Value of Function\n") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0)) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))



### with standardization
l<-lapply(ldun,function(j){
  xv<-names(j)[names(j)%in%vars]
  j[xv]<-lapply(j[xv],rescale,to=0:1) # with the 0:1 scaling
  #j[xv]<-lapply(j[xv],function(i){(i-mean(i))/sd(i)}) # with standardization within periods
  j
})
dats<-rbindlist(lapply(l,function(i){
  melt(i,measure.vars = setdiff(names(i),c("id","size","diversity","year")),variable.name = "func", value.name = "y")  
}))
dats[,funcs:=paste0(varnames[match(func,vars)],"\n",func)]
dats[,periods:=factor(year,levels=unique(year))]

ggplot(aes(x=diversity,y=y),data=dats)+geom_point(size=3)+facet_grid(~periods~funcs)+theme_bw(base_size=15)+stat_smooth(method="lm",colour="black",size=2)+xlab("\nDiversity")+ylab("Standardized Value of Function within each period\n") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0))# +
  #scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))


### with 0:1 accross all periods
dats<-rbindlist(ldun,fill=TRUE)
vs<-intersect(names(dats),vars)
dats[,(vs):=lapply(.SD,rescale,to=0:1),.SDcols=vs]
dats<-melt(dats,measure.vars = setdiff(names(dats),c("id","size","diversity","year")),variable.name = "func", value.name = "y")  
dats[,funcs:=paste0(varnames[match(func,vars)],"\n",func)]
dats[,periods:=factor(year,levels=unique(year))]

ggplot(aes(x=diversity,y=y),data=dats)+geom_point(size=3)+facet_grid(~periods~funcs)+theme_bw(base_size=15)+stat_smooth(method="lm",colour="black",size=2,se=TRUE)+xlab("\nDiversity")+ylab("Standardized Value of Function across all periods\n") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0)) +
  scale_y_continuous(breaks=c(0,0.5,1))


### raw
dats<-rbindlist(ldun,fill=TRUE)
vs<-intersect(names(dats),vars)
#dats[,(vs):=lapply(.SD,rescale,to=0:1),.SDcols=vs]
dats<-melt(dats,measure.vars = setdiff(names(dats),c("id","size","diversity","year")),variable.name = "func", value.name = "y")  
dats[,funcs:=paste0(varnames[match(func,vars)],"\n",func)]
dats[,periods:=factor(year,levels=unique(year))]

ggplot(aes(x=diversity,y=y),data=dats)+geom_point(size=3)+facet_grid(~funcs~periods,scales="free")+theme_bw(base_size=15)+stat_smooth(method="lm",colour="black",size=2,se=TRUE)+xlab("\nDiversity")+ylab("Real Reflected Value of Each Function\n") +
  theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0))


#####################
### diversity vs. function

vs<-intersect(unlist(lapply(ld,names)),vars)
cols<- qualitative_hcl(length(vs), palette = "dark3")
ltys<-rep(c(1,2,3),length.out=length(vs))

par(mfrow=c(2,2))
lapply(seq_along(ld),function(k){
  i<-ld[[k]]
  plot(0,0,xlim=range(i$diversity),ylim=0:1,type="n")
  lapply(intersect(vars,names(i)),function(j){
    abline(lm(as.formula(paste0(j,"~","diversity")),data=i),col=cols[match(j,vs)],lty=ltys[match(j,vs)],lwd=5)
  })
  mtext(names(ld)[k],3,line=-1.5)
})
legend("bottom",lwd=5,col=cols,lty=ltys,legend=vs,ncol=4,bty="n",seg.len=7)


