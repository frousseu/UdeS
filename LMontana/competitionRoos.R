
library(scales)
library(mgcv)
library(scam)
library(lme4)
library(visreg)
library(ggeffects)
library(gridExtra)
library(ggplot2)
library(DHARMa)
library(data.table)
library(classInt)
library(geoR)
library(ncf)
library(AICcmodavg)
library(grid)

load("~/UdeS/Consultation/LMontana/Data.RData")
all<-as.data.table(m)

#plot(density(m$UDOI,bw="SJ"))
#lines(density(m.py$UDOI,bw="SJ"),col="red")

brks<-seq(0,max(m$UDOI),length.out=50)
m<-m[m$PY%in%c(0,1),] # not sure if should keep just 0s or 1s as well
h1<-hist(m$UDOI,breaks=brks,plot=FALSE)
h2<-hist(m.py$UDOI,breaks=brks,plot=FALSE)
#h<-rbind(h2$density,h1$density)
d<-data.frame(udoi=h1$mids,p=h2$density/h1$density,w=h1$counts/max(h1$counts),counts=h1$counts)
d<-d[!is.nan(d$p),]

##############################################
### histograms of densities
##############################################

hist(m$UDOI,breaks=brks,freq=FALSE,ylim=c(0,8),border=alpha("darkgreen",0.1),col=alpha("darkgreen",0.35),yaxt="n",xaxt="n",main="",xlab="",ylab="")
par(new=TRUE)
hist(m.py$UDOI,breaks=brks,freq=FALSE,ylim=c(0,8),border=alpha("darkred",0.1),col=alpha("darkred",0.35),yaxt="n",xaxt="n",main="",xlab="",ylab="")
axis(1,main="",pos=c(0,0),col="grey60")
axis(2,las=2,main="",pos=c(0,0),col="grey60")
legend("top",legend=c("With Young","Whole population"),fill=alpha(c("darkred","darkgreen"),0.35),border=alpha(c("darkred","darkgreen"),0.2),bty="n",cex=2)
mtext("UDOI",1)
mtext("Density",2)
title("Density of UDOI values for the population and the pairs that produced a young")
#points(d$udoi,d$p,cex=1.5,pch=16,bg=gray(0.5),col=gray(0.5))

################################################
### plot of ratios
################################################

plot(d$udoi,d$p,cex=rescale(d$w^1,c(1,20)),pch=21,bg=gray(0.5,0.5),col=gray(0.5,0.5),xlab="UDOI",ylab="Ratio of densities")
text(d$udoi,d$p,label=d$counts,cex=0.65,col=gray(0,0.75),adj=c(0.5,-1.3))
abline(0,0)

##################
### DECAY models
##################

x<-seq(0,max(d$udoi),by=0.01)


### gam
#g<-gam(p~s(udoi),data=d,weights=d$w)
#pred<-predict(g,data.frame(udoi=d$udoi))
#lines(d$udoi,pred,lty=2)

### shape-constrained gam
#d$p2<-sqrt(d$p) # this constrains over 0
#decay<-scam(p~s(udoi,bs="micv",k=10),data=d,weights=d$w,family=gaussian)
#pred<-predict(decay,data.frame(udoi=x),type="response")
#lines(x,pred,lwd=2)
#abline(0,0,lty=3)

### generalized logistic
#glo<-function(x,A,K,C,Q,B,v,off){
#    A+((K-A)/(C+Q*exp(-B*(x-off)))^(1/v))
#}
#lines(x,glo(x,A=0,K=2.5,C=1,Q=1,B=5,v=0.00001,off=-2),type="l")
#start<-c(A=0,K=2.5,C=1,Q=1,B=5,v=0.00001,off=-2)
#decay<-nls(p~glo(udoi,A,K,C,Q,B,v,off),data=d,start=start,control=list(minFactor=1e-12,maxiter=5000),weights=d$w)
#pred<-predict(decay,data.frame(udoi=x))
#lines(x,pred,col="red")



### exponential growth
#f<-function(x,A,k){
#  A*(1-exp(-k*x))
#}
#lines(x,f(x,A=0.7,k=4),type="l",col="red")
#start<-c(A=0.6,k=4)
#decay<-nls(p~f(udoi,A,k),data=d,start=start,control=list(minFactor=1e-12,maxiter=500),weights=d$w)
#pred<-predict(decay,data.frame(udoi=x))
#lines(x,pred,col="red")


### exponential growth gompertz curve
f<-function(x,a,b,c){
  a*(exp(-b*exp(-c*x)))
}
start<-c(a=2.5,b=4.1,c=4.6)
#lines(x,do.call("f",c(list(x),as.list(start))),col="magenta",lwd=2)
decay<-nls(p~f(udoi,a,b,c),data=d,start=start,control=list(minFactor=1e-12,maxiter=5000),weights=d$w)
pred<-predict(decay,data.frame(udoi=x))
lines(x,pred,col="red",lwd=2)

### gaussian
f<-function(x,s,n,r){
  (s-n)*(1-exp(-3*((x/r)^2)))+n
}
start<-c(s=2.5,n=0.1,r=0.75)
#lines(x,do.call("f",c(list(x),as.list(start))),col="magenta",lwd=2)
decay<-nls(p~f(udoi,s,n,r),data=d,start=start,control=list(minFactor=1e-12,maxiter=5000),weights=d$w)
pred<-predict(decay,data.frame(udoi=x))
lines(x,pred,col="blue",lwd=2)

### extract coefficients and plot min/max
co<-coef(decay)
mm<-do.call("f",c(list(c(0,max(d$udoi))),as.list(co)))
sapply(mm,abline,b=0,lty=3)

legend("topleft",lwd=c(2,2,1,NA),legend=c("Gompertz Growth Curve","Gaussian Variogram","Maximum value observed and value predicted at UDOI = 0","Ratio of density for a given UDOI bin with nb. of observations in bin"),col=c("red","blue","black",gray(0.5,0.5)),pch=c(NA,NA,NA,21),pt.bg=gray(0.5,0.5),lty=c(1,1,3,NA),cex=1.5,bty="n",inset=c(0.1,0.1))



#######################
#######################
#######################

fe<-c("Dam","cohort")

ratio<-mm[1]/mm[2] # this is the ratio between the value predicted at x=0 and the maximum value observed
ratio ### the ratio leaves some possibility of mating for almost no overlap

threshold<-0.16 # chosen threshold for the Th method

dt<-data.table(m)
dt[,w:=do.call("f",c(list(UDOI),as.list(co)))]
dt[,w:=rescale(w,to=c(ratio,1))] 
dt[,compSc:=.N,by=fe]
dt[,compWe:=lapply(.SD,sum),by=fe,.SDcols="w"]
dt[,compTh:=lapply(.SD,function(i){sum(i>=threshold)}),by=fe,.SDcols="UDOI"]

lt<-data.table(l)
lt<-lt[,lapply(.SD,mean),by=c("ID","cohort"),.SDcols=c("X","Y")]

dt<-merge(dt,lt,by.y=c("ID","cohort"),by.x=c("Dam","cohort"),all.x=TRUE)
comps<-unique(dt[,c(fe,"X","Y","compSc","compWe","compTh"),with=FALSE])
plot(comps$X,comps$Y,cex=rescale(comps$compTh,c(0.5,5)))

mt<-data.table(l)
mt<-mt[,lapply(.SD,mean),by=c("ID","cohort"),.SDcols=c("X","Y")]
mt<-merge(mt,unique(dt[,c("Male","cohort","leg","mass"),with=FALSE]),by.x=c("ID","cohort"),by.y=c("Male","cohort"),all.x=TRUE)
#mt<-unique(dt[,c("Male","cohort","X","Y","leg","mass"),with=FALSE])


ggplot(comps,aes(X,Y,size=compWe))+
  geom_point(alpha=0.5)+
  facet_wrap(~cohort)


wmean<-function(x,w){
  sum(x*w)/(sum(w))
}

wsd<-function(x,w){
  mu<-sum(x*w)/(sum(w))
  sqrt(sum(w*(x-mu)^2)/sum(w))
}

ma<-c("Male","cohort")
morphos<-c("mass","leg")
measures<-unique(dt[,c(ma,morphos),with=FALSE])
measures[,c("massSc","legSc"):=lapply(.SD,function(i){i-mean(i)}),.SDcols=morphos][,(morphos):=NULL]
dt<-merge(dt,measures,by=ma)

dt[,c("massTh","legTh"):=lapply(.SD,function(i){i-mean(i[UDOI>=threshold])}),by=fe,.SDcols=morphos]
dt[,c("massWe","legWe"):=lapply(.SD,function(i){i-wmean(i,w)}),by=fe,.SDcols=morphos]

#######################
#######################
#######################

fit0<-glmer(PY~UDOI+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
fit<-glmer(PY~UDOI+mass+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
fitSc<-glmer(PY~UDOI+massSc+compSc+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
fitTh<-glmer(PY~UDOI+massTh+compTh+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
fitWe<-glmer(PY~UDOI+massWe+compWe+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)

ml<-list(fit0=fit0,fit=fit,fitSc=fitSc,fitTh=fitTh,fitWe=fitWe)
aictab(ml)

v<-seq(min(dt$massWe),max(dt$massWe),length.out=50)
p<-predict(fitWe,data.frame(massWe=v,UDOI=mean(dt$UDOI)),re.form=NA,type="response")
plot(dt$massWe,jitter(dt$PY),ylim=c(-0.2,1.2),pch=16,col=gray(0,0.1))
lines(v,p)
abline(0,0,lty=3)

par(mfrow=c(1,2))
visreg(fit,"mass",scale="response",ylim=c(0,1))
visreg(fit,"leg",scale="response",ylim=c(0,1))
par(mfrow=c(1,1))

g1<-ggpredict(fitTh,terms="massTh [all]")
g2<-ggpredict(fitTh,terms="compTh [all]")
g3<-ggpredict(fitTh,terms="UDOI [all]")
lims<-c(0,0.2)
grid.arrange(grobs=lapply(list(g1,g2,g3),plot,limits=lims),ncol=3)

#sims<-simulateResiduals(fitwe)
#plot(sims)


#########################################
### Check for all cohorts
#########################################

val<-"legWe"
temp<-dt[dt$cohort%in%2000:2017,c("Male","Dam","cohort",val,"w"),with=FALSE]
temp[,y:=get(val)]
temp<-as.data.frame(temp)
temp[1:2]<-lapply(temp[1:2],as.character)

ggplot(temp,aes(Dam,Male,fill=y)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
  ggtitle("Relative size of males for a given female")+
  facet_wrap(~cohort,scales="free")+
  theme_light(base_size = 8) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5))

#all[which(all$Male=="288" & all$Dam=="258"),]


############################################
### Compare three measurements
############################################

temp<-dt[dt$cohort%in%2017,]
temp<-as.data.frame(temp)
temp[c("Male","Dam")]<-lapply(temp[c("Male","Dam")],as.character)

g<-list()
g[[1]]<-ggplot(temp,aes(Dam,Male,fill=legSc)) + ggtitle("2017 Global Scaling")
g[[2]]<-ggplot(temp,aes(Dam,Male,fill=legTh)) + ggtitle("2017 Threshold Scaling")
g[[3]]<-ggplot(temp,aes(Dam,Male,fill=legWe)) + ggtitle("2017 Weighted Scaling")

g<-lapply(g,function(i){
  i + geom_tile() + 
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,limits=range(unlist(temp[c("legSc","legTh","legWe")]))) +  
    theme_light(base_size = 10) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5))
})  

grid.arrange(grobs=g,ncol=3,top = textGrob("Variation of leg size (centered) in relation to male competitors for a given female",gp=gpar(fontsize=20,font=1)))



######################################################
### Check spatial structure using correlograms
######################################################

males<-all
males<-unique(males[,c("Male","cohort","mass","leg","age"),with=FALSE])

lt<-data.table(l)
lt<-lt[,lapply(.SD,mean,na.rm=TRUE),by=c("ID","cohort"),.SDcols=c("X","Y")]

males<-merge(males,lt,by.x=c("Male","cohort"),by.y=c("ID","cohort"),all.x=TRUE)

par(mfrow=c(4,6))
for(i in unique(males$cohort)){
  x<-males[cohort==i,]
  coords <- cbind(x$X,x$Y)
  v<-variog(coords=coords,data=x$leg,breaks=seq(0,1500,length.out=25),bin.cloud=TRUE)
  plot(v, main = "Variogram",type="b",bin.cloud=FALSE)
  sc <- spline.correlog(x=x$X, y=x$Y,z=x$leg, resamp=100, quiet=TRUE)
  plot(sc)
  mtext(i,3)
}


#plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
#plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
#file.copy(from=plots.png.paths, to="C:/Users/rouf1703/Downloads")
#plots.png.detials <- file.info(plots.png.paths)
#plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
#sorted.png.names <- gsub(plots.dir.path, "path_to_your_dir", row.names(plots.png.detials), fixed=TRUE)
#numbered.png.names <- paste0("C:/Users/rouf1703/Downloads/", 1:length(sorted.png.names), ".png")
# Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
#file.rename(from=sorted.png.names, to=numbered.png.names)





