
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

load("~/UdeS/Consultation/LMontana/Data.RData")
all<-as.data.table(m)

plot(density(m$UDOI,bw="SJ"))
lines(density(m.py$UDOI,bw="SJ"),col="red")

brks<-seq(0,max(m$UDOI),length.out=50)
#brks<-classIntervals(m$UDOI)
#brks<-c(0,0.001,0.002,0.003,0.004,0.005,0.01,0.03,0.05,seq(0.05,max(m$UDOI),length.out=50))
m<-m[m$PY%in%c(0,1),] # not sure if should do this, remove PY=1 to get ratio of densities
h1<-hist(m$UDOI,breaks=brks,plot=FALSE)
h2<-hist(m.py$UDOI,breaks=brks,plot=FALSE)

h<-rbind(h2$density,h1$density)

barplot(h,beside=TRUE,col=c("darkred","darkgreen"),border=NA,names.arg=paste(head(brks,-1),brks[-1],sep=" - "))
legend("top",legend=c("With Young","Without Young"),fill=c("darkred","darkgreen"),bty="n",cex=2)


#d<-data.frame(udoi=h1$mids,p=h2$density/(h1$density+h2$density),w=h1$counts/max(h1$counts),counts=h1$counts)
d<-data.frame(udoi=h1$mids,p=h2$density/h1$density,w=h1$counts/max(h1$counts),counts=h1$counts)
#d<-data.frame(udoi=h1$mids,p=h2$counts/h1$counts,w=h1$counts/max(h1$counts),success=h2$counts,failure=h1$counts-h2$counts)
d<-d[!is.nan(d$p),]

plot(d$udoi,d$p,cex=rescale(d$w^(1/1),c(0.5,10)),pch=21,bg=gray(0.5,0.5),col=gray(0.5,0.5))
#plot(d$udoi,d$p,pch=21,bg=gray(0.5,0.5),col=gray(0.5,0.5))
abline(0,0,lty=3)
#points(d$udoi,d$p,pch=1,col=gray(0.5,0.5))
#text(d$udoi,d$p,label=d$counts,cex=0.65,col=gray(0,0.75),adj=c(0.5,-0.5))

### gam
#g<-gam(p~s(udoi),data=d,weights=d$w)
#pred<-predict(g,data.frame(udoi=d$udoi))
#lines(d$udoi,pred,lty=2)

### shape-constrained gam
#d$p2<-sqrt(d$p) # this constrains over 0
#g<-scam(p~s(udoi,bs="micv",k=10),data=d,weights=d$w,family=gaussian)
#vals<-seq(min(brks),max(brks),by=0.01)
#pred<-predict(g,data.frame(udoi=vals),type="response")
#lines(vals,pred,lwd=2)
#abline(0,0,lty=3)

### generalized logistic
#glo<-function(x,A,K,C,Q,B,v,off){
#    A+((K-A)/(C+Q*exp(-B*(x-off)))^(1/v))
#}
#x<-seq(-10,10,by=0.01)
#lines(x,glo(x,A=0,K=2.5,C=1,Q=1,B=5,v=0.00001,off=-2),type="l")
#start<-c(A=0,K=2.5,C=1,Q=1,B=5,v=0.00001,off=-2)
#mn<-nls(p~glo(udoi,A,K,C,Q,B,v,off),data=d,start=start,control=list(minFactor=1e-12,maxiter=500),weights=d$w)
#pred<-predict(mn,data.frame(udoi=x))
#lines(x,pred,col="red")

### exponential growth
#ex<-function(x,A,k){
#  A*(1-exp(-k*x))
#}
#x<-seq(0,10,by=0.01)
#lines(x,ex(x,A=0.7,k=4),type="l",col="red")
#start<-c(A=0.6,k=4)
#mn<-nls(p~ex(udoi,A,k),data=d,start=start,control=list(minFactor=1e-12,maxiter=500),weights=d$w)
#pred<-predict(mn,data.frame(udoi=x))
#lines(x,pred,col="red")

### dgamma
#x<-seq(0,10,by=0.01)
#lines(x,cumsum(dgamma(x,3,8))/40,type="l",col="blue")
#start<-c(shape=3,rate=8,C=40)
#mn2<-nls(p~cumsum(dgamma(udoi,shape,rate))/C,data=d,start=start,control=list(minFactor=1e-12,maxiter=500),weights=d$w)
#pred<-predict(mn2,data.frame(udoi=x))
#lines(x,pred,col="blue")

### exponential growth gompertz curve
gom<-function(x,a,b,c){
  a*(exp(-b*exp(-c*x)))
}
x<-seq(0,1.5,by=0.01)
#lines(x,gom(x,a=2.5,b=4.1,c=8.6),col="magenta",lwd=5)
start<-c(a=2.5,b=4.1,c=4.6)
mn<-nls(p~gom(udoi,a,b,c),data=d,start=start,control=list(minFactor=1e-12,maxiter=5000),weights=d$w)
pred<-predict(mn,data.frame(udoi=x))
lines(x,pred,col="magenta",lwd=2)


#######################
#######################
#######################

a<-coef(mn)[["a"]]
b<-coef(mn)[["b"]]
c<-coef(mn)[["c"]]

mm<-gom(c(0,max(d$udoi)),a=a,b=b,c=c)
sapply(mm,abline,b=0,lty=2)

fe<-c("Dam","cohort")

ratio<-mm[1]/mm[2]
ratio ### the ratio leaves some possibility of mating for almost no overlap

threshold<-0.16

dt<-data.table(m)
dt[,w:=gom(UDOI,a=a,b=b,c=c)]
dt[,w:=rescale(w,to=c(ratio,1))] 
dt[,compSc:=.N,by=fe]
dt[,compWe:=lapply(.SD,sum),by=fe,.SDcols="w"]
dt[,compTh:=lapply(.SD,function(i){sum(i>=threshold)}),by=fe,.SDcols="w"]

lt<-data.table(l)
lt<-lt[,lapply(.SD,mean),by=c("ID","cohort"),.SDcols=c("X","Y")]

dt<-merge(dt,lt,by.y=c("ID","cohort"),by.x=c("Dam","cohort"),all.x=TRUE)
comps<-unique(dt[,c(fe,"X","Y","compn","compwe","compth"),with=FALSE])
plot(comps$X,comps$Y,cex=rescale(comps$compth,c(0.5,5)))

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

#n<-100
#x<-runif(n);x<-c(1,2,3)
#w<-runif(n,0,100);w<-c(1000,1,1)
#wmean(x,w)
#weighted.mean(x,w)
#wsd(x,w)
#weighted.sd(x,w)
#rank((x-mean(x))/(sd(x)))-rank((x-wmean(x,w))/(wsd(x,w)))


#dt<-data.table(Male=c(1,2,3,1,2,3,1,2,3),Dam=c(33,33,33,66,66,66,99,99,99),UDOI=c(0.2,0.2,0.2,0.2,0.5,1,1,0.5,0.1),cohort=2011,leg=rep(c(100,150,200),3),mass=rep(c(100,150,200),3))
#dt[,w:=rescale(UDOI,0:1)]
#dt[,w:=1]

ma<-c("Male","cohort")
morphos<-c("mass","leg")
measures<-unique(dt[,c(ma,morphos),with=FALSE])
measures[,c("massSc","legSc"):=lapply(.SD,function(i){i-mean(i)}),.SDcols=morphos][,(morphos):=NULL]
dt<-merge(dt,measures,by=ma)

#dt[,c("massth","legth"):=lapply(.SD,function(i){(i-mean(i[UDOI>=0.16]))/sd(i[UDOI>=0.16])}),by=fe,.SDcols=morphos]
dt[,c("massTh","legTh"):=lapply(.SD,function(i){i-mean(i[UDOI>=threshold])}),by=fe,.SDcols=morphos]
#dt[,c("masswe","legwe"):=lapply(.SD,function(i){(i-wmean(i,w))/wsd(i,w)}),by=fe,.SDcols=morphos]
dt[,c("massWe","legWe"):=lapply(.SD,function(i){i-wmean(i,w)}),by=fe,.SDcols=morphos]
#dt[,c("masszz","legzz"):=lapply(.SD,function(i){scale(i)[,1]}),.SDcols=morphos]
#dt[,c("masszz","legzz"):=lapply(.SD,function(i){i-mean(i)}),.SDcols=morphos]

#morphos<-sort(names(dt)[grep("mass|leg",names(dt))])
#dt[,paste0("q",morphos):=lapply(.SD,function(i){rank(i)/length(i)}),by=fe,.SDcols=morphos]
#dt[,paste0("re",morphos):=lapply(.SD,function(i){rescale(i,0:1)}),by=fe,.SDcols=morphos]
#dt[order(Male,cohort)]


#######################
#######################
#######################

fit0<-glmer(PY~UDOI+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
fit<-glmer(PY~UDOI+mass+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
#fitzz<-glmer(PY~UDOI+masszz+(1|Dam)+(1|Male)+(1|cohort),family=binomial,data=dt)
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



val<-"legWe"
#temp<-as.data.frame(dt[dt$cohort=="2012",c("Male","Dam",val),with=FALSE])
temp<-dt[dt$cohort%in%2000:2018,c("Male","Dam","cohort",val,"w"),with=FALSE]
temp[,y:=legWe*1]
#temp[,y:=lapply(.SD,function(i){i-mean(i)}),by="Male",.SDcols="y"]
#temp<-temp[order(cohort,legwe),]
#temp[,Male:=factor(Male,levels=unique(Male))]
temp<-as.data.frame(temp)
temp[1:2]<-lapply(temp[1:2],as.character)

ggplot(temp,aes(Dam,Male,fill=y)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
  ggtitle("Relative size of males for a given female")+
  facet_wrap(~cohort,scales="free")

all[which(all$Male=="288" & all$Dam=="258"),]


###############
###############
###############
males<-all
males<-unique(males[,c("Male","cohort","mass","leg","age"),with=FALSE])

lt<-data.table(l)
lt<-lt[,lapply(.SD,mean,na.rm=TRUE),by=c("ID","cohort"),.SDcols=c("X","Y")]

males<-merge(males,lt,by.x=c("Male","cohort"),by.y=c("ID","cohort"),all.x=TRUE)
#males<-unique(males[,c("Male","cohort","X","Y","mass","leg","age"),with=FALSE])

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





