
library(readxl)
library(plyr)
library(scales)
library(vegan)
library(devtools)
library(ggbiplot)



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
d$time<-match(paste0(d$season,d$year),c("spring2009","fall2009","spring2010","fall2010","spring2011"))
keep<-c("plot","pair","sample","season","year","time")
d<-d[,c(keep,setdiff(names(d),keep))]
d[(length(keep)+1):ncol(d)]<-lapply(d[(length(keep)+1):ncol(d)],as.numeric)

n<-15
# ne pas oublier d'enlever le control
ag<-names(rev(sort(colSums(d[,-seq_along(keep)],na.rm=TRUE))))[1:n] # prendre les 15 ag les plus abondants
d<-d[,c(keep,ag)]
#d[(length(keep)+1):ncol(d)]<-lapply(d[(length(keep)+1):ncol(d)],scale)


x<-ddply(d,.(pair,season,year,time),function(i){
  w1<-which(i$plot%in%c("BI","SR"))
  w2<-which(i$plot%in%c("SI","CO")) 
  nn<-n #select all 15 ag (nn<-n) or not
  x<-as.matrix(i[w1,-seq_along(keep)])[,1:n][,1:nn]
  y<-as.matrix(i[w2,-seq_along(keep)])[,1:n][,1:nn]
  center<-colMeans(y)
  SigmaInv<-solve(var(d[,-seq_along(keep)][,1:n][,1:nn]))
  sqrt(apply(x,1,function(i){(i-center) %*% SigmaInv %*% (i-center)}))
  #sqrt(mahalanobis(x,center=unname(center),cov=SigmaInv,tol=1e-80,inverted=TRUE))
})
x<-x[order(x$pair,x$year,ifelse(x$season=="spring",1,2)),]
#x$time<-rep(1:5,2)
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
#visreg(m1)
#visreg(m2)

### pca
temp<-d[d$plot%in%c("SR","CO"),]
pca<-prcomp(temp[,7:ncol(d)],scale=TRUE)
#mod<-rda(temp[,7:ncol(d)],scale=TRUE)
par(mar=c(2,2,2,2))
biplot(pca,col=c("white","grey80"),expand=1.5,xlim=range(pca$rot[,1]),ylim=range(pca$rot[,2]))
box("plot",col="black")
#biplot(pca$rotation,pca$x,col=c("grey80","grey80"),expand=0.5)
text(pca$x[,1],pca$x[,2],temp$plot,col=alpha(temp$time+1,0.5),font=2,cex=1.25)
#text(mod$CA$v[,1],mod$CA$v[,2],temp$plot,col=temp$time+1,font=2)
#points(pca$x[,1],pca$x[,2],col=temp$time,cex=4,pch=ifelse(temp$plot%in%"SI",1,2))
legend("topleft",legend=unique(paste(temp$time,temp$season,temp$year,sep=" ")),fill=alpha(unique(temp$time)+1,0.5),cex=1.5,border=NA,bty="n")


### betadisper
temp<-d[d$plot%in%c("BI","SI"),]
dis <- dist(temp[,7:ncol(d)])
mod <- betadisper(dis, group=paste(temp$plot,temp$time),type="centroid")
mod
anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)

## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

## can also specify which axes to plot, ordering respected
plot(mod, axes = c(2,1), seg.col = "forestgreen", seg.lty = "dashed")

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)


### ggplot biplot pca
g<-ggbiplot(pca,choices=1:2,obs.scale=1,var.scale=1,groups=paste(temp$time),ellipse=TRUE,circle=FALSE,labels=temp$plot,labels.size=5)
g<-g+scale_color_discrete(name='YY')
g<-g+theme(legend.direction='horizontal',legend.position='top')
print(g)
#ggscreeplot(pca)



