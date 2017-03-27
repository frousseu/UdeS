library(readxl)
library(lme4)
library(visreg)
library(hglm)
library(plyr)
library(scales)
library(car)
library(AICcmodavg)


# prendre en considération le temps s'écoulant entre les trois saisons ou les périodes de mesure
# est-ce qu'il y a un intérêt par rapport à l'arbre individuel. La seule variable reliée à l'arbre est le DHP et ça semble plus vu comnme un contrôle...
# séparer l'analyse en deux temps avec la probabilité d'utilisation et après parmi les arbres utilisés, la superficie?
# s'intéresser à la proportion d'arbres utilisés plutôt qu'à l'arbre en tant que tel
# combien de temps ça prend pour être utilisé

d<-read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/SDufour/Doc/BD_LONG.xlsx",na="N/A")
d$SAISON<-factor(d$SAISON,levels=c("A2015","E2016","A2016"))
d<-d[order(d$SAISONNUM),]
#d$NBSEASON<-ifelse(d$SEQ>400,2,3)
d<-ddply(d,.(SEQ),function(i){
  i$NOMESURE<-1:nrow(i)
  i$NBSAISON<-nrow(i)
  w<-which(i$UTIL==1)
  if(any(w) && min(w)<nrow(i)){
    i$PREV<-c(0,i$UTIL[-nrow(i)])
    m<-match(1,i$PREV)
    i$PREV[m:nrow(i)]<-1
  }else{
    i$PREV<-rep(0,nrow(i))
  }
  i
})

d<-transform(d,PREV_SUPER=ave(SUPER,SEQ,FUN=function(i){cumsum(c(0,head(i,-1)))}))
d$PREV_SUPER_REL<-d$PREV_SUPER/((d$DHP/2)*pi*2*200)

d$NOMESURE<-as.factor(d$NOMESURE)
d$NBSAISON<-as.factor(d$NBSAISON)
d$PREV<-as.factor(d$PREV)


x<-ddply(d[!is.na(d$UTIL),],.(NBSAISON,SAISON,UTIL),nrow)
x<-x[order(x$SAISON,x$NBSAISON,x$UTIL),]

par(mfrow=c(1,2))
plot(SUPER~SAISONNUM,data=d,col=ifelse(d$SEQ>400,alpha("blue",0.3),alpha("red",0.3)),cex=0.6,xaxt="n")
axis(1,at=1:3,labels=unique(d$SAISON))
lapply(dlply(d,.(SEQ),function(j){j}),function(i){
  lines(i$SAISONNUM,i$SUPER,col=ifelse(i$SEQ>400,alpha("blue",0.3),alpha("red",0.3)))
})
plot(RELATIVE~SAISONNUM,data=d,col=ifelse(d$SEQ>400,alpha("blue",0.3),alpha("red",0.3)),cex=0.6,xaxt="n")
axis(1,at=1:3,labels=unique(d$SAISON))
lapply(dlply(d,.(SEQ),function(j){j}),function(i){
  lines(i$SAISONNUM,i$RELATIVE,col=ifelse(i$SEQ>400,alpha("blue",0.3),alpha("red",0.3)))
})


d$SEQ<-paste0("id",d$SEQ)


#mm1<-glmer(UTIL~TRAIT+DHP+SAISON+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=glmerControl(optCtrl=list(maxfun=20000) ))
#mm2<-glmer(UTIL~TRAIT+DHP+NOMESURE+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=glmerControl(optCtrl=list(maxfun=20000) ))
#mm3<-glmer(UTIL~TRAIT+DHP+SAISON+NOMESURE+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=glmerControl(optCtrl=list(maxfun=20000) ))
#mml<-list(mm1,mm2,mm3)
#aictab(mml)

m1<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+SAISON,family=binomial,data=d)
m2<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+NOMESURE,family=binomial,data=d)
m3<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+SAISON+NOMESURE,family=binomial,data=d)
m4<-glm(UTIL~TRAIT+DHP+SAISON+NOMESURE,family=binomial,data=d)

ml<-list(m1,m2,m3,m4)
aictab(ml)

visreg(m3,"PREV_SUPER_REL",scale="response",rug=FALSE,cond=list(SAISON="A2016"))
visreg(m3,"SAISON",scale="response",rug=FALSE,cond=list(NOMESURE=2))
visreg(m3,"NOMESURE",scale="response",rug=FALSE)
visreg(m3,"TRAIT",scale="response",rug=FALSE,cond=list(NOMESURE=3,SAISON="E2016"))

mm2<-lmer(log(SUPER)~SAISON+TRAIT+DHP+(1|BLOC/SEQ),data=d[d$SUPER>0,])
visreg(mm2,"SAISON",by="TRAIT",trans=exp,rug=FALSE)
plot(fitted(mm2),resid(mm2))
hist(resid(mm2))

d2<-d[d$SEQ>400,]
mm<-glmer(UTIL~SAISON+TRAIT+DHP+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=glmerControl(optCtrl=list(maxfun=100000) ))
#m<-glm(UTIL~SAISON*TRAIT+DHP,family=binomial,data=d2)
visreg(mm,"SAISON",by="TRAIT",scale="response",type="conditional",rug=FALSE)

boxplot(log(SUPER)~SAISON+TRAIT,data=d[d$SUPER>0,],las=2)


idseq<-unique(d$SEQ[which(d$SUPER>0)])
dpos<-d[d$SEQ%in%idseq,]

mm<-lmer(log(SUPER+1)~SAISON+TRAIT+DHP+(1|BLOC)+(1|SITE)+(1|SEQ),data=dpos)
plot(fitted(mm),resid(mm))


newdat<-data.frame(SAISON="E2015",TRAIT="PATCH")

m1<-glm(UTIL~SAISON+TRAIT+DHP,data=d,family=binomial(link=logit))
m2<-glm(SUPER~SAISON+TRAIT+DHP,data=d[d$SUPER>0,],family=Gamma(link=log))


mm2<-glmer(SUPER~SAISON+TRAIT+DHP+(1|BLOC)+(1|SITE)+(1|SEQ),data=d[d$SUPER>0,],family=Gamma(link=log))
plot(fitted(mm2),resid(mm2))

visreg(m2,"DHP",by="SAISON",scale="response",type="conditional",overlay=TRUE)
visreg(mm2,"DHP",by="SAISON",scale="response",type="conditional",overlay=TRUE)
points(d$DHP,d$SUPER,col=as.numeric(d$SAISON)+1,pch)



idseq<-sample(unique(d$SEQ),30)
boxplot(SUPER~SEQ,data=d[d$SEQ%in%idseq,])



vv <- visreg(mm, "DHP", by="SEQ", re.form=~(1|SEQ), plot=FALSE)
subSEQ <- sample(unique(d$SEQ[d$SUPER>0]), 10)
v[[1]]<-vv[[1]][vv[[1]]$SEQ %in% subSEQ,]
v[[2]]<-vv[[2]][vv[[2]]$SEQ %in% subSEQ,]
v[[3]]<-vv[[3]]
names(v)<-names(vv)
plot(v,layout=c(10,1))



m11<-hglm(fixed=SUPER~SAISON+TRAIT+DHP+NOMESURE,random = ~ 1|SEQ,family = Gamma(link = log),rand.family = Gamma(link = log),disp = ~ SAISON + TRAIT + DHP, data = model.frame(m3))










