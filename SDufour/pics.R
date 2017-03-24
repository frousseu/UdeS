library(readxl)
library(lme4)
library(visreg)

# prendre en considération le temps s'écoulant entre les trois saisons ou les périodes de mesure
# est-ce qu'il y a un intérêt par rapport à l'arbre individuel. La seule variable reliée à l'arbre est le DHP et ça semble plus vu comnme un contrôle...
# séparer l'analyse en deux temps avec la probabilité d'utilisation et après parmi les arbres utilisés, la superficie?
# s'intéresser à la proportion d'arbres utilisés plutôt qu'à l'arbre en tant que tel

d<-read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/SDufour/Doc/BD_LONG.xlsx",na="N/A")
d$SAISON<-factor(d$SAISON,levels=c("A2015","E2016","A2016"))
d$SEQ<-paste0("id",d$SEQ)

d1<-d[d$SEQ%in%1:400,]
mm<-glmer(UTIL~SAISON*TRAIT+DHP+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d1,control=glmerControl(optCtrl=list(maxfun=20000) ))
#<-glm(UTIL~SAISON*TRAIT+DHP,family=binomial,data=d)

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


m1<-glm(UTIL~SAISON+TRAIT+DHP,data=d,family=binomial(link=logit))
m2<-glm(SUPER~SAISON+TRAIT+DHP,data=d[d$SUPER>0,],family=Gamma(link=log))
mm2<-glmer(SUPER~SAISON+TRAIT+DHP+(1|BLOC)+(1|SITE)+(1|SEQ),data=d[d$SUPER>0,],family=Gamma(link=log))
plot(fitted(m2),resid(m2))

visreg(m2,"DHP",by="SAISON",scale="response",type="conditional",overlay=TRUE)
visreg(mm2,"DHP",by="SAISON",scale="response",type="conditional",overlay=TRUE)
points(d$DHP,d$SUPER,col=as.numeric(d$SAISON)+1)


vv <- visreg(mm, "DHP", by="SEQ", re.form=~(1|SEQ), plot=FALSE)
subSEQ <- sample(unique(d$SEQ[d$SUPER>0]), 10)
v[[1]]<-vv[[1]][vv[[1]]$SEQ %in% subSEQ,]
v[[2]]<-vv[[2]][vv[[2]]$SEQ %in% subSEQ,]
v[[3]]<-vv[[3]]
names(v)<-names(vv)
plot(v,layout=c(10,1))












