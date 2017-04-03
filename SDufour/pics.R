library(readxl)
library(lme4)
library(visreg)
library(hglm)
library(plyr)
library(scales)
library(car)
library(AICcmodavg)
library(multcomp)
library(betareg)
library(mgcv)
library(gamm4)

### questionnement
# prendre en considération le temps s'écoulant entre les trois saisons ou les périodes de mesure?
# est-ce qu'il y a un intérêt par rapport à l'arbre individuel. La seule variable reliée à l'arbre est le DHP et ça semble plus vu comnme un contrôle...
# séparer l'analyse en deux temps avec la probabilité d'utilisation et après parmi les arbres utilisés, la superficie?
# s'intéresser à la proportion d'arbres utilisés plutôt qu'à l'arbre en tant que tel
# combien de temps ça prend pour être utilisé?

d<-read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/SDufour/Doc/BD_LONG.xlsx",na="N/A")
d$SAISON<-factor(d$SAISON,levels=c("A2015","E2016","A2016"))
d$TRAIT<-factor(d$TRAIT)
d<-d[order(d$SAISONNUM),]
d<-ddply(d,.(SEQ),function(i){
  i$AGE<-1:nrow(i)
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

d$AGE<-as.factor(d$AGE) # 1ère, 2e ou 3e mesure?
d$BATCH<-ifelse(d$SEQ>400,"batch2","batch1") # batch de tiges 1 ou 2
d$SEQ<-paste0("id",d$SEQ)
d$PREV<-as.factor(d$PREV) # la tige a été utilisée précédemment ou non?
d$TRANS<-factor(paste(as.character(d$SAISON),as.character(d$BATCH),sep="_"),levels=c("A2015_batch1","E2016_batch1","E2016_batch2","A2016_batch1","A2016_batch2"))


##############################
### GRAPHIQUES
##############################

par(mfrow=c(1,2))
plot(SUPER~SAISONNUM,data=d,col=ifelse(d$BATCH=="batch2",alpha("blue",0.3),alpha("red",0.3)),cex=0.6,xaxt="n")
axis(1,at=1:3,labels=unique(d$SAISON))
invisible(lapply(dlply(d,.(SEQ),function(j){j}),function(i){
  lines(i$SAISONNUM,i$SUPER,col=ifelse(i$BATCH=="batch2",alpha("blue",0.3),alpha("red",0.3)))
}))
plot(RELATIVE~SAISONNUM,data=d,col=ifelse(d$BATCH=="batch2",alpha("blue",0.3),alpha("red",0.3)),cex=0.6,xaxt="n")
axis(1,at=1:3,labels=unique(d$SAISON))
invisible(lapply(dlply(d,.(SEQ),function(j){j}),function(i){
  lines(i$SAISONNUM,i$RELATIVE,col=ifelse(i$BATCH=="batch2",alpha("blue",0.3),alpha("red",0.3)))
}))


########################################
### MODÈLE DE PROBABILITÉ D'UTILISATION
########################################

m<-gamm(UTIL~TRAIT+DHP+PREV_SUPER_REL+SAISON,random=list(SEQ=~1,BLOC=~1,SITE=~1),family=binomial,data=d[d$BATCH=="batch1",],correlation=corCAR1(form=~SAISONNUM|SEQ))
m<-gamm4(UTIL~TRAIT+DHP+PREV_SUPER_REL+SAISON,random=~(1|BLOC/SITE/SEQ),family=binomial,data=d)
x<-seq(10,30,by=1)
newdat<-data.frame(DHP=x,TRAIT="PATCH",PREV_SUPER_REL=0.00,SAISON="A2015")
p<-predict(m$gam,newdata=newdat,type="response")
plot(x,p,type="l",ylim=0:1)
points(x,p)

### glm
m1<-glm(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT+,family=binomial,data=d)
m2<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+AGE,family=binomial,data=d)
m3<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+TRANS,family=binomial,data=d)
#m4<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+PREV+TRANS,family=binomial,data=d)
ml<-list(m1,m2,m3)
aictab(ml)

### glmer
control<-glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))
mm1<-glmer(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=control)
mm2<-glmer(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=control)
mm3<-glmer(UTIL~TRAIT+DHP+PREV_SUPER_REL+TRANS+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=d,control=control)
mm<-gamm(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,random=list(SEQ=~1,BLOC=~1,SITE=~1),family=binomial,data=d,correlation=corCompSymm(form=~SAISON|SEQ))
mml<-list(mm1,mm2,mm3)
aictab(mml)

### plot glm effects
par(mfrow=c(2,2))
visreg(mm1,scale="response",rug=FALSE)
p<-predict(mm1,newdata=newdat,type="response")
lines(x,p)

### comparaison entre les SAISON_BATCH
summary(glht(m3, mcp(TRANS="Tukey")))


##############################################
### MODÈLES DE SUPERFICIE/PROPORTION UTILISÉE
##############################################


### MODÈLE DE REGRESSION BETA
b<-betareg(RELATIVE~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,data=d[which(d$RELATIVE>0),])
par(mfrow=c(2,3))
v<-visreg(b,rug=FALSE,ylim=0:1,las=2)
#points(model.frame(b)$DHP,model.frame(mb2)$RELATIVE)

b<-gamm(RELATIVE~TRAIT+DHP+PREV_SUPER_REL+TRANS,random=list(SEQ=~1,BLOC=~1,SITE=~1),data=d[which(d$RELATIVE>0),],family=betar(link="logit"))

hglm(y=RELATIVE,fixed=SUPER~SAISON+TRAIT+DHP+NOMESURE,random = ~ 1|SEQ,family = Gamma(link = log),rand.family = Gamma(link = log),disp = ~ SAISON + TRAIT + DHP, data = model.frame(m3))

d$BLOC<-as.factor(d$BLOC)
d$SITE<-as.factor(d$SITE)
d$SEQ<-as.factor(d$SEQ)

b1<-glmmadmb(RELATIVE~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,random=~(1|BLOC/SITE/SEQ),family="beta",data=na.omit(d[which(d$RELATIVE>0),c("RELATIVE","TRAIT","DHP","DHPREP","PM5kF","PREV_SUPER_REL","SAISON","BLOC","SITE","SEQ")]))
b2<-glmmadmb(RELATIVE~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF,random=~(1|BLOC/SITE/SEQ),family="beta",data=na.omit(d[which(d$RELATIVE>0),c("RELATIVE","TRAIT","DHP","DHPREP","PM5kF","PREV_SUPER_REL","SAISON","BLOC","SITE","SEQ")]))

glmmadmb.AICc(list(b1,b2))
model.sel(list(b1,b2))
model.avg(list(b1,b2))

par(mfrow=c(2,3))
v<-visreg(b,rug=FALSE,ylim=0:1,las=2,trans=plogis)

### MODÈLE AVEC SUPERFICIE 
# peu concluant en raison de l'incertitude sur comment bien faire les analyses 
m<-glm(SUPER~TRAIT+DHP+PREV_SUPER_REL+TRANS,data=d[d$SUPER>0,],family=gaussian(link=log))
mm<-glmer(SUPER~TRAIT+DHP+PREV_SUPER_REL+TRANS+(1|BLOC)+(1|SITE)+(1|SEQ),data=d[d$SUPER>0,],family=gaussian(link=log),control=control)
m2<-glm(SUPER~TRAIT+DHP+PREV_SUPER_REL+TRANS,data=d[d$SUPER>0,],family=Gamma(link=log))
mm2<-glmer(SUPER~TRAIT+DHP+PREV_SUPER_REL+TRANS+(1|BLOC)+(1|SITE)+(1|SEQ),data=d[d$SUPER>0,],family=Gamma(link=log),control=control)
visreg(m,rug=FALSE,scale="response")
plot(fitted(mm),resid(mm))
hist(resid(mm))
boxplot(SUPER~TRANS,data=d[d$SUPER>0,])






#idseq<-sample(unique(d$SEQ),30)
#boxplot(SUPER~SEQ,data=d[d$SEQ%in%idseq,])
#vv <- visreg(mm, "DHP", by="SEQ", re.form=~(1|SEQ), plot=FALSE)
#subSEQ <- sample(unique(d$SEQ[d$SUPER>0]), 10)
#v[[1]]<-vv[[1]][vv[[1]]$SEQ %in% subSEQ,]
#v[[2]]<-vv[[2]][vv[[2]]$SEQ %in% subSEQ,]
#v[[3]]<-vv[[3]]
#names(v)<-names(vv)
#plot(v,layout=c(10,1))







