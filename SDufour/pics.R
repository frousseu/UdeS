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
library(glmmADMB)
library(MASS)

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

x<-seq(10,30,by=1)
newdat<-data.frame(DHP=x,TRAIT="PATCH",PREV_SUPER_REL=0.00,SAISON="A2015")
p<-predict(m$gam,newdata=newdat,type="response")
plot(x,p,type="l",ylim=0:1)
points(x,p)

### glm
m1<-glm(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,family=binomial,data=d)
m2<-glm(UTIL~SAISON*TRAIT+DHP+DHPREP+SAISON+PM5kF*TRAIT,family=binomial,data=d)
#m4<-glm(UTIL~TRAIT+DHP+PREV_SUPER_REL+PREV+TRANS,family=binomial,data=d)
ml<-list(m1,m2)
aictab(ml)

### glmer

d$BLOC<-as.factor(d$BLOC)
d$SITE<-as.factor(d$SITE)
d$SEQ<-as.factor(d$SEQ)
d2<-na.omit(d[,c("UTIL","RELATIVE","TRAIT","DHP","DHPREP","PM5kF","PREV_SUPER_REL","SAISON","BLOC","SITE","SEQ")])
dd<-na.omit(d[,c("UTIL","TRAIT","DHP","DHPREP","PM5kF","PREV_SUPER_REL","SAISON","BLOC","SITE","SEQ")])

control<-glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))
mm1<-glmer(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=dd,control=control)
mm2<-glmer(UTIL~SAISON*TRAIT+DHP+DHPREP+SAISON+PM5kF*TRAIT+(1|BLOC)+(1|SITE)+(1|SEQ),family=binomial,data=dd,control=control)
mm3<-glmmadmb(UTIL~SAISON+TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF+TRAIT,random=~(1|BLOC/SITE/SEQ),family="binomial",data=d2,admb.opts=admbControl(shess=FALSE,noinit=FALSE, impSamp=200,maxfn=1000,imaxfn=500,maxph=5))
mml<-list(mm1,mm2)
aictab(mml)

### gamm et glmmPQL (ne fonctionne pas et ne fait probablement de sens)
mm<-gamm(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,random=list(SEQ=~1,BLOC=~1,SITE=~1),family=binomial,data=d,correlation=corCompSymm(form=~SAISONNUM|SEQ))
mm<-glmmPQL(UTIL~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,random=~1|BLOC/SITE/SEQ,family=binomial,data=d,correlation=corCompSymm(form=~SAISON|BLOC/SITE/SEQ))

### plot glm (or glmm) effects
# ne pas oublier de considérer les interactions avec le by argument
par(mfrow=c(2,2))
visreg(mm1,scale="response",rug=FALSE)
p<-predict(mm1,newdata=newdat,type="response")
lines(x,p)

### comparaison entre les SAISON_BATCH
summary(glht(m3, mcp(TRANS="Tukey")))


##############################################
### MODÈLES DE PROPORTION UTILISÉE
##############################################


### MODÈLE DE REGRESSION BETA SANS EFFETS ALÉATOIRES
b<-betareg(RELATIVE~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,data=d[which(d$RELATIVE>0),])
par(mfrow=c(2,3))
v<-visreg(b,rug=FALSE,ylim=0:1,las=2)
#points(model.frame(b)$DHP,model.frame(mb2)$RELATIVE)


### MODÈLE DE REGRESSION BETA AVEC ADModelBuilder et EFFETS ALÉATOIRES
# ADModelBuilder est un programme externe pour faire des GLMM qui permet de spécifier plus de distributions et de la zero-inflation
# L'installation du package installe également le .exe du programme

# glmmadmb semble préférer les facteurs et l'absence de NA
d$BLOC<-as.factor(d$BLOC)
d$SITE<-as.factor(d$SITE)
d$SEQ<-as.factor(d$SEQ)
d2<-na.omit(d[which(d$RELATIVE>0),c("RELATIVE","TRAIT","DHP","DHPREP","PM5kF","PREV_SUPER_REL","SAISON","BLOC","SITE","SEQ")])

### MODÈLES
b1<-glmmadmb(RELATIVE~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF*TRAIT,random=~(1|BLOC/SITE/SEQ),family="beta",data=d2)
b2<-glmmadmb(RELATIVE~SAISON*TRAIT+DHP+DHPREP+SAISON+PM5kF*TRAIT,random=~(1|BLOC/SITE/SEQ),family="beta",data=d2)
b3<-glmmadmb(RELATIVE~SAISON*TRAIT+DHP+DHPREP+PREV_SUPER_REL+SAISON+PM5kF,random=~(1|BLOC/SITE/SEQ),family="beta",data=d2)

### Sélection de modèle avec MuMIn, AICcmodavg n'est pas (encore) défini pour la classe glmmadmb
model.sel(list(b1,b2,b3))
model.avg(list(b1,b2,b3))

par(mfrow=c(2,3))
v<-visreg(b1,rug=FALSE,ylim=0:1,las=2,trans=plogis)

# Les prédictions sont similaires à celles d'un glm















