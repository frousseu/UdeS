
library(lme4)
library(visreg)
library(plyr)
library(readxl)
library(reshape)

## c'est quoi season_per1, les id des individus et les dyads, variable même groupe ou pas?
## il semble qu'il faut vraiment scaler les données, sinon on se retrouve avec plein de warnings de non-convergence
## voir package asnipe qui calcule des valeurs de sri, mais qui ne conserve pas les valeurs des numérateurs et dénominateurs

###########################
### get data
###########################

#a<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLMOTHERcortada.csv")
PL<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLADRY-b.csv", colClasses= "character")
d<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/SRIModels_PA_Ago16_2.csv",stringsAsFactors=FALSE)
#d<-d[!is.na(d$sri),]

# convert values to 1 or NAs
ind<-names(PL)[which(nchar(names(PL))==2)]
PL[ind]<-lapply(PL[ind],function(i){ifelse(!i%in%c(1),NA,1)}) #transforms everything not 1 in NAs (including x)

# compute % of times observed in each scan
nbobs<-ddply(PL,.(period2),function(i){
 colSums(i[,ind],na.rm=TRUE)/nrow(i) 
})
nbobs<-melt(nbobs)
names(nbobs)[2]<-"id"


########################################
### read Braulio data
########################################

r<-read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/association_period 2009-2010-PourFrancois.xlsx",sheet=2,skip=4)
r[]<-lapply(r,function(i){if(all(is.na(i))){NULL}else{i}})
n<-unique(names(r))
n<-n[!n%in%c("",NA,"NA","TOTAL")]
names(r)<-c("id1","id2",paste(rep(n,each=3),rep(letters[1:3],times=length(unique(n))),sep="_"))
r<-r[-1,]
r<-lapply(n,function(i){
  x<-r[,c("id1","id2",names(r)[grep(i,names(r))])]
  names(x)[3:5]<-letters[1:3]
  x$period2<-i
  x
})
r<-do.call("rbind",r)
r$ab<-as.integer(r$c)
r$abnot<-as.integer(r$a)+as.integer(r$b)
r$sri<-r$ab/(r$ab+r$abnot)
r$sri<-ifelse(is.nan(r$sri),NA,r$sri)
r$dyad1<-paste(r$id1,r$id2,sep="_")
r$sri2<-r$sri
### eliminate cases by hand
# r is ordered according to the different periods and the elimination here depnds on this ordering

### commentaires par adriana ###
# "AI" était après May09_1. Alors, avant il serait NA.
# "BO" la même chose.
# "HI" a émigré en avril 19 du 2010. Donc NA avant "April10_1"
# "SR" serait éliminée  car elle était juvenile en 2009 et elle a émigré en janvier 2010.
### les deux premiers on n'a pas de May09_1 les données commencent après

# CLARIFIER QUOI FAIRE AVEC AI. AI ne sembla pas être présent dans May09_1
idpos<-(r$id1=="AI" | r$id2=="AI")
periodpos<-r$period2=="May09_1"
which(idpos & periodpos)
#r<-r[-w,]

idpos<-which(r$id1=="HI" | r$id2=="HI")
periodpos<-which(r$period2=="April10_1")
w<-idpos[idpos<min(periodpos)]
r<-r[-w,]

idpos<-which(r$id1=="SR" | r$id2=="SR")
r<-r[-w,]
r$value1<-NA # used later for verifications
r$value2<-NA # used later for verifications


########################
### create sri function 
########################
# d is an observation matrix and keep numerator and denominator values ({ab,abnot} for a and b together, a or b not together)
# e.g. d<-data.frame(a=c(0,0,0,0,1),b=c(1,1,1,0,1),c=c(0,1,1,1,1),d=c(1,0,0,0,1),e=c(0,0,0,1,1))
sri<-function(d){
  x<-setNames(data.frame(do.call("rbind",combn(names(d),2,simplify=FALSE)),stringsAsFactors=FALSE),c("id1","id2"))
  x[,c("ab","abnot")]<-t(apply(x,1,function(i){
    ab<-sum(rowSums(d[,c(i[1],i[2])])==2)
    abnot<-sum(rowSums(d[,c(i[1],i[2])])==1)
    cbind(ab,abnot)
  }))[,1:2]
  x$sri<-x$ab/(x$ab+x$abnot)  
  x$dyad1<-paste(x$id1,x$id2,sep="_")
  x
}

#####################################
### compute sri index for each period
#####################################
x<-ddply(PL,.(period2),function(i){
  m<-sapply(i[,match("AE",names(i)):match("VI",names(i))],as.integer)
  m[is.na(m)]<-0L
  m<-as.data.frame(m,stringsAsFactors=FALSE)
  sri(m)
})
x$sri2<-x$sri #second sri index to make sure the new is correctly computed
x<-x[,-match("sri",names(x))] #take out the sri one to make sure both are different columns

# build % of times observed for each individuals
nbobs1<-setNames(nbobs,paste0(names(nbobs),c("",1,1)))
nbobs2<-setNames(nbobs,paste0(names(nbobs),c("",2,2)))
x<-join(x,nbobs1,type="left")
x<-join(x,nbobs2,type="left")


######################################
### verifications and binding all data
######################################

# bind all data
x<-rbind(x,r[,names(x)])

# join dyad sri data with d data
d<-join(d,x,type="full")

### VERIF
#PL[PL$period2=="Dec14_1",c("TL","VE")]
#x[x$period2=="Dec14_1" & x$dyad1%in%c("TL_VE","VE_TL"),]


# compare old sri to new sri2
# some values not matching in three periods (Nov09-1,Nov13_2,Jan14_1)
#View(d[which(abs(d$sri-d$sri2)>0.0001),])
# verify
temp<-d[,c("period2","dyad1")]
nrow(temp)
temp<-unique(temp)
nrow(temp)
temp<-setNames(data.frame(temp$period2,do.call("rbind",lapply(temp$dyad1,function(i){sort(unlist(strsplit(i,"_")))})),stringsAsFactors=FALSE),c("period2","id1","id2"))
temp<-temp[,c("period2","id1","id2")]
nrow(temp)
nrow(unique(temp))
table(is.na(d$sri))
table(is.na(d$sri2))
temp<-d
temp$sri<-temp$sri==0
temp$sri2<-temp$sri2==0
temp$sum<-sum(temp$ab,temp$abnot)
ddply(temp,.(sri,sri2),nrow)
d$sri[which((d$ab+d$abnot)==0)]

#View(d)

#ddply(d,.(dyad1),function(i){all(i$ab==0)})

#####################################
### run model
#####################################

# sri is the old sri
# sri2 is the new sri

# for now, only eliminate based on sri2
d<-d[!is.na(d$sri2),]
# eliminate dyads with individual seen less than 10% of each period scans
d<-d[-which(d$value1<0.1 | d$value2<0.1),] #take out is temp


# give simpler names to variable and scale them
d$ficus<-d$ifa.ficus
d$brosimum<-d$ifa.brosimum
d$var<-d$variance.ft

# scale variable, otherwise model does not converge
d$ficus.s<-scale(d$ficus)[,1]
d$brosimum.s<-scale(d$brosimum)[,1]
d$var.s<-scale(d$var)[,1]

#m<-glmer(cbind(ab,abnot)~season*ficus+season*brosimum+season*var+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s+(1|dyad1)+(1|season_per1),data=d,family=binomial)

m<-glmer(cbind(ab,abnot)~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s+(1|dyad1)+(1|season_per1),data=d,family=binomial)
#m2<-glm(cbind(ab,abnot)~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s,data=d,family=binomial)

visreg(m,"brosimum.s",by="season",type="conditional",overlay=TRUE,scale="response")
visreg2d(m,"brosimum.s","ficus.s",type="conditional",scale="response")
