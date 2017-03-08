
library(lme4)
library(visreg)
library(plyr)
library(readxl)

## c,est quoi season_per1, les id des individus et les dyads, variable même groupe ou pas?
## il semble qu'il faut vraiment scaler les données, sinon on se retrouve avec plein de warnings de non-convergence
## check asnipe package, it computes sri values, but does not keep numerator and denominator
## need to determine which periods are used, because the number of observations in the dyad calculation does not match with the number of observations analysed

###########################
### get data
###########################

#a<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLMOTHERcortada.csv")
PL<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLADRY-b.csv", colClasses= "character")
PL<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLADRY-b.csv", colClasses= "character")
d<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/SRIModels_PA_Ago16_2.csv",stringsAsFactors=FALSE)
#d<-d[!is.na(d$sri),]

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


######################################
### verifications and binding all data
######################################

# bind all data
x<-rbind(x,r[,names(x)])
# join dyad sri data with d data
d<-join(d,x,type="left")
# compare old sri to new sri2
# some values not matching in three periods (Nov09-1,Nov13_2,Jan14_1)
View(d[which(abs(d$sri-d$sri2)>0.0001),])
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
View(d)

ddply(d,.(dyad1),function(i){all(i$ab==0)})

#####################################
### run model
#####################################

# for now, only eliminate based on sri2
d<-d[!is.na(d$sri2),]

# give simpler names toi variable and scale them
d$ficus<-d$ifa.ficus
d$brosimum<-d$ifa.brosimum
d$var<-d$variance.ft

# scale variable, otherwise model does not converge
d$ficus.s<-scale(d$ficus)[,1]
d$brosimum.s<-scale(d$brosimum)[,1]
d$var.s<-scale(d$var)[,1]

#m<-glmer(cbind(ab,abnot)~season*ficus+season*brosimum+season*var+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s+(1|dyad1)+(1|season_per1),data=d,family=binomial)

m<-glmer(cbind(ab,abnot)~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s+(1|dyad1)+(1|season_per1),data=d,family=binomial)

#ml<-lmer(sri~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s+(1|dyad1)+(1|season_per1),data=d)

#m2<-glm(cbind(ab,abnot)~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s,data=d,family=binomial)

visreg(m,"brosimum.s",by="ficus.s",scale="response",type="conditional",overlay=TRUE,gg=TRUE)

#visreg(ml,"brosimum.s",by="ficus.s",scale="response",type="conditional",overlay=TRUE)
visreg2d(m,"brosimum.s","ficus.s",type="conditional",scale="response")