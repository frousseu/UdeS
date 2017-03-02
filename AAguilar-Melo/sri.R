
library(lme4)
library(visreg)
library(plyr)

## c,est quoi season_per1, les id des individus et les dyads, variable même groupe ou pas?
## il semble qu'il faut vraiment scaler les données, sinon on se retrouve avec plein de warnings de non-convergence
## check asnipe package, it computes sri values, but does not keep numerator and denominator

###########################
### get data
###########################

#a<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLMOTHERcortada.csv")
PL<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/PLADRY-b.csv", colClasses= "character")
d<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/AAguilar/Doc/SRIModels_PA_Ago16_2.csv")
d<-d[!is.na(d$sri),]


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


#####################################
### run model
#####################################

# join dyad sri data with d data
d<-join(d,x,type="left")

# give simpler names toi variable and scale them
d$ficus<-d$ifa.ficus
d$brosimum<-d$ifa.brosimum
d$var<-d$variance.ft

# scale variable, otherwise model does not converge
d$ficus.s<-scale(d$ficus)[,1]
d$brosimum.s<-scale(d$brosimum)[,1]
d$var.s<-scale(d$var)[,1]


m<-glmer(cbind(ab,abnot)~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s+(1|dyad1)+(1|season_per1),data=d,family=binomial)

m2<-glm(cbind(ab,abnot)~season*ficus.s+season*brosimum.s+season*var.s+ficus.s*brosimum.s+ficus.s*var.s+brosimum.s*var.s,data=d,family=binomial)

visreg(m,"brosimum.s",by="ficus.s",scale="response",type="conditional",overlay=FALSE)
