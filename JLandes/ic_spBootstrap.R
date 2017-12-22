
library(icenReg)
library(foreach)
library(doParallel)

#################################################################
### Build example model from ic_sp help 
#################################################################

set.seed(1)

d<-simIC_weib(n=500,inspections=5,inspectLength=1) # données simulées
d$id<-as.factor(rep(1:(nrow(d)/10),length.out=nrow(d)))
myCluster<-makeCluster(4) # put number of cores
registerDoParallel(myCluster)
m<-ic_sp(Surv(l,u,type='interval2')~x1+x2,data=d,useMCores=TRUE,bs_samples=100) # un modèle à updater avec de nouvelles données
stopCluster(myCluster)

plot(m)

l<-getSCurves(m)$Tbull_ints[,"lower"]
u<-getSCurves(m)$Tbull_ints[,"upper"]
lines(l,getSCurves(m)$S_curves$baseline,col="red",type="s")
lines(u,getSCurves(m)$S_curves$baseline,col="green",type="s")

# m = model
# d = data
# d$id = frailty 

####################################################
### function to update an ic_sp model with new data
####################################################


update_icenReg<-function(object,...){
  call <- object$call
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0L) {
    ic_spa <- names(as.list(args(ic_sp)))
    names(extras) <- ic_spa[pmatch(names(extras), ic_spa[-length(ic_spa)])]
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

#####################################################
### run models
#####################################################

### bootstrap samples with the different ids
n<-100 # number of replication
ids<-unique(d$id) # list of unique frailties

w<-lapply(ids,function(i){which(d$id==i)})
names(w)<-ids

l<-lapply(seq_len(n),function(i){
  s<-sample(ids,length(ids),replace=TRUE)
  m<-match(s,names(w))
  d[unlist(w[m]),]
})

### run a model for each sample
registerDoParallel(4) #put number of cores
getDoParWorkers()

cl<-foreach(i=seq_along(l),.packages=c("icenReg")) %dopar% {
  newm<-update_icenReg(m,data=l[[i]],bs_samples=0)
  coef(newm) 
}

pr<-foreach(i=seq_along(l),.packages=c("icenReg")) %dopar% {
  newm<-update_icenReg(m,data=l[[i]],bs_samples=0)
  lo<-getSCurves(newm)$Tbull_ints[,"lower"]
  up<-getSCurves(newm)$Tbull_ints[,"upper"]
  ba<-getSCurves(newm)$S_curves$baseline
  cbind(lo,up,ba)
}

plot(m)
lapply(pr,function(i){
  lines(i[,"lo"],i[,"ba"],col=alpha("red",0.03),type="s",lwd=3)
  lines(i[,"up"],i[,"ba"],col=alpha("blue",0.03),type="s",lwd=3)
})

par(new=TRUE)
plot(m)


cl<-do.call("rbind",cl) # bind list of coefficients

co<-t(apply(cl,2,function(i){
  c(mean=mean(i),lowCI=quantile(i,0.025),uppCI=quantile(i,0.975),se=sd(i))  
})) # summary of coefficients

co
m





