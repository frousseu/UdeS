library(foreign)
library(data.table)
library(readxl)
library(raster)
library(doParallel)
library(foreach)
library(FRutils)
library(rasterVis)

#############################################
### temprature ##############################

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/Pieges_5km/Pieges_5km"
dos<-c("PRTOT","TMOY","TMIN","TMAX") # list folders and variables
d<-vector(mode="list",length=length(dos)) # create list to store data

for(j in seq_along(dos)){
  
  x<-list.files(file.path(path,dos[j]),pattern=".dbf",full.names=TRUE) # liste les fichiers dans le chemin
  l<-sapply(x,function(i){read.dbf(i)}) # lit tous les fichiers et les sauvegarde dans une liste

  g<-grep("030101",names(l[[1]])) # Ceci cherche le premier janvier et toutes les colonnes précédentes sont répétées
  keep<-names(l[[1]])[1:(g-1)] # Ce sont les colonnes à conserver et répéter (le buffer de 5km est sûrement à enlever)

  ld<-lapply(l,function(i){
    res<-melt(as.data.table(i),keep,variable.name="id") # put in long form
    res$id<-gsub("P|TY|TN|TX","",res$id) # eliminate prefix from date data
    names(res)[ncol(res)]<-dos[j] # assign variable name
    res
  })
  d[[j]]<-rbindlist(ld) # ramène tous dans un même data.frame à partir de la liste ld
  
}

d<-Reduce(function(...){merge(...,all=TRUE)},d) # merge four variables in a single data.frame
temp<-strsplit(d$id, "(?<=.{2})", perl = TRUE) # split date info to construct an ok date
d$date<-sapply(temp,function(i){paste(paste0(c("20","",""),i),collapse="-")}) # construct date as a character


### graph to check the time series
s<-split(d,d$Site_Seq)[1:10] # choose the first 10 sites for faster plotting
plot(as.Date(s[[1]]$date),rep(0,nrow(s[[1]])),ylim=c(-35,35),type="n")
invisible(lapply(s,function(i){
   lines(as.Date(i$date),i$TMOY) 
}))


##########################################
### epiweeks #############################

# dans le fichier epiweeks, on n'a aucune valeur de week pour le 2004-02-29, 2008-12-31, 2012-12-31

w<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/Pieges_5km/Pieges_5km/EpiWeeks_calculation.xls",skip=5,col_types="text"))
w<-melt(w,measure.vars=1:ncol(w))
w<-split(w,rep(1:3,each=nrow(w)/3))
w<-do.call("cbind",w)
w<-w[,grep("value",names(w))]
names(w)<-c("date","cdcweek","cdcweekcum")
w$cdcweek<-as.integer(w$cdcweek)
w$cdcweekcum<-as.integer(w$cdcweekcum)
w$date<-as.character(as.Date(as.integer(w$date),origin="1899-12-30")) # usually the origin is 1970-01-01, make sure the dates are ok when back in R, excel might treat them otherwise

###########################################
### full data with temp and cdc weeks #####
x<-merge(d,as.data.table(w),all.x=TRUE,by="date")

x<-x[!is.na(x$cdcweek)] # on élimine les lignes sans cdcweek (voir commentaire plus haut, penser aux conséquences sur les calculs basés sur les weeks)

###########################################
### LULC ##################################

### ceci devra être modifié pour accomoder plusieurs buffers

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc"
s<-read.dbf(file.path(path,"LULC2011_2km.dbf"))
vars<-names(s)[grep("HIST",names(s))]
res<-s[,vars]/apply(s[,vars],1,sum)
names(res)<-paste0(names(res),"p2km")
s<-cbind(s,res)

x<-merge(x,s[,setdiff(names(s),c(vars,"Buff"))],all.x=TRUE,by=c("Site_Seq","CodeSite"))

###########################################
### compute values ########################

### temp hebdo

lv<-list()

com<-c("Site_Seq","CodeSite","cdcweekcum")
temp<-c("TMOY","TMIN","TMAX")
lv[[1]]<-setnames(x[,lapply(.SD,mean),by=com,.SDcols=temp],temp,paste0(temp,"week"))
lv[[length(lv)+1]]<-setnames(x[,lapply(.SD,sum),by=com,.SDcols=c("PRTOT")],"PRTOT","PRTOTweek")

tday<-function(x,val,sign){
  if(sign){
    sum(x>val)
  }else{
    sum(x<val)  
  }
}

lv[[length(lv)+1]]<-setnames(x[,lapply(.SD,tday,val=32,sign=1),by=com,.SDcols=c("TMAX")],"TMAX","daysover32")
lv[[length(lv)+1]]<-setnames(x[,lapply(.SD,tday,val=-5,sign=0),by=com,.SDcols=c("TMIN")],"TMIN","daysbelow5")

### gdd

# on le veut par jour ou par semaine?

x$TMAXTMIN<-x$TMIN+x$TMAX

gdd<-function(x,tseuil){ # pas sûr de cette formule, on fait la somme sur la semaine?
  res<-(x/2)-tseuil
  sum(ifelse(res<0,0,res))
}

temp<-c(9,12,18,20,25)

for(i in seq_along(temp)){
  lv[[length(lv)+1]]<-setnames(x[,lapply(.SD,gdd,tseuil=temp[i]),by=com,.SDcols=c("TMAXTMIN")],"TMAXTMIN",paste0("gdd",temp[i]))
}

### bind everything
xx<-do.call("cbind",lapply(lv,function(i){
  i[,setdiff(names(i),com),with=FALSE]  
}))
xx<-cbind(lv[[1]][,com,with=FALSE],xx)

### merge with big database
x<-merge(x,xx,by=com)

### get only weekly data
xx<-x[!duplicated(paste(x$Site_Seq,x$cdcweekcum)),setdiff(names(x),c("date","id","Buff","PRTOT","TMOY","TMIN","TMAX")),with=FALSE]



###############################################
### MOSQUITO DATA #############################

m<-as.data.table(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/BD.xlsx"))
#setnames(m,c("Site","Day","Week"),c("CodeSite","date","cdcweek"))
#m$date<-substr(m$date,1,10) # check if the dates are ok since they are considered UTC by read_excel
#m<-merge(m,w,by=c("cdcweekcum"))
colSums(m[,names(m)[substr(names(m),1,2)%in%paste0("A",1:10)],with=FALSE])

m<-m[m$cdcweekcum<=679,] # meteo data not older

b<-intersect(names(xx),names(m))
m<-merge(xx,m[,c(b,"A9"),with=FALSE],by=b,all=TRUE)
#m$A29<-ifelse(is.na(m$A29),0,m$A29)

m<-m[order(m$CodeSite,m$cdcweekcum),]


lm<-split(m,m$CodeSite)

# clarify what is a season and remove lines for which the trap was not there


########################

#lm<-lapply(lm,function(i){
#  obs<-i$cdcweek[!is.na(i$A29)]
#  r<-range(obs,na.rm=TRUE)
#  i$season<-as.integer(i$cdcweek>=r[1] & i$cdcweek<=r[2])
#  i
#})



#n<-5
#d<-as.data.table(data.frame(v=1:n,r=runif(n)))


corlag<-function(d,v="v",r="r",lag=3){
  M<-expand.grid(0:lag,0:lag)
  m<-M[M[,1]>=M[,2],]
  me<-lapply(1:nrow(m),function(i){
    res<-Map(":",(1:nrow(d))-m[i,1],(1:(nrow(d))-m[i,2]))
    w<-which(!is.na(d[[r]]))
    ans<-sapply(w,function(j){
      mean(d[[v]][res[[j]]])
    })
    m<-rep(NA,nrow(d))
    m[w]<-ans
    m
  })
  co<-sapply(1:length(me),function(i){
    cor(me[[i]],d[[r]],use="complete.obs",method="spearman")  
  })
  ans<-cbind(m,co)
  ans<-merge(M,ans,all.x=TRUE)
  ans<-matrix(ans$co,nrow=lag+1,ncol=lag+1,byrow=TRUE)
  dimnames(ans)[[1]]<-0:(lag)
  dimnames(ans)[[2]]<-0:(lag)
  ans
}


lm<-lm[sapply(lm,function(i){!all(i$A9%in%c(0,NA))})] # remove traps with no observations

### add temperature normal
n<-lapply(lm,function(i){
  g<-gam(TMOYweek~s(cdcweek),data=i)
  val<-0:53
  p<-predict(g,data.frame(cdcweek=val))
  p[match(i$cdcweek,val)]
})

lm<-Map("cbind",lm,n)
lm<-lapply(lm,function(i){
  i$diff<-i$TMOYweek-i$V2    
  i
})


registerDoParallel(6) 
getDoParWorkers()

set.seed(1234)
lag<-20
do<-sample(1:length(lm),12)
#do<-261:262
cm<-foreach(i=do,.packages=c("data.table")) %dopar% {
  corlag(lm[[i]],v="diff",r="A9",lag=lag)
}



#cmat<-Reduce("+",cm)/length(cm)

r<-stack(lapply(cm,raster,xmn=0,xmx=lag,ymn=0,ymx=lag))
names(r)<-names(lm)[do]
#levelplot(r,col.regions=rasterTheme()$regions$col,cuts=99)
levelplot(r,col.regions=colo.scale(1:100,c("darkred","red","lightgoldenrod","blue","navyblue")),cuts=99)
plot(mean(r,na.rm=TRUE),col.regions=colo.scale(1:100,c("darkred","red","lightgoldenrod","blue","navyblue")),cuts=99)
#plot(r,axes=FALSE)
#axis(1,at=0:nrow(cmat)+0.5,labels=0:nrow(cmat))
#axis(2,at=0:nrow(cmat)-0.5,labels=rev(0:nrow(cmat)))

k<-262
plot(lm[[k]]$cdcweek,lm[[k]]$TMOYweek)
points(lm[[k]]$cdcweek,log(1+lm[[k]]$A9),pch=16,col="darkgreen",cex=1.5)
points(lm[[k]]$cdcweek[!is.na(lm[[k]]$A9)],lm[[k]]$TMOYweek[!is.na(lm[[k]]$A9)],col="red",pch=16)
lines(lm[[k]]$cdcweek,lm[[k]]$V2,pch=16)
#points(lm[[k]]$cdcweek,lm[[k]]$diff,pch=16)


k<-1:
dat<-rbindlist(lm[k])
g<-gam(A9~s(cdcweek,bs="cc")+s(diff),data=dat)
val<-15:45
p<-predict(g,data.frame(cdcweek=val),type="response")
plot(dat$cdcweek,dat$A9,ylim=c(0,10))
lines(val,p,type="l")







