library(foreign)
library(data.table)

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/Pieges_5km/Pieges_5km/TMOY" # new files
x<-list.files(path,pattern=".dbf") # liste les fichiers dans le chemin
l<-sapply(x,function(i){read.dbf(file.path(path,i))}) # lit tous les fichiers et les sauvegarde dans une liste


ld<-lapply(l,function(i){
  res<-i[rep(1:nrow(i),ncol(i)-3),c(1:3)] # répète les infos communes
  res$id<-rep(names(i)[-(1:3)],each=nrow(i))
  res$temp<-unlist(i[,4:ncol(i)],use.names=FALSE)
  res
})
d<-rbindlist(ld) # ramène tous dans un même data.frame à partir de la liste ld


### graph
s<-split(d,d$Site_Seq)
plot(1:365,ylim=c(-30,25),type="n")
lapply(s,function(i){
   lines(i$temp,col=gray(0,0.1)) 
})


### splitting info
s<-d$id
x<-strsplit(s, "(?<=.{2})", perl = TRUE)

d$type<-sapply(x,"[",1)
d$year<-sapply(x,"[",2)
d$month<-sapply(x,"[",3)
d$day<-sapply(x,"[",4)

