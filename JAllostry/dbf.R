library(foreign)
library(data.table)

#path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/Tmin_Tmax_25km"
path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/Doc/Pieges_5km/Pieges_5km/TMOY" # new files
x<-list.files(path,pattern=".dbf") # liste les fichiers dans le chemin
#x<-x[-grep(".pox",x)] # enlève les fichiers .pox
#x<-rep(x,10)
#x<-x[1:5]
l<-sapply(x,function(i){read.dbf(file.path(path,i))}) # lit tous les fichiers et les sauvegarde dans une liste

ld<-lapply(l,function(i){
  res<-i[rep(1:nrow(i),ncol(i)-3),c(1:3)] # répète les infos communes
  res$id<-rep(names(i)[-(1:3)],each=nrow(i))
  res$temp<-unlist(i[,4:ncol(i)],use.names=FALSE)
  #res$temp<-.Internal(unlist(i[,6:ncol(i)], FALSE, FALSE))
  #ddply(i,.(OBJECTID),function(j){ # décompose tous les fichiers
    #res<-j[rep(1,ncol(j)-5),1:5] # répète les infos communes
    #res$Temp<-j[1,6:ncol(j),drop=TRUE] # extrait les températures
    #res$Type<-substr(names(j)[-(1:5)],1,2) # extrait les dates et types de mesures
    #res$Year<-substr(names(j)[-(1:5)],3,4)
    #res$Month<-substr(names(j)[-(1:5)],5,6)
    #res$Day<-substr(names(j)[-(1:5)],7,8)
  res
})
d<-rbindlist(ld) # ramène tous dans un même data.frame à partir de la liste ld

s<-split(d,d$Site_Seq)
plot(1:365,ylim=c(-30,25),type="n")
lapply(s,function(i){
   lines(i$temp,col=gray(0,0.1)) 
})


#fwrite(d,"C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/GitHub/temp.csv")
#d<-fread("C:/Users/rouf1703/Documents/UdeS/Consultation/JAllostry/GitHub/temp.csv")


#d$type<-substr(d$id,1,2) # extrait les dates et types de mesures
#d$year<-substr(d$id,3,4)
#d$month<-substr(d$id,5,6)
#d$day<-substr(d$id,7,8)

s<-d$id
x<-strsplit(s, "(?<=.{2})", perl = TRUE)

d$type<-sapply(x,"[",1)
d$year<-sapply(x,"[",2)
d$month<-sapply(x,"[",3)
d$day<-sapply(x,"[",4)

