library(foreign)
library(data.table)

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





