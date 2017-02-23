
################################################################
### distance calculation by pairs
################################################################

library(zipcode)
library(sp)
library(rgeos)
library(leaflet)

### get example dummy data and take a small subset
data(zipcode)
n<-10
x<-zipcode[sample(1:nrow(zipcode),n),]
x<-x[!is.na(x$latitude) & !is.na(x$longitude),]


### convert to Spatial
coordinates(x)<-~longitude+latitude
proj4string(x)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### compute great circle distance
dist<-spDists(x, longlat=TRUE)

### add names to response matrix
rownames(d)<-x@data[,1]
colnames(d)<-x@data[,1]

### optionally plot the Spatial object with leaflet
#leaflet() %>% 
#  addTiles() %>% 
#  addCircleMarkers(data=x)


###################################################################################
### MORE COMPLEX BIRTH ORDERS #####################################################
###################################################################################
library(plyr)
# should make things faster with data.table

load("~/UdeS/Consultation/SEngelhardt/Doc/PRDH.RData")

# get simpler data.frame names
g<-PRDH
rm(PRDH)

# convert factors to characters for easier manipulations
g[]<-lapply(g, function(i){if(is.factor(i)){as.character(i)}else{i}})

# get the complete list of females
keep<-c("id","sex")
x<-unique(g[,keep])
# SOME individuals have more than one sex or have also a missing sex
x[duplicated(x$id),]
f<-unique(x$id[x$sex%in%c("f")])
#are all females in the list of mothers? NO
table(f%in%c(g$dam))


### for each individual in column dam in both databases, get a list of all childs
# assumes that all individuals that are mothers are found in the dam column in both datasets
temp<-dlply(g,.(dam),function(i){
  unique(i$id)
})
# take out dam NA values
temp<-temp[-which(names(temp)%in%c("NA",NA))]

### complete list of mothers-child pairs
p<-data.frame(mother=rep(names(temp),sapply(temp,length)),child=unlist(temp,use.names=FALSE),stringsAsFactors=FALSE)
p<-p[order(as.integer(p$mother),as.integer(p$child)),]



# create a function that gets the sex, year and month at birth for each individual for a specified generation
# the function first looks in the mentions3 database (m) and when the info is missing, it then looks in the PRDH database (g)
# nothing is done to check if there is a difference between both information
# x is the data.frame where the info should be added, gen is the generation for which we want information
getinfo<-function(x,gen="F1"){
  sex<-paste0("sex",gen)
  year<-paste0("byear",gen)
  month<-paste0("bmonth",gen)
  day<-paste0("bday",gen)
  date<-paste0("bdate",gen)
  approx<-paste0("approxbdate",gen)
  
  dyear<-paste0("dyear",gen)
  dmonth<-paste0("dmonth",gen)
  dday<-paste0("dday",gen)
  ddate<-paste0("ddate",gen)
  
  m<-match(x[,gen],g$id)
  
  x[,sex]<-g$sex[m]
  x[,year]<-g$yearb[m]
  x[,month]<-g$monthb[m]
  x[,month]<-ifelse(x[,month]%in%c(0,"NA"),NA,x[,month])
  x[,day]<-g$dayb[m]
  x[,day]<-ifelse(x[,day]%in%c(0,"NA"),NA,x[,day])
  
  x[,dyear]<-g$yeard[m]
  x[,dmonth]<-g$monthd[m]
  x[,dmonth]<-ifelse(x[,dmonth]%in%c(0,"NA"),NA,x[,dmonth])
  x[,dday]<-g$dayd[m]
  x[,dday]<-ifelse(x[,dday]%in%c(0,"NA"),NA,x[,dday])
  
  # get a approximate bmonth
  w1<-which(!is.na(x[,year]) & is.na(x[,month]))
  if(any(w1)){
    x[w1,month]<-6
    x[w1,day]<-ifelse(is.na(x[w1,day]),30,x[w1,day])
  }
  # get approximate bday
  w2<-which(!is.na(x[,year]) & !is.na(x[,month]) & is.na(x[,day]))
  if(any(w2)){
    x[w2,day]<-15
  }
  x[,date]<-ifelse(!is.na(x[,year]),paste(x[,year],formatC(x[,month],width=2,flag=0),formatC(x[,day],width=2,flag=0),sep="-"),NA)
  x[,ddate]<-ifelse(!is.na(x[,dyear]) & !is.na(x[,dmonth]) & !is.na(x[,dday]),paste(x[,dyear],formatC(x[,dmonth],width=2,flag=0),formatC(x[,dday],width=2,flag=0),sep="-"),NA)
  # flag approximate bdate
  x[,approx]<-"no"
  x[,approx][w1]<-"month"
  x[,approx][w2]<-"day"
  x
}

### function that computes the child number considering the total number of different child
seq_inc<-function(i){
  s<-rle(i)
  rep(cumsum(s[[1]]),s[[1]])
}


### get info for both F0 and F1 generations
x<-setNames(p,c("F0","F1"))
x<-getinfo(x,gen="F0") #missing sex for 37188 individuals?
x<-x[x$sexF0%in%"f",]
x<-getinfo(x,gen="F1")
x<-x[order(x$F0,x$bdateF1,x$F1),]


### give order to each reproductive event (order) and to each child (orderc)
# bdate with missing month will be given bmonth=06 and bday=30 if bday is NA and flagged "month" for bdate
# bdate with missing day will be given bday=15 and flagged "day" for bday
# any NAs in yearb will be placed sequentially in sequential >=5 years intervals between two reproductive event, one by one, and flagged "interval"
# the ordering of NA bdate is made using the id of individuals
# should use data.table or dplyr to make this faster

getorder<-function(x,gen="F1"){
  date<-paste0("bdate",gen)
  approx<-paste0("approxbdate",gen)
  pgen<-as.integer(gsub("F","",gen))
  if(pgen<1){
    stop("Trying to get order from F0 generation")
  }
  g<-paste0("F",pgen-1)
  ddply(x,g,function(i){
    wNA<-which(is.na(i[,date]))
    if(any(wNA)){
      s<-setdiff(seq_len(nrow(i)),wNA)
      if(length(s)>1){
        d<-as.integer(diff(as.Date(i[,date]))) 
        wsup<-which(d>=(5*365))
        if(any(wsup)){
          minv<-min(length(wNA),length(wsup))
          for(j in seq_len(minv)){
            #browser()
            i[,date][wNA[j]]<-as.character(as.Date(i[,date][wsup[j]])+d[wsup[j]]/2)  
            i[,approx][wNA[j]]<-"interval"
          }
        } 
      }
    }
    i<-i[order(i[,date],i[,gen]),]
    temp<-i[,date]
    temp<-paste(temp,cumsum(is.na(temp)))
    o<-as.integer(factor(temp))
    res<-cbind(i,o,seq_inc(o),stringsAsFactors=FALSE)
    names(res)[(ncol(res)-1):ncol(res)]<-paste0(c("order","orderc"),gen)
    res
  })
}

x<-getorder(x,gen="F1")
x<-x[order(x$F0,x$bdateF1,x$F1),]


### build second dataset by obtaining all F1 that are listed as mothers in the complete mothers-childs dataset (p)
x2<-p[p$mother%in%x$F1,]
x2<-setNames(x2,c("F1","F2"))
x2<-join(x,x2,type="left")
x2<-x2[x2$sexF1%in%c("f"),]
x2<-getinfo(x2,gen="F2")
x2<-x2[order(x2$F1,x2$bdateF2,x2$F2),]

x2<-getorder(x2,gen="F2")

x2<-x2[!is.na(match(x2$F1,p$mother)),]
x2<-x2[,c("F0",names(x2)[grep("F1|F2",names(x2))])]
x2<-x2[order(x2$byearF1,x2$F1,x2$bdateF2,x2$F2),]


####################################################################################
### SIMPLE BIRTH ORDERS AND NUMBERS OF CHILD OVERLAPPING FROM x and x2, not y and y2
####################################################################################

library(data.table)

int<-2*365

x$bdateF0<-as.Date(x$bdateF0)
x$bdateF1<-as.Date(x$bdateF1)
x2$bdateF1<-as.Date(x2$bdateF1)
x2$bdateF2<-as.Date(x2$bdateF2)

x2dt<-as.data.table(x2[,c("F0","F1","F2","bdateF1","bdateF2")])
setkey(x2dt,"F0")

f<-function(i){
  xx2<-x2dt[i$F0[1]]
  if(nrow(xx2)==0){
    nb_event<-NA
    nb_offspring<-NA
  }else{
    res<-xx2$F1[which(xx2$bdateF2>=(i$bdateF1-int) & xx2$bdateF2<=(i$bdateF1+int))]
    nb_event<-length(res)
    nb_offspring<-length(unique(res))
  }
  data.frame(i,nb_event,nb_offspring,stringsAsFactors=FALSE)
}

### this takes about 4 minutes on my computer
xdf<-x
system.time(y1<-setDT(xdf)[, f(.SD), by=c("F0","F1"), .SDcols=c("F0","F1","bdateF1")])


xdt<-as.data.table(x2[,c("F0","F1","bdateF1")])
setkey(xdt,"F0")

f2<-function(i){
  xx<-xdt[i$F0[1]]
  if(nrow(xx)==0){
    nb_event<-NA
  }else{
    nb_event<-sum(which(xx$bdateF1>=(i$bdateF1-int) & xx$bdateF1<=(i$bdatebF1+int)))
  }
  data.frame(i,nb_event,nb_offspring=NA,stringsAsFactors=FALSE)
}

### this takes about 4 minutes on my computer
x2df<-x2
system.time(y2<-setDT(x2df)[, f2(.SD), by=c("F1","F2"), .SDcols=c("F0","F1","bdateF1")])


### general checking
#check F0 141967 she has a F1 with NA but it is likely her first child and not her last based on her F2 child dates
id<-y1$F0[y1$nb_event==10]
yy1<-y1[y1$F0==id[1],]
xx<-x[x$F0==id[1],c("F0","F1","bdateF1")]
xx2<-x2[x2$F0==id[1],c("F0","F1","sexF1","F2","bdateF2")]





#xx<-y[y$F0==y$F0[300],]
#xx2<-y2[y2$F1%in%yy$F1,]
#x<-join(y,y2[,c("F1","F2","OrderF1","yearbF2")],type="full")
#<-y[,c("F0","OrderF0","F1","yearbF1","nb_event","nb_offspring")]
#y<-y[order(y$F0,y$F1,y$OrderF0),]
#y2<-y2[,c("F0","F1","F2","OrderF1","yearbF2","nb_event","nb_offspring")]
#y2<-y2[order(y2$F0,y2$F1,y2$F2,y2$OrderF1,y2$yearbF2),]










