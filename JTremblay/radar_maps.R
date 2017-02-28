
library(rgdal)
library(maptools)
library(sp)
library(rgeos)
library(dagR)
library(PBSmapping)
library(ggplot2)
library(grid)
library(plyr)
library(FRutils)
### need latest version of RTools and MIkTeX to install FRutils for colo.scale function (and package devtools)
# to install FRutils type:
# devtools:::install_github("frousseu/FRutils")

#### FONCTIONS diverses ========================================================

lim<-function(x=4){
  l<-locator(n=x)
  X<-c(min(l$x),max(l$x))
  Y<-c(min(l$y),max(l$y))
  return(list(x=X,y=Y))
}


# f_dir <- function(x, type, season) {
#               if(season=="fall"){
#                    return(switch(type,
#                           max.neg = as.numeric(x)+180,
#                           zero.to.s.E = as.numeric(x)+90,
#                           max.pos = as.numeric(x),
#                           zero.to.s.W = as.numeric(x)-90
#                           ))
#                    }
#               if(season=="spring"){
#                    return(switch(type,
#                           max.neg = as.numeric(x)-180,
#                           zero.to.n.E = as.numeric(x)-90,
#                           max.pos = as.numeric(x),
#                           zero.to.n.W = as.numeric(x)+90
#                           ))
#                    }
#               }


pang<-function (dep, angl, len){
  matrix(c(dep[,1] + (sin((angl/360)*2*pi) * len), dep[,2] + (cos((angl/360)*2*pi) * len)),ncol=2)
}


LL2UTM<-function(X,Y){
  Z<-ceiling((X+180)/6)
  d<-data.frame(X,Y,stringsAsFactors=F)
  n<-unique(Z)
  l<-NULL
  x<-rep(NA,nrow(d))
  y<-rep(NA,nrow(d))
  for(i in 1:length(n)){
    w<-which(Z==n[i])
    dd<-d[w,]
    attr(dd,"projection")<-"LL"
    attr(dd,"zone")<-as.character(n[i])
    ans<-convUL(dd,km=F)
    x[w]<-ans[,1]
    y[w]<-ans[,2]
  }
  return(data.frame(xutm=x,yutm=y,zoneutm=Z,stringsAsFactors=F))
}


quad<-function(x,y){
  if(x>0 && y>0){return(90)}
  if(x>0 && y<0){return(90)}
  if(x<0 && y<0){return(270)}
  if(x<0 && y>0){return(270)}
}

### function to build rounded squares and their labels in the same direction
rsquare<-function(xy=c(1,1),az=0,ll=100,l=5,a=4,lab=10000){
  # l = longeur en km
  # a = largeur en angle degr?
  # lab = 10000 label offset in the same direction in meters
  xy<-SpatialPoints(data.frame(x=xy[1],y=xy[2]))
  proj4string(xy)<-CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  p<-pang(coordinates(xy),az,ll*1000)
  p<-SpatialPoints(data.frame(x=p[,1],y=p[,2]))
  proj4string(p)<-CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  dis<-abs(gDistance(p,xy))
  c2<-gBuffer(xy,width=dis+((l*1000)/2),quadsegs=500)
  c1<-gBuffer(xy,width=dis-((l*1000)/2),quadsegs=500)
  pp<-rbind(xy,p)
  dif<-coordinates(p)[1,]-coordinates(xy)[1,]
  ang<-quad(dif[1],dif[2])-((atan(dif[2]/dif[1])/(2*pi))*360)
  x1<-pang(coordinates(xy),ang+(a/2),dis+(l*1000))
  x2<-pang(coordinates(xy),ang-(a/2),dis+(l*1000))
  tri<-SpatialPolygons(list(Polygons(list(Polygon(rbind(x1,x2,coordinates(xy),x1))),ID=1)))
  ring<-gDifference(c2,c1)
  proj4string(tri)<-proj4string(ring)
  ans<-gIntersection(tri,ring)
  plab<-pang(coordinates(xy),az,ll*1000+lab) #determines the position of the label in the same direction
  plab<-SpatialPoints(plab)
  list(ans,plab)
}


### load data
r<-readRDS("C:/Users/rouf1703/Documents/UdeS/Consultation/JTremblay/Doc/cartesFGagnon/Reflectivity_Comp_Moy_block.rds")
r$response<-as.numeric(gsub(",",".",r$response))
r$group.ltr<-gsub(" ","",r$group.ltr)
lak <- readOGR(dsn="C:/Users/rouf1703/Documents/UdeS/Consultation/JTremblay/Doc/cartesFGagnon", layer="lak")
d<-read.csv("C:/Users/rouf1703/Documents/UdeS/Consultation/JTremblay/Doc/cartesFGagnon/Block_coordinates.csv",header=T,stringsAsFactors=F)
l<-list(x=c(375013.9,696250.7),y=c(4711139,4959871))

### location of WKR radar
WKR<-c((-1)*(79+((34+25.97/60)/60)),43+((57+50.15/60)/60))
WKR<-SpatialPoints(data.frame(x=WKR[1],y=WKR[2]))
proj4string(WKR)<-CRS("+proj=longlat +datum=WGS84")
WKR<-spTransform(WKR,CRS(proj4string(lak)))

### location of WSO radar
WSO<-c((-1)*(81+((23+3.01/60)/60)),43+((22+13.01/60)/60))
WSO<-SpatialPoints(data.frame(x=WSO[1],y=WSO[2]))
proj4string(WSO)<-CRS("+proj=longlat +datum=WGS84")
WSO<-spTransform(WSO,CRS(proj4string(lak)))

### build rounded squares and their labels
temp<-dlply(d,names(d),function(i){
  x<-rsquare(xy=coordinates(get(i$radar))[1,],az=i$azimut,ll=i$range,l=5,a=4,lab=5000)
  ans<-x[[1]]
  ans<-SpatialPolygonsDataFrame(ans,data=i,match.ID=FALSE)
  ans2<-x[[2]]
  ans2<-SpatialPointsDataFrame(ans2,data=i,match.ID=FALSE)
  list(ans,ans2)
})
rs<-do.call("rbind",lapply(temp,function(i){i[[1]]}))
rslab<-do.call("rbind",lapply(temp,function(i){i[[2]]}))

### colors to be used
colors<-rev(c("white","yellow","orange","red","darkred"))

### iterate over seasons and times to build 4 figures
it<-unique(r[,c("season","t")])
for(i in 1:nrow(it)){
  x<-r[r$season==it$season[i] & r$t==it$t[i],]
  info<-join(rs@data,x[,c(intersect(names(rs),names(x)),"response","group.ltr")],type="left")[,c("response","group.ltr")] 
  
  png(file=paste0("C:/Users/rouf1703/Documents/UdeS/Consultation/JTremblay/Doc/cartesFGagnon/Carte_Blocs_2016_",paste(it$season[i],paste0("t",it$t[i]),sep="_"),".png"),width=7,height=6,units="in",res=700,pointsize=15)
  
  par(mar=c(0,0,0,0))
  plot(lak,ylim=range(l$y),xlim=range(l$x),bg="white",border="grey20")
  plot(gBoundary(gBuffer(WKR,width=c(80)*1000,quadsegs=100)),add=T,col=gray(3/8),lty=2)
  plot(gBoundary(gBuffer(WSO,width=c(80)*1000,quadsegs=100)),add=T,col=gray(3/8),lty=2)
  
  mi<-min(info$response,na.rm=TRUE)
  ma<-max(info$response,na.rm=TRUE)
  val<-1-((info$response-mi)/(ma-mi))
  val2<-ifelse(is.na(val),0,val)
  col<-colo.scale(val2,colors)
  col[is.na(val)]<-NA
  plot(rs,add=TRUE,col=col,lwd=0.3)
  text(coordinates(rs[is.na(val),])[,1],coordinates(rs[is.na(val),])[,2],"X",adj=c(0.5,0.5),cex=0.5)
  se<-seq(mi,ma,length.out=6)
  legend("bottomright",inset=c(0.025,0.07),legend=round(se,0),pch=22,pt.bg=colo.scale(1-((se-mi)/(ma-mi)),colors),title="Density",cex=0.7,pt.lwd=0.3,pt.cex=1.5)
  text(coordinates(rslab)[,1],coordinates(rslab)[,2],info$group.ltr,cex=0.4)
  
  ### title
  text(par("usr")[2],par("usr")[4],tools:::toTitleCase(paste(it$season[i],paste0("t",it$t[i]),sep="_")),adj=c(1.5,4),font=2)
  
  ### radar centers
  points(WKR@coords[,"x"],WKR@coords[,"y"],cex=4,lwd=2)
  text(WKR@coords[,"x"],WKR@coords[,"y"],label="WKR",col="black",cex=0.6,font=1,adj=c(0.5,0.45))
  points(WSO@coords[,"x"],WSO@coords[,"y"],cex=4,lwd=2)
  text(WSO@coords[,"x"],WSO@coords[,"y"],label="WSO",col="black",cex=0.6,font=1,adj=c(0.5,0.45))
  
  ### scale and north arrow
  sc<-list(x=371690.1,y=4907462)
  text(((l$x[2]-l$x[1])*0.05)+l$x[1],((l$y[2]-l$y[1])*0.905)+l$y[1],label="N",font=2,cex=2)
  arrows(x0=((l$x[2]-l$x[1])*0.05)+l$x[1],y0=((l$y[2]-l$y[1])*0.84)+l$y[1],x1=((l$x[2]-l$x[1])*0.05)+l$x[1],y1=((l$y[2]-l$y[1])*0.98)+l$y[1],lwd=2,length=0.1) 
  rect(xleft=sc$x, ybottom=sc$y-4000, xright=sc$x+20000, ytop=sc$y,col="white",border="black")
  rect(xleft=sc$x+20000, ybottom=sc$y-4000, xright=sc$x+40000, ytop=sc$y,col="black",border="black")
  text(sc$x,sc$y+7000,label="0",cex=0.5)
  text(sc$x+40000,sc$y+7000,label="40 km",adj=c(0.25,0.5),cex=0.5)
  
  ### locations
  #loc<-locator(5) #use locator function to sequentially determine position of labels
  loc<-list(x=c(389319,555584,694870,584445,629619),y=c(4877793,4947436,4818188,4726586,4835129))
  lab<-c("Lake\nHuron","Georgian Bay","Lake\nOntario","Lake\nErie","Toronto")
  for(j in seq_along(loc$x)){
    if(lab[j]=="Toronto"){
      points(loc$x[j],loc$y[j],pch=16,cex=0.3)
      text(loc$x[j],loc$y[j],lab[j],adj=c(-0.1,0.5),cex=0.5)
    }else{
      text(loc$x[j],loc$y[j],lab[j],adj=c(0.5,0.5),cex=0.5)
    }
  }
  
  dev.off()
  
}

