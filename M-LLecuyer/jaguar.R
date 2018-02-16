

library(splancs)
library(geoR)
library(gstat)
library(sp)
library(scales)
library(ncf)
library(car)
library(data.table)
library(readxl)
library(raster)
library(visreg)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(INLA)
library(velox)
library(geostatsp)
library(FRutils)
library(RColorBrewer)
library(tmap)
library(devtools)
library(doParallel)
library(foreach)
library(visreg)
library(plyr)
library(rasterVis)
library(INLAutils)
library(randomForest)
library(glmnet)
library(colorRamps)
library(fastshp)
library(sf)


####################################################################
### function to build waic table from a list of geostatsp models
####################################################################
waictab<-function(x){
  stopifnot(!is.null(names(x)))
  waic<-sapply(x,function(i){
    i$inla$waic$waic  
  })
  p.eff<-sapply(x,function(i){
    i$inla$waic$p.eff
  })
  dwaic<-waic-min(waic)
  w<-sapply(dwaic,function(i){
    exp((-1/2)*i)/sum(exp((-1/2)*dwaic))
  })
  d<-data.frame(model=names(x),p.eff=round(p.eff,2),waic=round(waic,2),dwaic=round(dwaic,2),w=round(w,2))
  #d<-d[order(d$dwaic),]
  d
}

################################################
### load data
################################################
setwd("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc")

#d<-as.data.frame(read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/LandEco_Complet_22_05.xlsx"),stringsAsFactors=FALSE)

d=as.data.frame(fread('LandEco_Complet_12_12_temp.txt', dec=','))
d[]<-lapply(d,function(i){
  if(any(grep(",",i))){
    as.numeric(gsub(",",".",i))
  }else{
    i
  }
})

source("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/GitHub/vif.R")

vif2<-function(formula,data){
  x<-model.matrix(formula,data=data)
  m<-match("(Intercept)",dimnames(x)[[2]])
  if(!is.na(m)){
    x<-x[,-m]  
  }
  vif(x)
}


### build shapefile of locations
ll<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
prj<-"+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
ds<-d

coordinates(ds)<- ~X+Y
proj4string(ds) <- ll
ds<-spTransform(ds,CRS(prj))

d$X<-coordinates(ds)[,1]
d$Y<-coordinates(ds)[,2]

d$x<-d$X
d$y<-d$Y

# ds is a spatial object that can be plotted

d$Cat_Typ<-as.factor(d$Cat_Typ)
plot(ds)


#########################################
### import raster and build pred grid
#########################################

### this section may not be relevant and is not used, but may be usefull to compare with data .txt

### raster landcover####

path<-"Landcover_2015_extended.tif"
code<-read.table("Legend.txt",header=TRUE,sep=",",stringsAsFactors=FALSE)
code=fread('Legend.txt',header=TRUE,sep=",",stringsAsFactors=FALSE)
code<-code[code$Code!=0,]
code$Category[which(code$Category=="Agriculture\t")]<-"Agriculture"
r <- raster(path)
rpred<-raster(ext=extent(r),ncols=25,nrows=25)
proj4string(rpred)<-proj4string(r)

### vérifier si cela donne les mêmes résultats si on prend juste le data.frame et non les calculs
### get proportions at observations
w<-1000
p<-gBuffer(SpatialPoints(coordinates(ds)),width=w,byid=TRUE)
proj4string(p)<-proj4string(r)
v<-velox(r)
eobs<-v$extract(p,fun=function(x){paste(x,collapse="_")})
ans<-lapply(eobs,function(i){
  x<-table(unlist(strsplit(i,"_")))
  x<-x[names(x)!=0] # ignore background values in raster
  s<-setdiff(code$Code,names(x))
  if(length(s)==0){
    x<-x[order(names(x))]
  }else{
    temp<-rep(0,length(s))
    names(temp)<-s
    x<-c(x,temp)
    x<-x[order(names(x))]
  }
  100*x/sum(x)
})
obs<-as.data.frame(do.call("rbind",ans))
names(obs)<-gsub(" ","_",code$Category[match(names(obs),code$Code)])
d<-cbind(d,obs,stringsAsFactors=FALSE)

ds@data<-cbind(ds@data,obs)

# g is a SpatialPixelsDataFrame

# indFor = Bajos
# SecFor = Secondary
# matFor= Selva baja + selva mediana + selva alta + Subcadifolia
# Past= Pasture
# Agr = Agriculture + milpa



###################################################################################
### add pop index#
###################################################################################

pop<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="poblacion")
pop<-read.shp("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/poblacion")

### for the obs

di<-gDistance(pop,ds,byid=TRUE)
f<-sqrt
ds$index_pop2<-rowSums(t(apply(di,1,function(i){exp(((log(0.1)/2000))*i)*f(pop$POBTOT)})))
ds$index_pop4<-rowSums(t(apply(di,1,function(i){exp(((log(0.1)/4000))*i)*f(pop$POBTOT)})))
ds$index_pop8<-rowSums(t(apply(di,1,function(i){exp(((log(0.1)/8000))*i)*f(pop$POBTOT)})))

d$index_pop2<-ds$index_pop2
d$index_pop4<-ds$index_pop4
d$index_pop8<-ds$index_pop8


### for the raster

p<-xyFromCell(rpred,1:ncell(rpred),spatial=TRUE)

di<-gDistance(pop,p,byid=TRUE)
f<-sqrt
rpop<-rpred
rpop1<-setValues(rpop,rowSums(t(apply(di,1,function(i){exp(((log(0.1)/2000))*i)*f(pop$POBTOT)}))))
rpop2<-setValues(rpop,rowSums(t(apply(di,1,function(i){exp(((log(0.1)/4000))*i)*f(pop$POBTOT)}))))
rpop3<-setValues(rpop,rowSums(t(apply(di,1,function(i){exp(((log(0.1)/8000))*i)*f(pop$POBTOT)}))))
                 
rpop<-stack(rpop1,rpop2,rpop3)
names(rpop)<-c("index_pop2","index_pop4","index_pop8")


###################################################################################
### add dist to road
###################################################################################

roads<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="carreteras")

### for the obs

di<-gDistance(ds,gLineMerge(roads),byid=TRUE)[1,]

d$Dist_Road<-di
ds$Dist_Road<-di

### import raster poblacion y carreteras###

# lire les shp
x1<-readOGR(".",layer="poblacion")
x2<-readOGR(".",layer="carreteras")

# lire le raster
rbis<-"Landcover_2015_extended.tif"
rbis <- stack(rbis)

# créer un raster de prédictions et le projeter
ras<-rpred
proj4string(ras)<-proj4string(rbis)

# faire deux rasters réponses
ras1<-ras
ras2<-ras

# extraire le centroide des pixels
coo<-SpatialPoints(coordinates(ras),proj4string=CRS(proj4string(ras)))

# calculer toutes les paires de distances
dis1<-gDistance(coo,x1,byid=TRUE)
dis2<-gDistance(coo,x2,byid=TRUE)

# trouver la valeur minimale pour chaque centroid
ans1<-apply(dis1,2,min)
ans2<-apply(dis2,2,min)

# écrire les valeurs dans les raster
ras1[]<-ans1
ras2[]<-ans2

# stacker les raster
rroad<-stack(ras1,ras2)

names(rroad)<-c("Dist_Pop","Dist_Road")

# faire une petite fonction pour visualiser les résultats
foo<-function(i){
  plot(x1,add=TRUE,pch=1)
  plot(x2,add=TRUE,lwd=2)
}

# visualiser les résultats
#plot(rroad,addfun=foo)

################################################################
### test
#summary(inla(Attack~Cat_Typ+Bridge_tg+Branch_tg,data=d,family="binomial",control.compute=list(waic=TRUE)))


#ds2<-ds[ds$Cat_Typ!="S",]
#ds2$Cat_Typ<-as.factor(as.character(ds2$Cat_Typ))

#Attack~P_Past_g
#Attack~Cat_Typ+P_Agr_g+P_Past_g+index_pop8+Dist_Road+Dens_Liv+Liv_Mgmt

#fit<-glgm(Attack~Cat_Typ+P_Past_g+index_pop8+Dist_Road,
#          data=ds,
#          grid=20,
#          covariates=NULL, 
#          family="binomial", 
#          buffer=10000,
#          shape=1,
#          priorCI=list(sd=c(0.5,3),range=c(5000,20000)),
#          control.compute=list(waic=TRUE,dic=TRUE,mlik=TRUE),#,config=TRUE),
#          #control.predictor=list(compute=TRUE,link=1)
#)

#summary(fit$inla)

#inla.rerun(fit$inla)$waic

#m<-fit

#par(mfrow=c(1,2))

#ma<-max(c(m$parameters$sd$prior[,2],m$parameters$sd$posterior[,2]))
#plot(m$parameters$sd$prior,type="l",xlab='standard deviation', ylab='density',lty=2,xlim=c(0,5),ylim=c(0,ma))
#lines(m$parameters$sd$posterior,lty=1)
#legend("topright", lty=2:1, legend=c("prior","posterior"))

#ma<-max(c(m$parameters$range$prior[,2],m$parameters$range$posterior[,2]))
#plot(m$parameters$range$posterior,type="l",xlim = c(0,50*1000),xlab='range (m)', ylab='density',lty=1,ylim=c(0,ma))
#lines(m$parameters$range$prior,lty=2)
#legend("topright", lty=2:1, legend=c("prior","posterior"))


#################selection variable VIF###########################@@
##Colinearity in our data using VIF#####

source('vif.R')
summary(d)
head(d)

#VIF at the 0,5 km scale

VifLEp=model.matrix(~d$Liv_Mgmt+d$WatPA+d$Cat_Typ+d$Dens_Liv+d$P_matFor_p+d$P_indFor_p+d$P_SecFor_p+d$P_Past_p+d$P_Agr_g+d$Nombre_MinPS_p+d$Den__MinPS_p+d$Perforation_p+d$Bridge_p+d$Branch_p+d$Human_dens_p+d$Dist_HabP+d$Dist_Road)


VifLEp2=VifLEp[,-1] ##enlever l'intercept

vif(VifLEp2)

#drop sec forest

VifLEp3=VifLEp2[,-11]

VifLEp2=VifLEp[,-1] 

vif(VifLEp3)

##we take away the Nombre_MinPS_S

VifLEp4=VifLEp3[,-12]

vif(VifLEp4)

## VIF analysis at the 2,5km scale
VifLEm=model.matrix(~d$Liv_Mgmt+d$WatPA+d$Cat_Typ+d$Dens_Liv+d$P_matFor_m+d$P_indFor_m+d$P_SecFor_m+d$P_Past_m+d$P_Agr_m+d$Nombre_MinPS_m+d$Den__MinPS_m+d$Perforation_m+d$Bridge_m+d$Branch_m+d$Human_dens_m+d$Dist_HabP+d$Dist_Road)


VifLEm2=VifLEm[,-1] 

vif(VifLEm2)

#We do not consider the score of pasture in this part of this sutdy and just decided to drop the variable P_SecFor_p

VifLEm3=VifLEm2[,-11]

vif(VifLEm3)

#When sec for is drop as well, VIF score of mat forest go under 3

VifLEm4=VifLEm3[,-10]

vif(VifLEm4)

##We drop nombre of min patch

VifLEm5=VifLEm4[,-11]

vif(VifLEm5)

## VIF analysis at the 5km scale

VifLEg=model.matrix(~d$Liv_Mgmt+d$WatPA+d$Cat_Typ+d$Dens_Liv+d$P_matFor_g+d$P_indFor_g+d$P_SecFor_g+d$P_Past_g+d$P_Agr_g+d$Nombre_MinPS_g+d$Den__MinPS_g+d$Perforation_g+d$Bridge_g+d$Branch_g+d$Human_dens_g+d$Dist_HabP+d$Dist_Road)


VifLEg2=VifLEg[,-1] 

vif(VifLEg2)

#We drop the P_SecFor_g variable

VifLEg3=VifLEg2[,-10]

vif(VifLEg3)

#When pasture is drop as well, 

VifLEg4=VifLEg3[,-10]
vif(VifLEg4)

##We drop nombre of min patch

VifLEg5=VifLEg4[,-11]

vif(VifLEg5)

##We drop branch of min patch

VifLEg6=VifLEg5[,-14]

vif(VifLEg6)

## VIF analysis at the 10km scale

VifLEtg=model.matrix(~d$Liv_Mgmt+d$WatPA+d$Cat_Typ+d$Dens_Liv+d$P_matFor_tg+d$P_indFor_tg+d$P_SecFor_tg+d$P_Past_tg+d$P_Agr_tg+d$Nombre_MinPS_tg+d$Den__MinPS_tg+d$Perforation_tg+d$Bridge_tg+d$Branch_tg+d$Human_dens_tg+d$Dist_HabP+d$Dist_Road)

VifLEtg2=VifLEtg[,-1] 

vif(VifLEtg2)

#We drop the P_SecFor_g variable

VifLEtg3=VifLEtg2[,-10]

vif(VifLEtg3)

#When pasture is drop as well, VIF score of mat forest is still higher then 3 but considering the other model, we decide to keep the variables left

VifLEtg4=VifLEtg3[,-10]

vif(VifLEtg4)

##We drop nombre of min patch

VifLEtg5=VifLEtg4[,-11]

vif(VifLEtg5)

##We drop branch of min patch

VifLEtg6=VifLEtg5[,-14]

vif(VifLEtg6)

######################################################################
### build model formulas
######################################################################


####################@
###################@@
##################@

#Natural habitat hypothesis## BISSSSSS on test avec le fragmentation fans F###

##(SC1) The amount of forest favors the occurrence of attacks:%matFor, %indFor, %SecFor##

##Colinearity in our data using VIF#####

source('vif.R')
summary(d)
head(d)

#VIF SC1bisp

VifLESC1p=model.matrix(~d$Cat_Typ+d$P_matFor_p+d$P_indFor_p+d$P_SecFor_p)

VifLESC1p2=VifLESC1p[,-1] ##enlever l'intercept

vif(VifLESC1p2)

#VIF SC1bism

VifLESC1m=model.matrix(~d$Cat_Typ+d$P_matFor_m+d$P_indFor_m+d$P_SecFor_m)

VifLESC1m2=VifLESC1m[,-1] ##enlever l'intercept

vif(VifLESC1m2)

#We drop the P_SecFor_m variable

VifLESC1m3=VifLESC1m2[,-3]

vif(VifLESC1m3)

#VIF SC1bisg

VifLESC1g=model.matrix(~d$Cat_Typ+d$P_matFor_g+d$P_indFor_g+d$P_SecFor_g)

VifLESC1g2=VifLESC1g[,-1] ##enlever l'intercept

vif(VifLESC1g2)

#We drop the P_SecFor_tg variable

VifLESC1g3=VifLESC1g2[,-3]

vif(VifLESC1g3)

#VIF SC1bistg

VifLESC1tg=model.matrix(~d$Cat_Typ+d$P_matFor_tg+d$P_indFor_tg+d$P_SecFor_tg)

VifLESC1tg2=VifLESC1tg[,-1] ##enlever l'intercept

vif(VifLESC1tg2)

#We drop the P_SecFor_tg variable

VifLESC1tg3=VifLESC1tg2[,-3]

vif(VifLESC1tg3)

#############################################
### MODEL LIST
###############################################


ml<-list()



# Model determination and model selection 

###Forêt actuelle

SC1bis0<-Attack~Cat_Typ
SC1bis1<-Attack~Cat_Typ+P_matFor_p+P_indFor_p
SC1bis2<-Attack~Cat_Typ+P_matFor_m+P_indFor_m
SC1bis3<-Attack~Cat_Typ+P_matFor_g+P_indFor_g
SC1bis4<-Attack~Cat_Typ+P_matFor_tg+P_indFor_tg

SC1bis<-list(SC1bis0=SC1bis0,SC1bis1=SC1bis1,SC1bis2=SC1bis2,SC1bis3=SC1bis3,SC1bis4=SC1bis4)
ml[[length(ml)+1]]<-SC1bis
sapply(SC1bis,vif2,data=d)


## echelle temporelle 

TempSC1bis0<-Attack~Cat_Typ
TempSC1bis1<-Attack~Cat_Typ+P_matFor_p_2000+P_indFor_p_2000
TempSC1bis2<-Attack~Cat_Typ+P_matFor_m_2000+P_indFor_m_2000
TempSC1bis3<-Attack~Cat_Typ+P_matFor_g_2000+P_indFor_g_2000
TempSC1bis4<-Attack~Cat_Typ+P_matFor_tg_2000+P_indFor_tg_2000

TempSC1bis<-list(TempSC1bis0=TempSC1bis0,TempSC1bis1=TempSC1bis1,TempSC1bis2=TempSC1bis2,TempSC1bis3=TempSC1bis3,TempSC1bis4=TempSC1bis4)
ml[[length(ml)+1]]<-TempSC1bis
sapply(TempSC1bis,vif2,data=d)


####Temp
Temp0<-Attack~Cat_Typ
Temp1<-Attack~Cat_Typ+P_matFor_g+P_indFor_g
Temp2<-Attack~Cat_Typ+P_matFor_g_2000+P_indFor_g_2000
Temp3<-Attack~Cat_Typ+P_matFor_g_2000+P_indFor_g_2000+P_matFor_g+P_indFor_g

Temp<-list(Temp0=Temp0,Temp1=Temp1,Temp2=Temp2,Temp3=Temp3)
ml[[length(ml)+1]]<-Temp
sapply(Temp,vif2,data=d)

####(SC2) Perforation dans le territoire.

V2SC2bis0<-Attack~Cat_Typ
V2SC2bis1<-Attack~Cat_Typ+Perforation_p
V2SC2bis2<-Attack~Cat_Typ+Perforation_m
V2SC2bis3<-Attack~Cat_Typ+Perforation_g
V2SC2bis4<-Attack~Cat_Typ+Perforation_tg

V2SC2bis<-list(V2SC2bis0=V2SC2bis0,V2SC2bis1=V2SC2bis1,V2SC2bis2=V2SC2bis2,V2SC2bis3=V2SC2bis3,V2SC2bis4=V2SC2bis4)

ml[[length(ml)+1]]<-V2SC2bis
sapply(V2SC2bis,vif2,data=d)

### Selection de l'échelle tg mais écart moyen avec les autres ! ###


###(SC3) Presence or absence of water (WatPA) in the pasture area favors attack on livestock##,
##(SC4) combined effect of the structural characteristics of the natural habitat are important to explain occurrence of attacks##

##Structural characteristics (SC): jaguar attack is directly conditioned by the presence of natural habitat nearby; 
#previous studies have shown that the amount of forest and the presence of water have influenced the occurrence of attacks.##

#VIF pour SC

VifLESC=model.matrix(~d$Cat_Typ+d$WatPA+d$Cat_Typ+d$P_matFor_tg+d$P_indFor_tg+d$Perforation_tg)

VifLESC2=VifLESC[,-1] ##enlever l'intercept

vif(VifLESC2)

# Model determination and model selection 

V2SC0<-Attack~Cat_Typ
V2SC1<-Attack~Cat_Typ+P_matFor_g+P_indFor_g+P_matFor_g_2000+P_indFor_g_2000
V2SC2<-Attack~Cat_Typ+Perforation_tg
V2SC3<-Attack~Cat_Typ+WatPA
V2SC4<-Attack~Cat_Typ+P_matFor_g+P_indFor_g+WatPA+Perforation_tg+P_matFor_g_2000+P_indFor_g_2000

V2SC<-list(V2SC0=V2SC0,V2SC1=V2SC1,V2SC2=V2SC2, V2SC3=V2SC3,V2SC4=V2SC4)
ml[[length(ml)+1]]<-V2SC
sapply(V2SC,vif2,data=d)

##Selection de V2SC4 meilleure que les autres###

##Functional characteristics (FC): Landscape fragmentation can also influence jaguar movement and jaguar presence and in consequences, 
##potentially influence the occurrence of attacks##

###(FC4) Local aggregation of forest influence the presence of prey and can influence the occurrence of attacks##

#VIF SC2p

VifLESC2p=model.matrix(~d$Bridge_p+d$Branch_p)

VifLESC2p2=VifLESC2p[,-1] ##enlever l'intercept

vif(VifLESC2p2)

#VIF SC2m

VifLESC2m=model.matrix(~d$Bridge_m+d$Branch_m)

VifLESC2m2=VifLESC2m[,-1] ##enlever l'intercept

vif(VifLESC2m2)

#VIF SC2g

VifLESC2g=model.matrix(~d$Bridge_g+d$Branch_g)

VifLESC2g2=VifLESC2g[,-1] ##enlever l'intercept

vif(VifLESC2g2)

#VIF SC2tg

VifLESC2tg=model.matrix(~d$Bridge_tg+d$Branch_tg)

VifLESC2tg2=VifLESC2tg[,-1] ##enlever l'intercept

vif(VifLESC2tg2)


# Model determination and model selection 

SC2bis0<-Attack~Cat_Typ
SC2bis1<-Attack~Cat_Typ+Bridge_p+Branch_p
SC2bis2<-Attack~Cat_Typ+Bridge_m+Branch_m
SC2bis3<-Attack~Cat_Typ+Bridge_g+Branch_g
SC2bis4<-Attack~Cat_Typ+Bridge_tg+Branch_tg

SC2bis<-list(SC2bis0=SC2bis0,SC2bis1=SC2bis1,SC2bis2=SC2bis2,SC2bis3=SC2bis3,SC2bis4=SC2bis4)
ml[[length(ml)+1]]<-SC2bis
sapply(SC2bis,vif2,data=d)

## Selection de l'échelle tg avec assez bonne différence avec le modèle suivant###

##(FC1) The presence of minimum patch size (2km2) influence jaguar movement

##(FC1) The presence of minimum patch size (2km2) (Dens-MinP) influence jaguar movement##

FC1bis0<-Attack~Cat_Typ
FC1bis1<-Attack~Cat_Typ+Den__MinPS_p
FC1bis2<-Attack~Cat_Typ+Den__MinPS_m
FC1bis3<-Attack~Cat_Typ+Den__MinPS_g
FC1bis4<-Attack~Cat_Typ+Den__MinPS_tg

FC1bis<-list(FC1bis0=FC1bis0,FC1bis1=FC1bis1,FC1bis2=FC1bis2,FC1bis3=FC1bis3,FC1bis4=FC1bis4)
ml[[length(ml)+1]]<-FC1bis
sapply(FC1bis,vif2,data=d)

## Selection de Den__MinPS_g mais la différence n'est pas grande 

# Model determination  and model selection FC

#VIF FC

VifLEFC=model.matrix(~d$Cat_Typ+d$Dist_HabP+d$Den__MinPS_g+d$Bridge_tg+d$Branch_tg)

VifLEFC2=VifLEFC[,-1] ##enlever l'intercept

vif(VifLEFC2)

# Model determination and model selection 

V2FC0<-Attack~Cat_Typ
V2FC1<-Attack~Cat_Typ+Den__MinPS_g
V2FC2<-Attack~Cat_Typ+Dist_HabP
V2FC3<-Attack~Cat_Typ+Den__MinPS_g+Dist_HabP+Bridge_tg+Branch_tg
V2FC4<-Attack~Cat_Typ+Bridge_tg+Branch_tg

V2FC<-list(V2FC0=V2FC0,V2FC1=V2FC1,V2FC2=V2FC2,V2FC3=V2FC3,V2FC4=V2FC4)
ml[[length(ml)+1]]<-V2FC
sapply(V2FC,vif2,data=d)


### Selection de V2FC3 

# Model determination  and model selection for the natural habitat hypothesis 

#VIF for N

V2VifLEN=model.matrix(~d$Cat_Typ+d$P_matFor_g+d$P_indFor_g+d$WatPA+d$Perforation_p+d$Bridge_p+d$Branch_p+d$P_SecFor_p)

V2VifLEN2=V2VifLEN[,-1] ##enlever l'intercept

vif(V2VifLEN2)

##Model determination and model selection 

V2N0<-Attack~Cat_Typ
V2N1<-V2SC1
V2N2<-V2FC3
V2N3<-Attack~Cat_Typ+Den__MinPS_g+Dist_HabP+Bridge_g+Branch_g+P_matFor_tg+P_indFor_g+WatPA+Perforation_tg+P_matFor_g_2000+P_indFor_g_2000 ## + SC + FC

V2N<-list(V2N0=V2N0,V2N1=V2N1,V2N2=V2N2,V2N3=V2N3)
ml[[length(ml)+1]]<-V2N
sapply(V2N,vif2,data=d)


#Selection de N1 m^ême si N3 proche mais bien meilleure que 2

###Human population and activity (HPA): human activity has been shown to influence jaguar distribution
##and the occurrence of attacks while human settlement has been shown to have a negative effect on jaguar 
##presence and the occurrence of attacks

##(HPA 1) The amount of land dedicated to agricultural (%Agr) activity influence the presence and movement of jaguar and decrease the occurrence of attacks

HPA1bis0<-Attack~Cat_Typ
HPA1bis1<-Attack~Cat_Typ+P_Agr_p
HPA1bis2<-Attack~Cat_Typ+P_Agr_m
HPA1bis3<-Attack~Cat_Typ+P_Agr_g
HPA1bis4<-Attack~Cat_Typ+P_Agr_tg

HPA1bis<-list(HPA1bis0=HPA1bis0,HPA1bis1=HPA1bis1,HPA1bis2=HPA1bis2,HPA1bis3=HPA1bis3,HPA1bis4=HPA1bis4)
ml[[length(ml)+1]]<-HPA1bis
sapply(HPA1bis,vif2,data=d)


## Selection de l'échelle g écart moyen avec les autres.

##(HPA2) Jaguar attack will be favored by the presence of a minimum amount of pasture (%Past) but will also decrease in an area dominated by pasture, (

HPA2bis0<-Attack~Cat_Typ
HPA2bis1<-Attack~Cat_Typ+P_Past_p
HPA2bis2<-Attack~Cat_Typ+P_Past_m
HPA2bis3<-Attack~Cat_Typ+P_Past_g
HPA2bis4<-Attack~Cat_Typ+P_Past_tg
HPA2bis<-list(HPA2bis0=HPA2bis0,HPA2bis1=HPA2bis1,HPA2bis2=HPA2bis2,HPA2bis3=HPA2bis3,HPA2bis4=HPA2bis4)
ml[[length(ml)+1]]<-HPA2bis
sapply(HPA2bis,vif2,data=d)


# selection de l'echelle p mais ecar moyen avec echelle tg 

##(HPA3) Human presence index pop

HPA3bis0<-Attack~Cat_Typ
HPA3bis1<-Attack~Cat_Typ+index_pop2
HPA3bis2<-Attack~Cat_Typ+index_pop4
HPA3bis3<-Attack~Cat_Typ+index_pop8

HPA3bis<-list(HPA3bis0=HPA3bis0,HPA3bis1=HPA3bis1,HPA3bis2=HPA3bis2,HPA3bis3=HPA3bis3)
ml[[length(ml)+1]]<-HPA3bis
sapply(HPA3bis,vif2,data=d)


# selection de l'echelle 8

# Model determination and model selection for the human effect, HPA

##VIF de HPA

VifLEHPA=model.matrix(~d$Cat_Typ+d$P_Agr_g+d$P_Past_tg+d$index_pop8)

VifLEHPA2=VifLEHPA[,-1] ##enlever l'intercept

vif(VifLEHPA2)

#Selection

HPA0<-Attack~Cat_Typ
HPA1<-Attack~Cat_Typ+P_Agr_g
HPA2<-Attack~Cat_Typ+P_Past_g
HPA3<-Attack~Cat_Typ+index_pop8+Dist_Road
HPA4<-Attack~Cat_Typ+P_Agr_g+P_Past_g+index_pop8+Dist_Road

HPA<-list(HPA0=HPA0,HPA1=HPA1,HPA2=HPA2,HPA3=HPA3,HPA4=HPA4)
ml[[length(ml)+1]]<-HPA
sapply(HPA,vif2,data=d)


#Selection de HP4 pas mal au dessus des autres


##Livestock production (LP): Attack of jaguar depend on the type of cattle and the management practice used by ranchers

##VIF de LP

VifLELP=model.matrix(~d$Cat_Typ+d$Dens_Liv+d$Liv_Mgmt)

VifLELP2=VifLELP[,-1] ##enlever l'intercept

vif(VifLELP2)

# Model determination and model selection for the human effect, livestock production hypothesis LP
LP0<-Attack~Cat_Typ
LP1<-Attack~Cat_Typ+Dens_Liv
LP2<-Attack~Cat_Typ+Liv_Mgmt
LP3<-Attack~Cat_Typ+Dens_Liv+Liv_Mgmt

LP<-list(LP0=LP0,LP1=LP1,LP2=LP2,LP3=LP3)
ml[[length(ml)+1]]<-LP
sapply(LP,vif2,data=d)


## Selection de LP3 bien meilleure 

# Model determination and model selection for the human effect, H

##VIF de H

VifLEH=model.matrix(~d$Cat_Typ+d$P_Agr_g+d$P_Past_tg+d$Human_dens_g+d$Dist_Road+d$Dens_Liv+d$Liv_Mgmt)

VifLEH2=VifLEH[,-1] ##enlever l'intercept

vif(VifLEH2)

##Model determination

H0<-Attack~Cat_Typ
H1<-Attack~Cat_Typ+P_Agr_g+P_Past_g+index_pop8+Dist_Road
H2<-Attack~Cat_Typ+Dens_Liv+Liv_Mgmt
H3<-Attack~Cat_Typ+P_Agr_g+P_Past_g+index_pop8+Dist_Road+Dens_Liv+Liv_Mgmt

H<-list(H0=H0,H1=H1,H2=H2,H3=H3)
ml[[length(ml)+1]]<-H
sapply(H,vif2,data=d)


##Selection de H1 meilleure que les autre mais ecart moyen

# Final model selection between nature and human effect model 


## VIF for general model

VifLEF=model.matrix(~d$Cat_Typ+d$Perforation_p+d$Bridge_p+d$Branch_p+d$P_Agr_g+d$P_Past_tg+d$Human_dens_g+d$Dist_Road)

VifLEF2=VifLEF[,-1] ##enlever l'intercept

vif(VifLEF2)

#Model determination
F0<-Attack~Cat_Typ
F1<-Attack~Cat_Typ+Den__MinPS_g+Dist_HabP+Bridge_p+Branch_p
F2<-Attack~Cat_Typ+Cat_Typ+P_Agr_g+P_Past_g+index_pop8+Dist_Road+Dens_Liv+Liv_Mgmt
F3<-Attack~Cat_Typ+Cat_Typ+Den__MinPS_g+Dist_HabP+Bridge_tg+Branch_tg+P_Agr_g+P_Past_g+index_pop8+Dist_Road+Dens_Liv+Liv_Mgmt

F<-list(F0=F0,F1=F1,F2=F2,F3=F3)
ml[[length(ml)+1]]<-F
sapply(F,vif2,data=d)


#############################
### run models
#############################

########################
### exploration techniques

### random Forest
vars<-unique(unlist(lapply(ml,function(i){lapply(i,function(j){all.vars(i[[3]])})})))
rf<-randomForest(as.factor(Attack)~.,data=d[,c("X","Y",vars)],ntree=10000)
varImpPlot(rf)

### LASSo
X<-as.matrix(d[,setdiff(vars,c("Attack","Cat_Typ"))])
X<-cbind(X,Cat_TypS=ifelse(d$Cat_Typ=="S",1,0))
X<-cbind(X,Cat_TypB=ifelse(d$Cat_Typ=="B",1,0))
X<-cbind(X,Cat_TypC=ifelse(d$Cat_Typ=="C",1,0))
X<-apply(X,2,scale)
#corrplot:::corrplot(cor(X), order = "hclust")
Y=d$Attack
lasso<-glmnet(X,Y,family="binomial",alpha=1,maxit=10^6)
lassocv<-cv.glmnet(X,Y,family="binomial",alpha=1,maxit=10^6)
coef(lassocv, s = "lambda.min",alpha=1)


####################################
### glgm
registerDoParallel(6) 
getDoParWorkers()

#ds$Attack[nrow(ds)]<-NA
#ds$Attack[nrow(ds)]<-NA
#ds<-ds[rep(1:nrow(ds),each=2),]
#ds$Attack[(nrow(ds)/2):nrow(ds)]<-NA

### Ce code était pour recompiler le package en hackant la fonction glgm pour qu'elle n'enlève pas les lignes avec les NA dans la variable réponse. Avec la version 1.6.0, ceci est corrigé comparativement à la version d,avant, 1.5.4

#untar("C:/Users/rouf1703/Downloads/geostatsp_1.6.0.tar.gz",exdir="C:/Users/rouf1703/Downloads",list=TRUE)  ## check contents
#untar("C:/Users/rouf1703/Downloads/geostatsp_1.6.0.tar.gz",exdir="C:/Users/rouf1703/Downloads")

#lf<-list.files("C:/Users/rouf1703/Downloads/geostatsp/R",full.names=TRUE)
#sapply(lf,source)

#build("C:/Users/rouf1703/Downloads/geostatsp","C:/Users/rouf1703/Downloads")
#install.packages("C:/Users/rouf1703/Downloads/geostatsp_1.6.0.tar.gz",repos=NULL)


cvar<-combn(vars[-(1:2)],3,simplify=FALSE)
ml2<-lapply(cvar,function(i){
  as.formula(paste0("Attack~Cat_Typ+",paste(i,collapse="+")),env = .GlobalEnv) 
})



m<-foreach(i=1:length(ml),.packages=c("raster","sp","geostatsp")) %dopar% {
#m<-for(i in 1:length(ml)){
     lapply(ml[[i]],function(j){
       glgm(j, 
         data=ds,
         grid=20,
         covariates=NULL, 
         family="binomial", 
         buffer=10000,
         shape=1,
         priorCI=list(sd=c(0.4,4),range=c(2000,50000)),
         control.compute=list(waic=TRUE,dic=TRUE,mlik=TRUE),
         num.threads=1 # to try and get more stable results, since foreach is already parallel
       )
     }) 
}

aic<-lapply(m,waictab)
#aic<-lapply(m,function(i){i$inla$waic$waic})

aic<-lapply(seq_along(aic),function(i){
  var<-gsub("Attack ~| ","",as.character(unlist(ml[[i]])))
  res<-cbind(aic[[i]],mod=names(ml[[i]]),var)
  res<-res[order(res$dwaic),]
  res
})

lapply(ml,function(i){lapply(i,vif2,data=d)})

summary(m[[16]][[2]]$inla)
plot(d$P_indFor_g_2000,d$P_indFor_g)


######################################################################
### prediction rasters stacks
######################################################################

### build grid to extract pixel values around a buffer where predictions are going to be made
### 
### check what to do with background values and make sure
### 
### check if proportions sum to 1

path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Landcover_2015_extended.tif"
r <- raster(path)
code<-read.table("Legend.txt",header=TRUE,sep=",",stringsAsFactors=FALSE)
code$Category[which(code$Category=="Agriculture\t")]<-"Agriculture"

code$Category[which(code$Category%in%c("Agriculture","Milpa"))]<-"P_Agr"
code$Category[which(code$Category%in%c("Pasture"))]<-"P_Past"
code$Category[which(code$Category%in%c("Selva baja","Selva mediana","Selva alta mediana","Subcaducifolia"))]<-"P_matFor"
code$Category[which(code$Category%in%c("Secondary"))]<-"P_secFor"
code$Category[which(code$Category%in%c("Bajos"))]<-"P_indFor"

l<-sapply(split(code,code[,"Category"]),function(i){i[,"Code"]})
mn<-unlist(l,use.names=FALSE)
lab<-rep(names(l),times=sapply(l,length))
labi<-as.integer(factor(lab))
mat<-cbind(mn,labi)
r <- reclassify(r, mat)

ras <- rpred
ras[] <- runif(ncell(ras)) # this is to give data so that g becomes a SPDataFrame
g<-as(ras,"SpatialPixelsDataFrame")
proj4string(g)<-proj4string(r)

w<-c(500,2500,5000,10000) # this is the buffer width, values over 2000 start to take a long time to cumpute
wn<-c("_p","_m","_g","_tg")

pl<-lapply(w,function(i){
  spChFIDs(gBuffer(SpatialPoints(coordinates(g)),width=i,byid=TRUE),paste(1:nrow(g),i,sep="_"))
})

p<-do.call("rbind",pl)
proj4string(p)<-proj4string(r)

v<-velox(r)
e<-v$extract(p,fun=function(x){paste(x,collapse="_")})

ans<-lapply(e,function(i){
  x<-table(unlist(strsplit(i,"_")))
  x<-x[names(x)!=1] # ignore background values in raster
  s<-setdiff(labi,names(x))
  if(length(s)==0){
    x<-x[order(names(x))]
  }else{
    temp<-rep(0,length(s))
    names(temp)<-s
    x<-c(x,temp)
    x<-x[order(names(x))]
  }
  100*x/sum(x)
})
land<-as.data.frame(do.call("rbind",ans))
names(land)<-gsub(" ","_",lab[match(names(land),labi)])

la<-split(land,rep(1:length(w),each=nrow(g)))
la<-lapply(1:length(w),function(i){
  temp<-la[[i]]
  names(temp)<-paste0(names(temp),wn[i])
  temp
})
la<-do.call("cbind",la)

g@data<-la
g$Cat_Typ<-"C"

rland<-stack(g)
levelplot(rland,col.regions=rev(terrain.colors(100)),cuts=99)


####################################################################################
### add bridge branch et al
####################################################################################

###
### check how are frag metrics calculated in relation to background
### 

mspa<-"MSPA_for_mature_forest_include_bajos_secondary.tif"
mspa <- raster(mspa)

l<-list(Core=c(117,17),Islet=c(109,9),Perforation=c(105,5),Edge=c(103,3),Loop=c(165,65),Loop=c(167,67),Loop=c(169,69),Bridge=c(133,33),Bridge=c(135,35),Bridge=c(137,37),Branch=c(101,1),Background=c(100,0),Missing=129)

mn<-unlist(l,use.names=FALSE)
lab<-rep(names(l),times=sapply(l,length))
labi<-as.integer(factor(lab))
mat<-cbind(mn,labi)
mspa <- reclassify(mspa, mat)

w<-c(500,2500,5000,10000) # this is the buffer width, values over 2000 start to take a long time to cumpute
wn<-c("_p","_m","_g","_tg")
pl<-lapply(w,function(i){
  spChFIDs(gBuffer(SpatialPoints(coordinates(g)),width=i,byid=TRUE),paste(1:nrow(g),i,sep="_"))
})
p<-do.call("rbind",pl)
proj4string(p)<-proj4string(r)

v<-velox(mspa)

### check what to do with NAs in calculations

e<-v$extract(p,fun=function(x){paste(x,collapse="_")})
co<-lapply(e,function(i){
  table(unlist(strsplit(i,"_")),useNA="always")
})
n<-unique(unlist(sapply(co,names)))
co2<-lapply(co,function(i){
  val<-setdiff(n,names(i)) 
  if(length(val)>0){
    res<-c(i,rep(0,length(val)))
    names(res)[(length(i)+1):length(res)]<-val 
  }else{
    res<-i  
  }
  res<-res[order(names(res))]
  res
})
co2<-do.call("rbind",co2)

frag<-co2/rowSums(co2)
dimnames(frag)[[2]]<-lab[match(dimnames(frag)[[2]],labi)]
frag<-as.data.frame(frag)

la<-split(frag,rep(1:length(w),each=nrow(g)))
la<-lapply(1:length(w),function(i){
  temp<-la[[i]]
  names(temp)<-paste0(names(temp),wn[i])
  temp
})
la<-do.call("cbind",la)

g@data<-la
g$Cat_Typ<-"C"

rfrag<-stack(g)
levelplot(rfrag,col.regions=rev(terrain.colors(100)),cuts=99)

#rfrag<-ras[[1]]
#lr<-lapply(1:ncol(frag),function(i){
#  rfrag<-ras[[1]]
#  rfrag<-setValues(rfrag,frag[,i])
#  rfrag
#})
#rfrag<-stack(lr)
#names(rfrag)<-dimnames(frag)[[2]]
#rfrag<-subset(rfrag,1:8) # no missing values
#plot(rfrag)

#e<-raster:::extract(mspa,b[1:10])

#########################################
### dist and dens to patch
#########################################

# on prend sf, car plus rapide et il n'y a pas de bugs d'indexage

outshp<-st_read("C:/Users/rouf1703/Documents",layer="outshp")
dis<-as.vector(st_area(outshp))/1000/1000

largepatch<-outshp[dis>1500,]

xy<-xyFromCell(rpred,1:ncell(rpred),spatial=TRUE)

g1<-apply(st_distance(largepatch,st_as_sf(xy)),2,min)
ds$dis<-g1

rdens<-rpred
rdens<-setValues(rdens,g1)
names(rdens)<-"Dist_HabP"

#########################################################################
### Plot predictions
#########################################################################

### build total prediction raster

rvar<-stack(rland,rfrag,rroad,rpop,rdens)

### run chosen model with raster covariates

mform<-ml[[16]][[2]]
mform<-Attack~Cat_Typ+Den__MinPS_g+Dist_HabP+Bridge_p+Branch_p
mform<-Attack~Cat_Typ+Branch_p+Dist_HabP+Bridge_p+Den__MinPS_g
#mform<-Attack~Cat_Typ+Bridge_p+Branch_p
#mform<-Attack~Cat_Typ+Dist_Road

fit<-glgm(mform, 
            data=ds,
            grid=50,
            #grid=rpred,
            covariates=rvar, 
            family="binomial", 
            buffer=10000,
            shape=1,
            priorCI=list(sd=c(0.4,4),range=c(2000,50000))
)

lo<-calc(fit$raster$predict.0.025quant,fun=inla.link.invlogit)
me<-calc(fit$raster$predict.mean,fun=inla.link.invlogit)
up<-calc(fit$raster$predict.0.975quant,fun=inla.link.invlogit)

pr<-stack(lo,me,up)
names(pr)<-c("lower","mean","upper")


levelplot(pr,col.regions=colo.scale(1:100,c("black","darkred","red","gold","white")),at=seq(0,1,by=0.01))#+spplot(ds,"Attack",pch=1,col.regions=c("blue","green"))

#tmap_mode("view")
#tm_shape(me)+tm_raster(alpha=0.7,palette=colo.scale(1:10,c("black","darkred","red","gold","white")),n=7)+tm_shape(ds)+tm_dots("Attack",palette=c("green","blue"))




### visualize predictions (static)

mm<-glm(mform,data=ds@data,family=binomial) 
p<-predict(mm,data.frame(as.matrix(rvar),Cat_Typ="B"),type="response")
tempr<-setValues(rpred,p)
#plot(tempr)

### 
bp<-brewer.pal(5,"RdBu")
bp[3]<-"#4dac26"
cols<-colo.scale(1:100,bp)
par(mfrow=c(3,2))
plot(tempr,col=cols,main="glm")
plot(ds,add=TRUE,pch=1)
plot(fit$raster[["predict.invlogit"]],col=cols,main="glgm")
plot(ds,add=TRUE,pch=1)
plot(calc(fit$raster[["predict.mean"]],inv.logit),col=cols,main="glgm")
plot(calc(fit$raster[["predict.0.5quant"]],inv.logit),col=cols,main="glgm")
plot(calc(fit$raster[["predict.0.025quant"]],inv.logit),col=cols,main="glgm")
plot(calc(fit$raster[["predict.0.975quant"]],inv.logit),col=cols,main="glgm")

test<-stack(
  setExtent(tempr,extent(fit$raster[["predict.invlogit"]])),
  fit$raster[["predict.invlogit"]],
  calc(fit$raster[["predict.mean"]],inla.link.invlogit),
  calc(fit$raster[["predict.0.5quant"]],inla.link.invlogit),
  calc(fit$raster[["predict.0.025quant"]],inla.link.invlogit),
  calc(fit$raster[["predict.0.975quant"]],inla.link.invlogit)
)
levelplot(test,col.regions=colo.scale(100,c("white","lightgoldenrod","orange","red","darkred")),cuts=99)

### visualisation prediction (dynamic)
tmap_mode("view")
tm_shape(fit$raster[["predict.invlogit"]])+tm_raster(alpha=0.6,palette=colo.scale(1:10,bp),n=7)+tm_shape(ds)+tm_dots("Attack")



##############################################
####### prediction graphs
##############################################

### temp pred

p<-stack(calc(m[[13]][[1]]$raster$predict.0.025quant,inla.link.invlogit),
         calc(m[[13]][[1]]$raster$predict.mean,inla.link.invlogit),
         calc(m[[13]][[1]]$raster$predict.0.975quant,inla.link.invlogit))

levelplot(p,col.regions=terrain.colors(100),cuts=99)

ds2<-ds@data
ds2$x<-coordinates(ds)[,1]
ds2$y<-coordinates(ds)[,2]

#mform<-ml[[16]][[2]]
mform<-Attack~Cat_Typ+Den__MinPS_g+Dist_HabP+Bridge_p+Branch_p
nvar<-"Branch_p"
xvar<-seq(min(ds2[,nvar]),max(ds2[,nvar]),length.out=50)
allvars<-all.vars(mform)
vars<-allvars[!allvars%in%c("Cat_Typ","Attack",nvar)]
newdat<-with(ds2,data.frame(Cat_Typ="B",as.data.frame(as.list(colMeans(ds2[,vars,drop=FALSE])))))
newdat<-newdat[rep(1,length(xvar)),]
newdat[,nvar]<-xvar

newdat$x<-233173.6
newdat$y<-2037153

#plot(r)
#plot(ds,add=TRUE)
#points(newdat$x[1],newdat$y[2],pch=16,cex=5,col="yellow")

### Explore variogram from glm residuals from model3

#plot(calc(fit$raster$predict.mean-fit$raster$random.mean,fun=inla.link.invlogit))

ds2<-rbind.fill(ds2,newdat)
coordinates(ds2)<-~x+y
proj4string(ds2)<-proj4string(ds)


fit<-glgm(mform, 
          data=ds2,
          grid=20,
          covariates=NULL, 
          family="binomial", 
          buffer=10000,
          shape=1,
          priorCI=list(sd=c(0.4,4),range=c(2000,50000)),
          control.predictor=list(compute=TRUE,link=1)
)

#summary(fit$inla)
p<-inla.link.invlogit(fit$inla$summary.linear.predictor[102:nrow(ds2),])
p<-fit$inla$summary.fitted.values[102:nrow(ds2),]
plot(xvar,p$mean,ylim=0:1,type="l",xlab=nvar,ylab="Probability of attack")
lines(xvar,p$"0.025quant",lty=2)
lines(xvar,p$"0.975quant",lty=2)

glm1<-glm(mform,data=d,family="binomial")
summary(glm1)
visreg(glm1,scale="response")

### posteriors param
autoplot(fit$inla,which=1)


######################################################################
### Explore variogram from glm residuals from model3
######################################################################

glm1<-glm(Attack~Secondary+Cat_Typ+Bridge_p+Branch_p+index_pop8+Dist_Road,data=d,family="binomial")

coords<-as.matrix(d[,c("X","Y")])
v<-variog(coords=coords,data=resid(glm1),breaks=seq(0,50000,by=1000))
fitv<-variofit(v,ini.cov.pars=c(2,5000),cov.model="matern",fix.nugget=TRUE,nugget=0.7,fix.kappa=TRUE,kappa=1)
plot(v)
lines(fitv)

### explore patch size dist habp


patch<-raster("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Mature_forest_2015_include_bajos_secondary.tif")
par(mfrow=c(1,2))
plot(patch)
plot(ds,cex=5*ds$Dist_HabP/max(ds$Dist_HabP),add=TRUE,pch=1,col="white")
plot(patch,col=c("white","lightgreen"))
plot(ds,cex=5*ds$Den__MinPS_tg/max(ds$Den__MinPS_tg),add=TRUE,pch=1,col="darkred")

patch<-aggregate(patch,50,fun=mean)
test<-rasterToPolygons(patch,dissolve=TRUE)
test<-rasterToContour(patch)

