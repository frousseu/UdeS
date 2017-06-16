library(sp)
library(rgdal)
library(rgeos)

# lire les shp
x1<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="poblacion")
x2<-readOGR("C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc",layer="carreteras")

# lire le raster
r<-"C:/Users/rouf1703/Documents/UdeS/Consultation/M-LLecuyer/Doc/Landcover_2015_extended.tif"
r <- stack(r)

# créer un raster de prédictions et le projeter
ras <- raster(xmn=bbox(r)[1,1],xmx=bbox(r)[1,2],ymn=bbox(r)[2,1],ymx=bbox(r)[2,2],ncols=100,nrows=100)
proj4string(ras)<-proj4string(r)

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
ras<-stack(ras1,ras2)

# faire une petite fonction pour visualiser les résultats
foo<-function(i){
  plot(x1,add=TRUE,pch=1)
  plot(x2,add=TRUE,lwd=2)
}

# visualiser les résultats
plot(ras,addfun=foo)

