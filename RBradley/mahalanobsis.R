
library(readxl)



d<-read_excel("C:/Users/rouf1703/Documents/UdeS/Consultation/RBradley/Doc/PLFA data.xls",sheet=2,skip=6)

### format data
n<-names(d)
n<-n[!n%in%c("Peak name",NA,"NA")]
n<-gsub(" ","",gsub(" - ","_",n))
n<-sapply(strsplit(n,"_"),function(i){paste(i[2],i[1],sep="_")})
n<-rep(n,each=5)
names(d)<-c("peak",n)

l<-lapply(2:ncol(d),function(i){
  if(is.na(names(d)[i])){
    res<-NULL
  }else{
    temp<-unlist(strsplit(names(d)[i],"_"))
    season<-temp[1]
    year<-temp[2]
    plot<-substr(d[,i][1],1,2)
    replicate<-substr(d[,i][1],4,4)
    res<-data.frame(peak=d[-1,1],season,year,plot,replicate,val=d[-1,i],stringsAsFactors=FALSE)
  }
  res
})
l<-l[-((length(l)-3):(length(l)))]
d<-do.call("rbind",l)



### help file
x <- matrix(rnorm(10*3), ncol = 3)
stopifnot(mahalanobis(x, 0, diag(ncol(x))) == rowSums(x*x))
##- Here, D^2 = usual squared Euclidean distances
Sx <- cov(x)
D2 <- mahalanobis(x, colMeans(x), Sx)




