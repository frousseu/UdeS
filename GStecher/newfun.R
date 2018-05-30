

toseq<-function(x,n=100){
  if(is.numeric(x)){
    r<-range(x,na.rm=TRUE)
    seq(r[1],r[2],length.out=n)
  }else{
    sort(unique(x))  
  }
}




newdata<-function(x,v=names(x),n=100,fun=mean,list=FALSE){
  
  # returns the mean or the levels
  mm<-function(y){
      if(is.numeric(y)){
        fun(y)
      }else{
        names(rev(sort(table(y))))[1]
      }
  }
  ## possibly a bug with the function did not save the last part
  ans<-lapply(v,function(i){
      if(n==1L){
        val<-mm(x[,i])
      }else{
        val<-toseq(x[,i],n=n)
      }
      l<-lapply(x[,setdiff(names(x),i),drop=FALSE],mm)
      res<-data.frame(val,as.data.frame(l),stringsAsFactors=FALSE)
      names(res)[1]<-i
      if(list){ # inefficient, should not be turned to data.frame if list=TRUE
        as.list(res)
      }else{
        res
      } 
  })  
  names(ans)<-v
  
  if(length(v)==1L){
    unlist(ans,recursive=FALSE)
  }else{
    ans
  }

}


#d<-data.frame(x=1:10,y=runif(10),f=letters[c(1,1,1,2,2,2,2,2,2,2)])

l<-newdata(d)

