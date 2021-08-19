

######################################################################
### function to create a sequence from the range of values in a vector
toseq<-function(x,n=100){
  if(is.numeric(x)){
    r<-range(x,na.rm=TRUE)
    seq(r[1],r[2],length.out=n)
  }else{
    sort(unique(x))  
  }
}

#####################################################################################
### function to create a data.frame for predictions for each variable in a data.frame

newdata<-function(x,v=names(x),n=100,fun=mean,list=FALSE,factors=TRUE){
  
  # returns the mean or the levels
  mm<-function(y){
    if(is.numeric(y)){
      fun(y)
    }else{
      names(rev(sort(table(y))))[1]
      #factor(f,levels=)
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
    if(length(l)==0){
      res<-data.frame(val,stringsAsFactors=TRUE)
    }else{
      res<-data.frame(val,as.data.frame(l),stringsAsFactors=TRUE)
    }
    names(res)[1]<-i
    res<-res[,names(x),drop=FALSE] # put in the same order as
    if(factors){
      w<-which(sapply(x,function(k){is.character(k) | is.factor(k)}))
      res[w]<-lapply(w,function(j){
        factor(res[[j]],levels=levels(x[[j]]))
      })  
    }
    if(list){ # inefficient, should not be turned to data.frame if list=TRUE
      as.list(res)
    }else{
      res
    } 
  })  
  names(ans)<-v
  #browser()
  #if(length(v)==1L){
  #  unlist(ans,recursive=FALSE)
  #}else{
  ans
  #}
}
