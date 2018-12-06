
#####################################################################################
### function to create a data.frame for predictions for each variable in a data.frame

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
    if(length(l)==0){
      res<-data.frame(val,stringsAsFactors=FALSE)
    }else{
      res<-data.frame(val,as.data.frame(l),stringsAsFactors=FALSE)
    }
    names(res)[1]<-i
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