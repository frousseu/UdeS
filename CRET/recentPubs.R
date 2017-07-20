library(rcrossref)
library(dplyr)
library(stringi)
library(stringr)
library(RCurl)


#aut<-c("Marc B\u00E9lisle")

recentPubs<-function(authors=NULL,bold=NULL,keyword=NULL,first=30,date="2016-01-01",html=TRUE){
  
  ### this function turns strings to names with capital letters or to abbreviated names in the scientific paper style
  string2name<-function(x,abbrev=TRUE){
    xl<-stri_trans_totitle(tolower(x))
    if(abbrev){
      xl<-strsplit(xl," ")
      first<-sapply(xl,function(i){
        ii<-unlist(strsplit(i[1],"-"))
        paste(paste0(substr(ii,1,1),"."),collapse="-")
      })
      second<-sapply(xl,function(i){
        paste(i[2:length(i)],collapse=" ")
      })
      paste0(second,", ",first)
    }else{
      xl  
    }
  }
  
  if(is.null(authors)){
    authors<-c("Marc B\u00E9lisle","Marc Belisle","Mark Vellend","Marco Festa-Bianchet","Fanie Pelletier","Robert Bradley","Bill Shipley","Dominique Gravel","Dany Garant","Sophie Calm\u00E9","Sophie Calme")  
  }
  if(is.null(bold)){ # if no names are given to bold, the list of students is taken from the list
    g<-getURL("https://docs.google.com/spreadsheets/d/1E6yTZJlZqvwxFESTn0xheB1KS3a7L9PXAFfHc8QXP7w/pub?gid=0&single=true&output=csv",.encoding="UTF-8")
    bold<-read.csv(text=g,header=TRUE,stringsAsFactors=FALSE)$Nom
  }else{
    bold<-authors # tkae this out redundant  
  }
  
  l<-lapply(authors,function(i){
    a<-unlist(strsplit(i," "))
    tx<-cr_works(query=keyword,filter=c(from_pub_date=date),flq=c(query.author=a[2]),limit=200,sort='published',order="desc")
    x<-as.data.frame(tx$data)
    if(nrow(x)==0L){
      stop("No matches found")
    }
    k<-sapply(x[,"author"],function(j){ ### searching for a specific pattern here cause app does fuzzy matching
      if(any(!c("given","family")%in%names(j))){ # sometimes some names are missing and not consistant
        return(TRUE)
      }
      g1<-grep(a[1],pull(j,"given")) # dplyr ugly column extraction... 
      g2<-grep(a[2],pull(j,"family"))
      g<-base:::intersect(g1,g2)
      if(any(g)){
        TRUE
      }else{
        FALSE  
      }
    })
    x<-x[k,]
    if(nrow(x)){
      x<-x[1:min(first,nrow(x)),]
    }
    x
  })
  
  lnames<-Reduce(intersect,lapply(l,names)) # do.call not working !?
  l<-lapply(l,function(i){i[,lnames]}) # for some reason, columns retained are not always consistant
  df<-do.call("rbind",l)
  df<-df[rev(order(df$created)),]
  if(html){
    refs_orig<-unlist(cr_cn(dois = df$DOI, format = "text", style = "apa"))
    prof<-string2name(authors)
    stud<-string2name(bold)
    name<-unique(c(prof,stud))
    bname<-paste0("<b>",name,"</b>")
    names(bname)<-name
    refs<-str_replace_all(refs_orig,bname)
    stopifnot(nrow(df)==length(refs)) # safety check, not sure if cr_cn gets everything each time
    elim<-!duplicated(refs)
    refs<-refs[elim]
    url<-df$URL[elim]
    refs<-strsplit(refs,". doi")
    refs<-unique(Map(c,refs,url))
    refs<-sapply(refs,function(i){
      paste0("<p>",i[1]," <a href=",i[3],">","doi",i[2],"</a>","</p>")
    })
    invisible(sapply(refs,function(i){
      cat(i,"\n\n")
    }))
    html<-refs
  }else{
    html<-NULL  
  }
  ans<-list(df=df,html=html)
  ans
}

x<-recentPubs(date="2016-01-01")














