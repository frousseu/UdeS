library(rcrossref)
library(dplyr)
library(stringi)
library(stringr)
library(RCurl)


#aut<-c("Marc B\u00E9lisle")
# when a - is not present in the first name, it should not be added (e.g. F. Guillaume Blanchet)

recentPubs<-function(authors=NULL,bold=NULL,keyword=NULL,first=30,date="2018-01-01",html=TRUE){
  
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
  
  #authors<-c(authors,bold)
  
  l<-lapply(authors,function(i){
    a<-unlist(strsplit(i," "))
    tx<-cr_works(query=keyword,filter=c(from_pub_date=date),flq=c(query.author=a[2]),limit=200,sort='published',order="desc")
    x<-as.data.frame(tx$data)
    if(nrow(x)==0L){
      warning(paste("No matches found for",i))
      #browser()
      x
    }else{
      k<-sapply(x[,"author"],function(j){ ### searching for a specific pattern here cause app does fuzzy matching
        
        #if(i=="Fanie Pelletier"){browser()}
        
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
    }
  })
  
  print(paste(authors,sapply(l,nrow)))
  
  l<-l[sapply(l,nrow)!=0L] # Sometimes, a df with 0 rows, 0 cols is returned! take it out here to get common column names
  
  lnames<-Reduce(intersect,lapply(l,names)) # do.call not working !?
  l<-lapply(l,function(i){i[,lnames]}) # for some reason, columns retained are not always consistant
  df<-do.call("rbind",l)
  df<-df[!duplicated(df$DOI),] # remove because of authors together
  df<-df[rev(order(df$created)),]
  if(html && nrow(df)){
    refs_orig<-lapply(df$DOI,function(i){
      m<-tryCatch(
        cr_cn(dois = i, format = "text", style = "apa")
        ,error=function(j){TRUE}
      )
      if(!isTRUE(m)){
        m
      }else{
        NULL
      }
    }) # cr_cn gets nothing out if it can't find a DOI, so use sapply to know which one is missing
    refs_orig<-unlist(lapply(refs_orig,function(i){if(is.null(i)){NA}else{i}}),use.names=FALSE)
    #print missing refs
    cat("MISSING REFS \n\n\n")
    print(df[is.na(unlist(refs_orig)),])
    cat("\n\n\n\n\n\n")
    
    prof<-string2name(authors)
    stud<-string2name(bold)
    name<-unique(c(prof,stud))
    bname<-paste0("<b>",name,"</b>")
    names(bname)<-name
    refs<-str_replace_all(refs_orig,bname)
    stopifnot(nrow(df)==length(refs)) # safety check, not sure if cr_cn gets everything each time (returns nothing if nothing found)
    #elim<-!duplicated(refs)
    #refs<-refs[elim]
    #url<-df$URL[elim]
    refs<-strsplit(refs,". doi:")
    refs<-unique(Map(c,refs,df$URL))
    refs<-sapply(refs,function(i){
      paste0("<p>",i[1]," <a href=",i[3],">",i[2],"</a>","</p>")
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

x<-recentPubs(date="2018-01-01")





