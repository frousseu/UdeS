
#### Aims of prog:
## Reshape the data
#		- Detect at which lines starts and ends each coumpound dataset
#		- Save them in separated files
#
## Derive statistics
#		- mean, sd for each compound
#		- mean, sd for 3 samplings among ...

#### Load packages
library(data.table)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
options(max.print = 500)

#### Set working directory
setwd("C:/Users/rouf1703/Documents/UdeS/Consultation/LLahens")
path<-"C:/Users/rouf1703/Documents/UdeS/Consultation/LLahens/Fichiers_txt"

#### Tool functions
## Split sample text to 3 parts to create groups of measurements
splitFct = function(strVec)
{
	tsn_dt = data.table(lake = character(length(strVec)),
		depth = character(length(strVec)), id = integer(length(strVec))) # tsn = Taxonomic Serial Number

	tsn_dt[, lake := stri_sub(str = strVec,
		to = stri_locate_first(str = strVec, regex = "[:alpha:]")[,1] - 2)]

	tsn_dt[, depth := stri_sub(str = strVec,
		from = stri_locate_first(str = strVec, regex = "[:alpha:]")[,1],
		to = stri_locate_last(str = strVec, regex = "[:alpha:]")[,1])]

	tsn_dt[, id := stri_sub(str = strVec,
		from = stri_locate_last(str = strVec, regex = "[:alpha:]")[,1] + 1)]

	return (tsn_dt)
}

## Check if indices there is the same number of row for each compound.
checkIndices = function(dt)
{
	if(length(unique(dt[, endLine - startLine])) > 1)
	{
		print("*** Warning ***: Non regular file")
		return (list(check = FALSE, nrow = unique(dt[, endLine - startLine])+1))
	}
	return (list(check = TRUE, nrow = unique(dt[, endLine - startLine])+1))
}

## To change some names automatically. Needs to be in the same order
rename = function(df, colNamesToChange, newColNames)
{
	if (length(colNamesToChange) != length(newColNames))
		stop("*** Error (from rename) ***: colNamesToChange and newColNames should be the same size")
	if (FALSE %in% (colNamesToChange %in% names(df)))
		warning("*** Warning (from rename) ***: there are some unknown colNames")

	names(df)[names(df) %in% colNamesToChange] = newColNames

	return (names(df))
}

#####################################################################
#####################################################################
### Declare variables
## File to open

files<-list.files(path)
res<-vector(mode="list",length=length(files))


for(j in seq_along(files)){

  filepath = file.path(path,files[j])

  l<-readLines(filepath)

  g1<-grep("Compound",substr(l,1,8))
  #g2<-grep("Name",substr(l,1,8))
  g2<-g1+2 # finds the header section after each compound
  
  nbCompound = length(g1) # determine number of Compounds assuming the first 8 chars represents a compound block
  indices = data.table(compound = character(nbCompound))

  name = stri_sub(l[g1], from = stri_locate_first(l[g1], regex = ":")[,1])
  name = stri_replace_all(name, replacement = "", regex = ":  ")
  fend = tail(l,1)=="" | all(strsplit(tail(l,1),"")[[1]] %in% "\t") # some files end with an empty line either "" or "\t\t\t..."
  indices[,compound:=name][,startLine:=g2+1][,endLine:=c(g1[-1]-2,length(l)+ifelse(fend,-1,0))] # assumes that there is a single empty line before and after each compound name
  
  #print(checkIndices(indices)[["nrow"]])

  if (checkIndices(indices)[["check"]])
  {
  	n = checkIndices(indices)[["nrow"]]
  	comps<-vector(mode="list",length=nbCompound)
  	for (i in 1:nbCompound)
  	{
  		data = fread(filepath, skip = indices[i, startLine - 2], nrows = n,
  			select = c("Sample Text", "Type", "ng/L", "Conc. Dev. Flagged", "Sig/Noise Flag", "CD Flag"))

  		names(data) = rename(df = data,
  			colNamesToChange = c("Sample Text", "ng/L", "Conc. Dev. Flagged", "Sig/Noise Flag", "CD Flag"),
  			newColNames = c("sampleText", "ppb", "concentrationFlag", "noiseFlag", "cdFlag"))
  
  		data[stri_detect(data[, sampleText], regex = "^[:digit:][:digit:]\\-[:digit:][:digit:][:digit:]"),
  			c("lake", "depth", "id") := splitFct(sampleText)]
  
  		data[, group := ifelse(is.na(lake) | is.na(depth), NA, paste0(lake, "_", depth))]
  		data[!is.na(group), sampleInGroup := .N, by = group]
  
  		data[!is.na(group), c("mean_ppb", "sd_ppb") := .(mean(ppb, na.rm = TRUE), sd(ppb, na.rm = TRUE)), by = group]
  		data[!is.na(group), nonNA_sample := sum(!is.na(ppb)), by = group]
  
  		results = unique(data[!is.na(group), .(group, mean_ppb, sd_ppb, sampleInGroup, nonNA_sample)])
  		
  		### add CB CH (assumes CB and CH are not found anywhere in sample names!!!
  		results[, CB := data$concentrationFlag[grep("CB",substr(data$sampleText,1,2))[which.max(1/(grep("CB",substr(data$sampleText,1,2))-match(group,data$group)))]], by=group]
  		results[, CH := data$concentrationFlag[grep("CH",substr(data$sampleText,1,2))[which.max(1/(grep("CH",substr(data$sampleText,1,2))-match(group,data$group)))]], by=group]
  		

      
  		results[,compound:=name[i]]
  		results[,file:=files[j]]
  		
  		comps[[i]]<-results
  	}
  }
  res[[j]]<-rbindlist(comps)
  print(j)
  
}

d<-rbindlist(res)
ss<-strsplit(d$file,"_")
d[,date:=as.Date(sapply(ss,"[",1),"%Y%m%d")]
d[,type:=sapply(ss,"[",3)]
d[,changed:=ifelse(grepl("Dchanged",file),"yes","no")]
d[,compound:=gsub("\t","",compound)] # some "\t" are staying in compound names

dups<-duplicated(d[,c("group","compound")]) | duplicated(d[,c("group","compound")],fromLast=TRUE)
table(dups)

d<-d[which(!dups | (dups & d$changed=="yes")),]

fwrite(d,"data_compounds.csv",row.names=FALSE)




