cap <- function(x){
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
}

getAnnotation <- function(species=NULL, release=84, version=NULL, annotationFolder=NULL, type="gtf"){
  if(is.null(version)){
    temp <- hoardeR::species
    version <- temp$Ensembl.Assembly[temp$Scientific.name==species]
  }
  
# fread zip support is OS dependend
  os <- "linux"
  if(grepl("Windows", sessionInfo()$running)) os <- "windows"
  
  if(is.null(annotationFolder)) annotationFolder <- getwd()
  
  dir.create(annotationFolder, showWarnings=FALSE)
  species.int <- gsub(" ", "_",tolower(species))
  ensemblURL <- paste("ftp://ftp.ensembl.org/pub/release-",release,"/gtf/",species.int,"/",sep="")
  fileName <-  paste(cap(species.int),".",gsub(" ","",version),".",release,".chr.",type,".gz",sep="")    
  .file = file.path(annotationFolder, fileName)
  
  message("Check if file ",.file," exists ...")
  
  # download file
  if(!file.exists(.file)){
    message("... file wasn't found. Try to download it from Ensembl ftp server.")
    download.file(paste(ensemblURL,fileName,sep=""), .file)
  } else {
    message("... found!") 
  }
  
# Now create the input string, depending on the os
  if(os=="linux"){
   inputString <-  paste('zcat',.file)
  } else if(os=="windows"){
    inputString <- paste("gzip -dc",.file)
  }
  
  if(type=="gtf"){
    temp <- fread(input = inputString, skip=5, colClasses = c("character",
                                                              "character",
                                                              "character",
                                                              "integer",
                                                              "integer",
                                                              "character",
                                                              "character",
                                                              "character",
                                                              "character"))
  }else if(type=="gff"){
    temp <- fread(input = inputString, skip=0, colClasses = c("character",
                                                              "character",
                                                              "character",
                                                              "integer",
                                                              "integer",
                                                              "character",
                                                              "character",
                                                              "character",
                                                              "character"))
  }
  temp
}