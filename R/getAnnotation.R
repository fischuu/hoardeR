cap <- function(x){
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
}

getAnnotation <- function(species=NULL, release=84, version=NULL, annotationFolder, type="gtf"){
  if(is.null(version)){
    temp <- hoardeR::species
    version <- temp$Ensembl.Assembly[temp$Scientific.name==species]
  }
  
  dir.create(annotationFolder, showWarnings=FALSE)
  species.int <- gsub(" ", "_",tolower(species))
  ensemblURL <- paste("ftp://ftp.ensembl.org/pub/release-",release,"/gtf/",species.int,"/",sep="")
  fileName <-  paste(cap(species.int),".",version,".",release,".chr.",type,".gz",sep="")    
  .file = file.path(annotationFolder, fileName)
  
  # download file
  if(!file.exists(.file)) download.file(paste(ensemblURL,fileName,sep=""), .file)
  
  if(type=="gtf"){
    temp <- fread(input = paste('zcat',.file), skip=5, colClasses = c("character",
                                                                      "character",
                                                                      "character",
                                                                      "integer",
                                                                      "integer",
                                                                      "character",
                                                                      "character",
                                                                      "character",
                                                                      "character"))
  }else if(type=="gff"){
    temp <- fread(input = paste('zcat',.file), skip=0, colClasses = c("character",
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