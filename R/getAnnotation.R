cap <- function(x){
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
}

getAnnotation <- function(species=NULL, assembly=NULL, annotationFolder=NULL, type="gff3", verbose=TRUE){

  species.DF <- hoardeR::species
  species.cur <- species.DF[species.DF$Organism.Name==species,]
  if(nrow(species.cur)>1) species.cur <- species.cur[1,]
  
  if(is.null(assembly)){
    assembly <- species.cur$Assembly.Name[1]
    if(verbose) cat("No assembly version provided, use the default:",assembly,"\n")
  }
      
# fread zip support is OS dependend
  os <- "linux"
  if(grepl("Windows", sessionInfo()$running)) os <- "windows"
  
  if(is.null(annotationFolder)){
    if(verbose) cat("No directory with annotations files given! Use the working directory: \n", getwd(),"\n")
    annotationFolder <- getwd()
  } else {
    if(substrRight(annotationFolder,1)=="/") annotationFolder <- removeRight(annotationFolder,1)
  }
  
  
  dir.create(annotationFolder, showWarnings=FALSE)
  NCBI.URL <- paste(species.cur$NCBI.Url,"/GFF/ref_",assembly,"_top_level.gff3.gz",sep="")
  NCBI.assemblyTable <- paste(species.cur$NCBI.Url,"/Assembled_chromosomes/chr_accessions_",assembly,sep="")
  fileName <-  paste("ref_",assembly,"_top_level.gff3.gz",sep="")    
  fileName.assemblyTable <- paste("chr_accessions_",assembly,sep="")
  .file = file.path(annotationFolder, fileName)
  .file.assemblyTable = file.path(annotationFolder, fileName.assemblyTable)
  
  if(verbose) cat("Check if file ",.file," exists ... \n")
  
  # download file
  if(!file.exists(.file)){
    if(verbose) cat("... file wasn't found. Try to download it from NCBI ftp server.\n")
    download.file(NCBI.URL, .file, quiet=TRUE)
  }
  
  if(file.size(.file)<1){
    if(verbose) cat("File could not be downloaded. If you want to use the assembly", speciesAssembly, "for", species,"please provide the file: ", .file,"\n")       
  } else {
    if(verbose) cat("... found!\n") 
  }
  
  if(verbose) cat("Check if file ",.file.assemblyTable," exists ... \n")
  
  # download file
  if(!file.exists(.file.assemblyTable)){
    if(verbose) cat("... file wasn't found. Try to download it from NCBI ftp server.\n")
    download.file(NCBI.assemblyTable, .file.assemblyTable, quiet=TRUE)
  }
  
  if(file.size(.file.assemblyTable)<1){
    if(verbose) cat("File could not be downloaded. If you want to use the assembly", speciesAssembly, "for", species,"please provide the file: ", .file.assemblyTable,"\n")       
  } else {
    if(verbose) cat("... found!\n") 
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
  } else if(type=="gff3"){
    accessions <-  read.table(.file.assemblyTable, sep="\t", header=TRUE, comment.char = "", stringsAsFactor=FALSE)
    temp <- importGFF3(.file, chromosomes=accessions[,2])
    for(i in 1:nrow(accessions)){
      temp[V1 == accessions$RefSeq.Accession.version[i], V1 := accessions$X.Chromosome[i]]  
    }
  }
  temp
}