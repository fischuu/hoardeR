importBed <- function(file, header=FALSE){

  .Deprecated("GenomicTools.fileHandler::importBed", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
  if(header){
    out <- read.table(file=file, row.names=FALSE, header=TRUE, sep="\t")    
  } else {
    out <- read.table(file=file, row.names=FALSE, header=FALSE, sep="\t")    
  }

  out
}