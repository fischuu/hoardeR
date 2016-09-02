importBed <- function(file, header=FALSE){

  if(header){
    out <- read.table(file=file, row.names=FALSE, header=TRUE, sep="\t")    
  } else {
    out <- read.table(file=file, row.names=FALSE, header=FALSE, sep="\t")    
  }

  out
}