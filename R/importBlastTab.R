importBlastTab <- function(file){

  .Deprecated("GenomicTools.fileHandler::importBlastTab", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
  out <- read.table(file,stringsAsFactors=FALSE, sep="\t")
  out
}