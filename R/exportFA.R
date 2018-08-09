exportFA <- function(fa, file=NULL){

  .Deprecated("GenomicTools.fileHandler::exportFA", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
  faNames <- names(fa)
  if(is.null(file)){
    file <- deparse(substitute(fa))
    file <- paste(file,".fa",sep="")
    cat("No file name (option: 'file') given, use the variable name instead:", file, "\n")  
  }
  con <- file(file, "w")
  writeLines(mixVectors(faNames,fa),con)
  close(con)
}