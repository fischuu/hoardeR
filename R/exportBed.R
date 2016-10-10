exportBed <- function(x, file=NULL, header=FALSE){

  if(is.null(file)){
    file <- deparse(substitute(x))
    file <- paste(file,".bed",sep="")
    cat("No file name (option: 'file') given, use the variable name instead:", file, "\n")  
  }
  
  if(header){
    cat("write out the bed-file, using the following column names:\n")
    cat(colnames(x))
    write.table(x, file=file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)    
  } else {
    write.table(x, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)    
  }

}