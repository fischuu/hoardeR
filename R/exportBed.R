exportBed <- function(x, file, header=FALSE){

  if(header){
    cat("write out the bed-file, using the following column names:\n")
    cat(colnames(x))
    write.table(x, file=file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)    
  } else {
    write.table(x, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)    
  }

}