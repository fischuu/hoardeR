exportFA <- function(fa, file=NULL){
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