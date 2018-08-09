summary.fq <- function(object, ...){

  .Deprecated("GenomicTools.fileHandler::summary.fq", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
   nCharObj <- nchar(object)
   cat("Summary of fq object\n")
   cat("--------------------\n")
   cat("Sequences      :",length(object),"\n")
   cat("Minimum length :",min(nCharObj),"\n")
   cat("1st quartile   :",quantile(nCharObj, 0.25),"\n")
   cat("Median length  :",median(nCharObj),"\n")
   cat("Average length :",mean(nCharObj),"\n")
   cat("3rd quartile   :",quantile(nCharObj, 0.75),"\n")
   cat("Maximum length :",max(nCharObj),"\n")
   invisible(object)

} 
