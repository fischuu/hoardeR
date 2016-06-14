summary.fa <- function(object, ...){

   nCharObj <- nchar(object)
   cat("Summary of fa object\n")
   cat("---------------\n")
   cat("Sequences      :",length(object),"\n")
   cat("Minimum length :",min(nCharObj),"\n")
   cat("1st quartile   :",quantile(nCharObj, 0.25),"\n")
   cat("Median length  :",median(nCharObj),"\n")
   cat("Average length :",mean(nCharObj),"\n")
   cat("3rd quartile   :",quantile(nCharObj, 0.75),"\n")
   cat("Maximum length :",max(nCharObj),"\n")
   invisible(object)

} 
