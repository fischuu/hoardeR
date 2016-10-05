exportFA <- function(fa, file){
  faNames <- names(fa)
  con <- file(file, "w")
  writeLines(mixVectors(faNames,fa),con)
  close(con)
}