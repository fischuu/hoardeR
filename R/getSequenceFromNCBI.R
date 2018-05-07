getSequenceFromNCBI <- function(id, file=NULL){
  out <- c()
  for(idRun in 1:length(id)){
    out[idRun] <- getURL(paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",id[idRun],"&rettype=fasta&retmode=text", sep=""))    
    if(idRun%%10==0) message("Processed",idRun,"queries.")
    if(idRun==length(id)) message("All queries processed!")
  }

  if(!is.null(file)){
    fileConn<-file(file)
    writeLines(out, fileConn)
    close(fileConn)
  } else {
    out    
  }
}