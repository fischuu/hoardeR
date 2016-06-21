tableSpecies <- function(xml, species=NULL, type="chr", minOutput=TRUE, exclude="", locations=FALSE){
  species.int <- species
  if(is.null(species.int)){
    species.int <- hoardeR::species$Scientific.name
  }
  xmlOne <- xml[[1]]
  for(i in 1:length(xml)){
    xmlOne <- rbind(xmlOne,xml[[i]])
  }
  if(type=="chr") xmlOne <- xmlOne[!is.na(xmlOne$hitChr),]
  
  if(locations){
    res <- xmlOne[grepl(species.int,xmlOne[,1]),]
  } else {
    xmlOne <- xmlOne[,1]
    
    res <- c()
    for(i in 1:length(species.int)){
      res[i] <- sum(grepl(species.int[i],xmlOne))
    }
    names(res) <- species.int
    if(minOutput)res <- res[res>0]
    if(!is.null(exclude)){
      remThese <- which(names(res)==exclude)
      if(length(remThese)>0) res <- res[-remThese]
    }
  }
  res
}