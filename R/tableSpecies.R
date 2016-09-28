tableSpecies <- function(xml, species=NULL, type="chr", minOutput=TRUE, exclude="", locations=FALSE){

  species.int <- species
  if(is.null(species.int)){
    species.int <- hoardeR::species$Scientific.name
  }

  if(locations){
    xmlOne <- xml[[1]]
    origLoc <- names(xml)[1] 
    xmlOne$origChr <- gsub(">","",strsplit(origLoc,":")[[1]][1])
    xmlOne$origStart <- strsplit(strsplit(origLoc,":")[[1]][2],"-")[[1]][1]
    xmlOne$origEnd <- strsplit(strsplit(origLoc,":")[[1]][2],"-")[[1]][2]
    
    if(length(xml)>=2){
      for(i in 2:length(xml)){
        tempXML <- xml[[i]]
        origLoc <- names(xml)[i] 
        if(nrow(tempXML)>0){
          tempXML$origChr <- gsub(">","",strsplit(origLoc,":")[[1]][1])
          tempXML$origStart <- strsplit(strsplit(origLoc,":")[[1]][2],"-")[[1]][1]
          tempXML$origEnd <- strsplit(strsplit(origLoc,":")[[1]][2],"-")[[1]][2]
          xmlOne <- rbind(xmlOne,tempXML)
        }
      }
    }
    
    if(type=="chr") xmlOne <- xmlOne[!is.na(xmlOne$hitChr),]
    
    res <- xmlOne[grepl(species.int,xmlOne[,1]),]
  } else {
    xmlOne <- xml[[1]]
    if(length(xml)>=2){
      for(i in 2:length(xml)){
        xmlOne <- rbind(xmlOne,xml[[i]])
      }
    }
    
    if(type=="chr") xmlOne <- xmlOne[!is.na(xmlOne$hitChr),]
    
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