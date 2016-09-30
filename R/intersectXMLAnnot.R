intersectXMLAnnot <- function(tabSpecies, annot, level="gene", flanking=NULL){
  res <- list()
  for(i in 1:nrow(tabSpecies)){
    res[[i]] <- intersectXMLAnnot.internal(tabSpecies=tabSpecies[i,], annot=annot, level=level, flanking=flanking)
  }
  res <- do.call(rbind, res)
  res
}

intersectXMLAnnot.internal <- function(tabSpecies, annot, level, flanking){
  start <- min(c(tabSpecies$hitStart,tabSpecies$hitEnd))
  end <-  max(c(tabSpecies$hitStart,tabSpecies$hitEnd))
  if(!is.null(flanking)){
    flanking <- 1000 * flanking
    if(length(flanking)==1) flanking <- c(flanking, flanking)
    start <- start - flanking[1]
    end <- end + flanking[2]
  }
  tempAnnot <- annot[annot$V3==level,]
  tempAnnot <- tempAnnot[tempAnnot$V1==tabSpecies$hitChr]
  tempAnnot <- tempAnnot[(tempAnnot$V4>=start) & (tempAnnot$V5<=end),]
  tempAnnot$origChr <- tabSpecies$origChr
  tempAnnot$origStart <- tabSpecies$origStart
  tempAnnot$origEnd <- tabSpecies$origEnd
  tempAnnot$hitChr <- tabSpecies$hitChr
  tempAnnot$hitStart <- tabSpecies$hitStart
  tempAnnot$hitEnd <- tabSpecies$hitEnd
  tempAnnot
}