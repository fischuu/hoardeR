intersectXMLAnnot <- function(tabSpecies, annot, level="gene", flanking=NULL){
  start <- min(c(tabSpecies$hitStart,tabSpecies$hitEnd))
  end <-  max(c(tabSpecies$hitStart,tabSpecies$hitEnd))
  if(!is.null(flanking)){
    if(length(flanking)==1) flanking <- c(flanking, flanking)
    start <- start - flanking[1]
    end <- end + flanking[2]
  }
  tempAnnot <- annot[annot$V3==level,]
  tempAnnot <- tempAnnot[tempAnnot$V1==tabSpecies$hitChr]
  tempAnnot <- tempAnnot[(tempAnnot$V4>=start) & (tempAnnot$V5<=end),]
  tempAnnot$V10 <- tabSpecies$origChr
  tempAnnot$V11 <- tabSpecies$origStart
  tempAnnot$V12 <- tabSpecies$origEnd
  tempAnnot
}