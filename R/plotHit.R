plotHit <- function(hits, flanking=1, window=NULL, annot=TRUE, coverage=FALSE,
                    smoothPara=NULL, diagonal=0.25, verbose=TRUE, output=FALSE,
                    hitSpecies=NULL, hitSpeciesVersion=NULL, origSpecies=NULL,
                    origSpeciesVersion=NULL, fastaFolder=NULL, release=84,
                    origAnnot=NULL, hitAnnot=NULL, nTick=5, which=NULL, figureFolder=NULL,
                    figurePrefix=NULL, indexOffset=0, bamFolder=NULL, bamFiles=NULL, groupIndex=NULL, groupColor=NULL,
                    countWindow=NULL){
  
# Store the original values
  windowOrig <- window
  hitsOrig <- hits
  origAnnotOrig <- origAnnot
  hitAnnotOrig <- hitAnnot
  
# Input checks
#  is(!is.numeric(diagonal)) stop("The option 'diagonal' needs to be numeric between 0 and 1.")
  
  if(!is.null(bamFolder)){
    if(is.null(bamFiles)){
      fls <- list.files(bamFolder, recursive=TRUE, pattern="\\.bam$", full.names=TRUE)
    } else {
      fls <- paste(bamFolder, bamFiles, sep="")
    }
  }
  
  if(coverage){
    if(is.null(groupIndex)) groupIndex <- rep(1,length(fls))
    if(is.null(groupColor)) groupColor <- 1:max(groupIndex)
  }
  
  if(!is.null(which)) hits <- hits[which,]
  for(hitRun in 1:nrow(hits)){
  # Restore for each run the original values
    window <- windowOrig
    origAnnot <- origAnnotOrig
    hitAnnot <- hitAnnotOrig
    flankingLeft.orig <- flanking
    flankingRight.orig <- flanking
    flankingLeft.hit <- flanking
    flankingRight.hit <- flanking
    
    if(grepl("gene_name", hits[hitRun]$V9)){
      tempGeneName <-gsub(' \"', '', strsplit(strsplit(hits[hitRun]$V9,"gene_name")[[1]][2],'\";')[[1]][1])
    } else if(grepl("gene=", hits[hitRun]$V9)){
      tempGeneName <- strsplit(strsplit(hits[hitRun]$V9,"gene=")[[1]][[2]],";")[[1]][[1]]
    } else if(grepl("gene_id", hits[hitRun]$V9)){
      tempGeneName <- gsub(' \"', '', strsplit(strsplit(hits[hitRun]$V9,"gene_id")[[1]][2],'\";')[[1]][1])
    } else {
      stop("Cannot resolve the hit organism gene name from the annotation file.")
    }

    if(!is.null(figureFolder)) png(filename=file.path(figureFolder,figurePrefix,hitRun+indexOffset,"-",tempGeneName,".png"), width=2000, height=1400)
    if(verbose) cat("Start to create figure for",tempGeneName,"(",date(),")\n")
  # Import the fasta files
  # ADD HERE STILL THE OPTION FOR AN OWN FASTA FILE!!!
    if(is.null(fastaFolder)){
      message("No directory with fasta file given! Use the working directory:\n", getwd())    
      fastaFolder <- getwd()
    }
      if(is.null(origSpeciesVersion)){
        temp <- hoardeR::species
        origSpeciesVersion <- temp$Ensembl.Assembly[temp$Scientific.name==origSpecies]
      }
      if(verbose) cat("Original species version:", origSpeciesVersion,"\n")
      if(is.null(hitSpeciesVersion)){
        temp <- hoardeR::species
        hitSpeciesVersion <- temp$Ensembl.Assembly[temp$Scientific.name==hitSpecies]
      }
      if(verbose) cat("Hit species version:", hitSpeciesVersion,"\n")
      # First the required original fasta file
        dir.create(fastaFolder, showWarnings=FALSE)
        origSpecies.int <- gsub(" ", "_",tolower(origSpecies))
        ensemblURL <- paste("ftp://ftp.ensembl.org/pub/release-",release,"/fasta/",origSpecies.int,"/dna/",sep="")
        fileName <-  paste(cap(origSpecies.int),".",origSpeciesVersion,".dna.chromosome.",hits$origChr[hitRun],".fa.gz",sep="")
        .file = file.path(fastaFolder, fileName)
        if(!file.exists(.file)) download.file(paste(ensemblURL,fileName,sep=""), .file)
        if(verbose) cat("Read original fasta file:\n   ", .file ,"\n")
        seqOrig <- read.fasta(.file,seqtype="DNA")
      # And then the hit species    
        hitSpecies.int <- gsub(" ", "_",tolower(hitSpecies))
        ensemblURL <- paste("ftp://ftp.ensembl.org/pub/release-",release,"/fasta/",hitSpecies.int,"/dna/",sep="")
        fileName <-  paste(cap(hitSpecies.int),".",hitSpeciesVersion,".dna.chromosome.",hits$hitChr[hitRun],".fa.gz",sep="")
        .file = file.path(fastaFolder, fileName)
        if(!file.exists(.file)) download.file(paste(ensemblURL,fileName,sep=""), .file)
        if(verbose) cat("Read hit fasta file:\n   ", .file ,"\n")
        seqHit <- read.fasta(.file,seqtype="DNA")
    
    # Getting the constants
    origChr <- hits$origChr[hitRun]
    origStart <- as.numeric(hits$origStart[hitRun])
    origEnd <- as.numeric(hits$origEnd[hitRun])
    hitChr <- hits$hitChr[hitRun]
    hitStart <- as.numeric(hits$hitStart[hitRun])
    hitEnd <- as.numeric(hits$hitEnd[hitRun]) 
    if(verbose){
      cat("Original loci:\n   ", origChr,":",
          prettyNum(as.character(origStart), big.mark=","),"-",
          prettyNum(as.character(origEnd), big.mark=","),"\n")
      cat("Hit loci:\n   ", hitChr,":",
          prettyNum(as.character(hitStart), big.mark=","),"-",
          prettyNum(as.character(hitEnd), big.mark=","),"\n")
    }
    
    
    # Getting the flanking coordinates
    # Check first, that the flanking isn't getting to the negative side or leaves the chromosome area.
      if(is.null(window)) window <- origEnd-origStart
      if(verbose) cat("Using window size:", window,"\n")
      ORIGcoord <- c(origStart,origEnd)
      HITcoord <- c(hitStart,hitEnd)
      ORIGneg <- FALSE
      HITneg <- FALSE    
      if(ORIGcoord[1]>ORIGcoord[2]){
        ORIGneg <- TRUE
        if(origStart + 1000 * flanking > length(seqOrig[[1]])){    # The original start flankinf site is outside the chromosome boarder
          flankingLeft.orig <- floor(origStart/window)             # Adjust the flanking factor
        } else if(origEnd - 1000 * flanking <= 1){
          flankingRight.orig <- floor(origEnd/window)
        }
        origStartFlank <- origStart + 1000 * flankingLeft.orig
        origEndFlank <- origEnd - 1000 * flankingRight.orig
        origStartFlank.forPlotting <- origStart + 1000 * flanking
        origEndFlank.forPlotting <- origEnd - 1000 * flanking
      } else {
        if(origStart - 1000 * flanking <=0){                       # The original sequence start is outside the chromsome border (with flanking)
          flankingLeft.orig <- floor(origStart/window)
        } else if(origEnd + 1000 * flanking > length(seqOrig[[1]])){
          flankingRight.orig <- floor(origEnd/window)
        }
        origStartFlank <- origStart - 1000 * flankingLeft.orig
        origEndFlank <- origEnd + 1000 * flankingRight.orig
        origStartFlank.forPlotting <- origStart - 1000 * flanking
        origEndFlank.forPlotting <- origEnd + 1000 * flanking
        
      }
      if(HITcoord[1]>HITcoord[2]){
        HITneg <- TRUE
        if(hitStart + 1000 * flanking > length(seqHit[[1]])){
          flankingLeft.hit <- floor(hitStart/window)
        } else if(hitEnd - 1000 * flanking <= 1){
          flankingRight.hit <- floor(hitEnd/window)
        }
        hitStartFlank <- hitStart + 1000 * flankingLeft.hit
        hitEndFlank <- hitEnd - 1000 * flankingRight.hit  
        hitStartFlank.forPlotting <- hitStart + 1000 * flanking
        hitEndFlank.forPlotting <- hitEnd - 1000 * flanking
        
      } else{
        if(hitStart - 1000 * flanking <=0){
          flankingLeft.hit <- floor(hitStart/window)
        } else if(hitEnd + 1000 * flanking > length(seqHit[[1]])){
          flankingRight.hit <- floor(hitEnd/window)
        }
        hitStartFlank <- hitStart - 1000 * flankingLeft.hit
        hitEndFlank <- hitEnd + 1000 * flankingRight.hit        
        hitStartFlank.forPlotting <- hitStart - 1000 * flanking
        hitEndFlank.forPlottinig <- hitEnd + 1000 * flanking
        
      }
    
      
    # Get the sequence coordiantes

      origLeft <- seq(origStart,origStartFlank,-window)
      origRight <- seq(origEnd,origEndFlank,window)
      origLeft.forPlotting <- seq(origStart,origStartFlank.forPlotting,-window)
      origRight.forPlotting <- seq(origEnd,origEndFlank.forPlotting,window)
      
      if(verbose) message("Slots (width/window) in area:", length(origLeft) + length(origRight) + 1)
      
    # Get the plotting coordinates
      xRange.forPlotting <- range(c(origLeft.forPlotting,origRight.forPlotting))
      xRange <- range(c(origLeft,origRight))
      yRange <- c(-1.25,1.25)
      if(annot) yRange <- yRange + c(-0.5,0.5)
      if(coverage) yRange <- yRange + c(0,1.5)
    
    if(HITneg){
      hitLeft <- seq(hitEnd,hitEndFlank,-window)
      hitRight <- seq(hitStart,hitStartFlank,window)      
    } else {
      hitLeft <- seq(hitStart,hitStartFlank,-window)
      hitRight <- seq(hitEnd,hitEndFlank,window)  
    }
    
    # Calculate the matching scores
      scoresLeft <- rep(NA,length(origLeft)-1)
      scoresRight <- rep(NA,length(origRight)-1)
      scoresHit <- NA
      

      ORIGwinSeqLeft <- c()
      HITwinSeqLeft <- c()
      ORIGwinSeqRight <- c()
      HITwinSeqRight <- c()

      if(length(scoresLeft)>1){
        if(verbose) message("Get the left sequences")
        for(i in 1:length(scoresLeft)){
          if(length(origLeft>1)) ORIGwinSeqLeft[i] <- paste(seqOrig[[1]][origLeft[length(origLeft)+1-i]:origLeft[length(origLeft)-i]],collapse="")
          if(length(hitLeft>1)) HITwinSeqLeft[i] <- paste(seqHit[[1]][hitLeft[length(hitLeft)+1-i]:hitLeft[length(hitLeft)-i]] ,collapse="")
          if(HITneg) HITwinSeqLeft[i] <- paste(rev(comp(strsplit(HITwinSeqLeft[i],"")[[1]])),collapse="")
        }
      } else {
        if(verbose) message("No left sequences, original sequence is on chromosome border")
      }
  
      if(length(scoresRight)>1){
         if(verbose) message("Get the right sequences")
         for(i in 1:length(scoresRight)){
           if(length(origRight>1)) ORIGwinSeqRight[i] <- paste(seqOrig[[1]][origRight[i]:origRight[1+i]], collapse="")  
           if(length(hitRight>1)) HITwinSeqRight[i] <- paste(seqHit[[1]][hitRight[i]:hitRight[i+1]],collapse="")
           if(HITneg) HITwinSeqRight[i] <- paste(rev(comp(strsplit(HITwinSeqRight[i],"")[[1]])),collapse="")
         }
      } else {
        if(verbose) message("No right sequences, original sequence is on chromosome border")
      }
      origSequence <- paste(seqOrig[[1]][origStart:origEnd],collapse="")
      hitSequence <- paste(seqHit[[1]][hitStart:hitEnd],collapse="")
        
      if(HITneg) hitSequence <- paste(comp(strsplit(hitSequence,"")[[1]]),collapse="")
        
      ORIGwinSeq <- c(ORIGwinSeqLeft,origSequence,ORIGwinSeqRight)
      HITwinSeq <- c(HITwinSeqLeft,hitSequence, HITwinSeqRight)
      scoreMatrix <- matrix(-1,nrow=length(ORIGwinSeq), ncol=6)
      colnames(scoreMatrix) <- c("ORIGstart","ORIGend", "HITstart", "HITwidth", "score", "hitIndex")
      mat <- nucleotideSubstitutionMatrix(match = 1, 
                                          mismatch = -1,
                                          baseOnly = FALSE)
      if(verbose) cat("Calculate the scores of",nrow(scoreMatrix),"times",nrow(scoreMatrix),"=",nrow(scoreMatrix)^2,"combinations \n")
      if(verbose) pb   <- txtProgressBar(1, nrow(scoreMatrix), style=3, width=60)
     for(origIndRun in 1:nrow(scoreMatrix)){
       bestHIT <- 0
      for(hitIndRun in 1:nrow(scoreMatrix)){
        evalThis <- pairwiseAlignment(toupper(ORIGwinSeq[origIndRun]),
                                      toupper(HITwinSeq[hitIndRun]),
                                      gapOpening=-1,
                                      gapExtension=-1,
                                      substitutionMatrix = mat,
                                      type="local")
        if(evalThis@score/(window)>scoreMatrix[origIndRun,5]){
          best <- evalThis
          bestHit <- hitIndRun
          scoreMatrix[origIndRun,5] <- evalThis@score/(window)
        }
      }
        
      scoreMatrix[origIndRun,1] <- xRange[1] + (origIndRun-1)*window + best@pattern@range@start
      scoreMatrix[origIndRun,2] <- scoreMatrix[origIndRun,1] + best@pattern@range@width
      scoreMatrix[origIndRun,3] <- xRange[1] + (bestHit-1)*window + best@subject@range@start
      scoreMatrix[origIndRun,4] <- scoreMatrix[origIndRun,3] + best@subject@range@width
      scoreMatrix[origIndRun,6] <- bestHit
      if(verbose) setTxtProgressBar(pb, origIndRun)
      
  }

      scoresHit <- pairwiseAlignment(toupper(origSequence), 
                                     toupper(hitSequence), 
                                     gapOpening=-1, 
                                     gapExtension=-1, 
                                     substitutionMatrix = mat,
                                     type="local")@score
      sch <- scoresHit / (window)
      
    plot(-100,-100, ylim=yRange, xlim=c(xRange.forPlotting), xaxt="n", yaxt="n", ylab="", xlab="Chromosomal positions of original and target organism")
    
    
    redGreenFun <- colorRampPalette(c("red", "green"))
    redGreen <- redGreenFun(99)
    
    # Plot the hit polygons
    plotDiag <- scoreMatrix[scoreMatrix[,5]>diagonal,]
    if(nrow(as.matrix(plotDiag))==0){
      warning("No matches found for 'diagonal' filter setting!")
    } else {
      if(is.vector(plotDiag)){
        polygon(c(plotDiag[1], plotDiag[2] , plotDiag[4],plotDiag[3]),c(1,1,-1,-1),col=redGreen[min(99,round((plotDiag[5])*100))], border=NA)
      } else {
        for(diagRun in 1:nrow(plotDiag)){
          polygon(c(plotDiag[diagRun,1], plotDiag[diagRun,2] , plotDiag[diagRun,4],plotDiag[diagRun,3]),c(1,1,-1,-1),col=redGreen[min(99,round((plotDiag[diagRun,5])*100))], border=NA)
        }        
      }
    }
      
      
  # Plot the hitbox
    rect(origStart,-1,origEnd,1,col=redGreen[min(99, round((sch * 100)))], border=NA)
  # Plot the other hits in the window
  # Get first the candidates to plot
    otherHits <- hitsOrig[(hitsOrig$hitChr==hitChr) & (hitsOrig$origChr==origChr),]
    otherHits <- otherHits[-which((otherHits$origStart==origStart) & (otherHits$origEnd==origEnd)),]
    otherHits <- otherHits[((otherHits$hitStart <= hitEndFlank) & (otherHits$hitStart >= hitStartFlank)) | ((otherHits$origStart <= origEndFlank) & (otherHits$origStart >= origStartFlank)) ,]
    otherHits <- otherHits[((otherHits$hitEnd <= hitEndFlank) & (otherHits$hitEnd >= hitStartFlank)) | ((otherHits$origEnd <= origEndFlank) & (otherHits$origEnd >= origStartFlank)),]
  # And then adjust the coordinates for the hit, so that they can be plotted
    if(nrow(otherHits)>0){
      otherHits$hitStart <- otherHits$hitStart - hitStartFlank + origStartFlank
      otherHits$hitEnd <- otherHits$hitEnd - hitStartFlank + origStartFlank
      for(otherRun in 1:nrow(otherHits)){
        polygon(c(otherHits$origStart[otherRun], otherHits$origEnd[otherRun], otherHits$hitEnd[otherRun], otherHits$hitStart[otherRun]),c(1,1,-1,-1),col="black", border=NA)
      }
    }
          
    lines(c(xRange[1],xRange[2]),c(1,1))
    if(xRange[1]>=xRange.forPlotting[1]){
      origTicks <- seq(xRange[1], xRange[2],length.out = nTick)      
    } else {
      origTicks <- seq(xRange[1], xRange[2],length.out = nTick/2+1)
    }

    for(tick in 1:length(origTicks)){
      lines(c(origTicks[tick],origTicks[tick]),c(1,1.05))
    }
    origTicksNice <- prettyNum(as.character(round(origTicks,0)), big.mark=",")
    text(origTicks,1.1,paste(origChr,origTicksNice,sep=":"))
    
    lines(c(xRange.forPlotting[1],xRange.forPlotting[2]),c(-1,-1))
    hitTicks <- seq(hitStartFlank, hitEndFlank,length.out = nTick)
    for(tick in 1:nTick){
      lines(c(origTicks[tick],origTicks[tick]),c(-1,-1.05))
    }
    hitTicksNice <- prettyNum(as.character(round(hitTicks,0)), big.mark=",")
    text(origTicks,-1.1,paste(hitChr,hitTicksNice,sep=":"))
    
    lines(c(origStart, origEnd),c(1,1),lwd=5)
    lines(c(origStart, origEnd),c(-1,-1),lwd=5) 

    if(HITneg){
      axis(2,at=c(-1,1),labels=c("Neg. strand","Pos. strand"))
    } else {
      axis(2,at=c(-1,1),labels=c("Pos. strand","Pos. strand"))
    }
    
    # In case the annot flag is set, add it to the plot
    if(annot){  
      origAnnot <- origAnnot[origAnnot$V1==origChr,]
      origAnnot <- origAnnot[((origAnnot$V4<=xRange[2]) & (origAnnot$V4>=xRange[1]))  ,]
      origAnnot <- origAnnot[((origAnnot$V5<=xRange[2]) & (origAnnot$V5>=xRange[1]))  ,]
      origAnnot <- origAnnot[origAnnot$V3=="gene",]
      if(nrow(origAnnot)[1]>0){
        for(i in 1:nrow(origAnnot)){
          lines(c(origAnnot$V4[i],origAnnot$V5[i]), c(1.5,1.5), lwd=5)
          tempGeneName <-gsub(' \"', '', strsplit(strsplit(origAnnot$V9[i],"gene_name")[[1]][2],'\";')[[1]][1])
          text(mean(c(origAnnot$V4[i],origAnnot$V5[i])), 1.6,tempGeneName) 
        }
      }
  
      hitAnnot <- hitAnnot[hitAnnot$V1==hitChr,]
      if(HITneg){
        hitAnnot <- hitAnnot[((hitAnnot$V4>=hitEndFlank) & (hitAnnot$V4<=hitStartFlank))  ,]
        hitAnnot <- hitAnnot[((hitAnnot$V5>=hitEndFlank) & (hitAnnot$V5<=hitStartFlank))  ,]
        hitAnnot <- hitAnnot[hitAnnot$V3=="gene",]
        hitRange <- range(c(hitLeft,hitRight))
      } else {
        hitAnnot <- hitAnnot[((hitAnnot$V4<=hitEndFlank) & (hitAnnot$V4>=hitStartFlank))  ,]
        hitAnnot <- hitAnnot[((hitAnnot$V5<=hitEndFlank) & (hitAnnot$V5>=hitStartFlank))  ,]
        hitAnnot <- hitAnnot[hitAnnot$V3=="gene",]
        hitRange <- range(c(hitLeft,hitRight))
      }
      
      if(nrow(hitAnnot)>0){

        for(i in 1:nrow(hitAnnot)){
          lines(c(hitAnnot$V4[i]+ xRange[1] - hitRange[1], hitAnnot$V5[i] + xRange[1]-hitRange[1]), c(-1.5,-1.5), lwd=5)

          if(grepl("gene_name", hitAnnot$V9[i])){
            tempGeneName <-gsub(' \"', '', strsplit(strsplit(hitAnnot$V9[i],"gene_name")[[1]][2],'\";')[[1]][1])
          } else if(grepl("gene=", hitAnnot$V9[i])){
            tempGeneName <- strsplit(strsplit(hitAnnot$V9[i],"gene=")[[1]][[2]],";")[[1]][[1]]
          } else if(grepl("gene_id", hitAnnot$V9[i])){
            tempGeneName <- gsub(' \"', '', strsplit(strsplit(hitAnnot$V9[i],"gene_id")[[1]][2],'\";')[[1]][1])
          } else {
            stop("Cannot resolve the hit organism gene name from the annotation file.")
          }

          text(mean(c(hitAnnot$V4[i]+ xRange[1] - hitRange[1], hitAnnot$V5[i] + xRange[1]-hitRange[1])), -1.6, tempGeneName) 
          
        }
      }
    } # if(annot)

    if(coverage){

      if(is.null(countWindow)) countWindow <- window
      if(verbose) message("\nWe use",length(fls),"bamfiles to calculate the average coverage.")
      tmpStart <- seq(xRange[1],xRange[2],countWindow)
      tmpStart <- tmpStart[-length(tmpStart)]
      features <- GRanges( seqnames = origChr,
                           ranges = IRanges(tmpStart,
                                            width=countWindow))

      counts <- matrix(0, ncol=length(tmpStart), nrow=length(fls))
      for(bamRun in 1:length(fls)){
        counts[bamRun,] <- bamCount(fls[bamRun], features, verbose=FALSE)        
      }

      countsGroups <- list()
      for(groupRun in 1:max(groupIndex)){
        countsGroups[[groupRun]] <- counts[groupIndex==groupRun,]
      }
      
      groupMeans <- matrix(0, ncol=length(tmpStart), nrow=max(groupIndex))
      for(groupRun in 1:max(groupIndex)){
        if(is.null(nrow(countsGroups[[groupRun]]))){
          groupMeans[groupRun,] <- countsGroups[[groupRun]]
        } else {
          groupMeans[groupRun,] <- apply(countsGroups[[groupRun]],2,mean)
        }
      }
      
      maxCounts <- max(groupMeans)
      groupMeans <- groupMeans/maxCounts
      groupMeans <- groupMeans + 2

#      centers <- as.matrix(features@ranges)
      centers <- (2*features@ranges@start + features@ranges@width)/2 
      
      if(is.numeric(smoothPara)){
        for(groupRun in 1:max(groupIndex)){
          tmpCurve <- loess(groupMeans[groupRun,]~centers, span=smoothPara)
          lines(centers,predict(tmpCurve), col=groupColor[groupRun], lwd=2)
        }
      } else {
        for(groupRun in 1:max(groupIndex)){
          lines(centers,groupMeans[groupRun,], col=groupColor[groupRun], lwd=2)
        }
      }
      
      lines(xRange,c(2,2))
      # lines(xRange,c(3,3))
      # Y-Axis Label
      if(HITneg){
        axis(2,at=c(2,3),labels=c("0",maxCounts))
      } else {
        axis(2,at=c(2,3),labels=c("0",maxCounts))
      }
    }
    
    if(output) list(SHit=scoresHit, SL=scoresLeft,SR=scoresRight, scoreMat=scoreMatrix)
    if(!is.null(figureFolder)) dev.off()
    if(verbose) message("\nFigure created (",date(),")")
  }
}
