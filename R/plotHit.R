plotHit <- function(hits, flanking=1, window=NULL, annot=TRUE, coverage=FALSE, labels=FALSE,
                    smoothPara=NULL, exonMatch=NULL, diagonal=NULL, verbose=TRUE, output=FALSE,
                    annotFile=NULL, hitSpecies=NULL, hitSpeciesVersion=NULL, origSpecies=NULL,
                    origSpeciesVersion=NULL, fastaFolder=NULL, release=84,
                    origAnnot=NULL, hitAnnot=NULL, nTick=5, which=NULL, figureFolder=NULL,
                    figurePrefix=NULL){
  
# Store the original values
  windowOrig <- window
  hitsOrig <- hits
  origAnnotOrig <- origAnnot
  hitAnnotOrig <- hitAnnot
  if(!is.null(which)) hits <- hits[which,]
  for(hitRun in 1:nrow(hits)){
  # Restore for each run the original values
    window <- windowOrig
    origAnnot <- origAnnotOrig
    hitAnnot <- hitAnnotOrig
    tempGeneName <-gsub(' \"', '', strsplit(strsplit(hits[hitRun]$V9,"gene_name")[[1]][2],'\";')[[1]][1])
    # If the previous version hasn't found a gene name, try another way
    if(is.na(tempGeneName)){
      tempGeneName <- strsplit(strsplit(hits[hitRun]$V9,"gene=")[[1]][[2]],";")[[1]][[1]]
    }

    if(!is.null(figureFolder)) png(file=paste(figureFolder,figurePrefix,hitRun,"-",tempGeneName,".png",sep=""), width=2000, height=1400)
    if(verbose) cat("Start to create figure for",tempGeneName,"(",date(),")\n")
  # Import the fasta files
  # ADD HERE STILL THE OPTION FOR AN OWN FASTA FILE!!!
    if(is.null(fastaFolder)){
      stop("No directory with fasta file given!")    
    } else {
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
    }
    
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
      ORIGcoord <- c(origStart,origEnd)
      HITcoord <- c(hitStart,hitEnd)
      ORIGneg <- FALSE
      HITneg <- FALSE    
      if(ORIGcoord[1]>ORIGcoord[2]){
        ORIGneg <- TRUE
        origStartFlank <- origStart + 1000 * flanking
        origEndFlank <- origEnd - 1000 * flanking
      } else {
        origStartFlank <- origStart - 1000 * flanking
        origEndFlank <- origEnd + 1000 * flanking
      }
      if(HITcoord[1]>HITcoord[2]){
        HITneg <- TRUE
        hitStartFlank <- hitStart + 1000 * flanking
        hitEndFlank <- hitEnd - 1000 * flanking  
      } else{
        hitStartFlank <- hitStart - 1000 * flanking
        hitEndFlank <- hitEnd + 1000 * flanking        
      }
    
      
    # Get the sequence coordiantes
      if(is.null(window)) window <- origEnd-origStart
      if(verbose) cat("Using window size:", window,"\n")
      origLeft <- seq(origStart,origStartFlank,-window)
      origRight <- seq(origEnd,origEndFlank,window)
      if(verbose) cat("Slots (width/window) in area:", length(origLeft) + length(origRight) + 1,"\n")
      
    # Get the plotting coordinates
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
      if(is.numeric(diagonal)){
        ORIGwinSeqLeft <- c()
        HITwinSeqLeft <- c()
        ORIGwinSeqRight <- c()
        HITwinSeqRight <- c()
        if(verbose) cat("Get the left sequences\n")
        if(verbose) pb   <- txtProgressBar(1, length(scoresLeft), style=3, width=60)
        for(i in 1:length(scoresLeft)){
          ORIGwinSeqLeft[i] <- paste(seqOrig[[1]][origLeft[length(origLeft)+1-i]:origLeft[length(origLeft)-i]],collapse="")
          HITwinSeqLeft[i] <- paste(seqHit[[1]][hitLeft[length(hitLeft)+1-i]:hitLeft[length(hitLeft)-i]] ,collapse="")
          if(HITneg) HITwinSeqLeft[i] <- paste(rev(comp(strsplit(HITwinSeqLeft[i],"")[[1]])),collapse="")
          if(verbose) setTxtProgressBar(pb, i)
        }
        if(verbose) cat("\nGet the right sequences\n")
        if(verbose) pb   <- txtProgressBar(1, length(scoresRight), style=3, width=60)
        for(i in 1:length(scoresRight)){
          ORIGwinSeqRight[i] <- paste(seqOrig[[1]][origRight[i]:origRight[1+i]], collapse="")  
          HITwinSeqRight[i] <- paste(seqHit[[1]][hitRight[i]:hitRight[i+1]],collapse="")
          if(HITneg) HITwinSeqRight[i] <- paste(rev(comp(strsplit(HITwinSeqRight[i],"")[[1]])),collapse="")
          if(verbose) setTxtProgressBar(pb, i)
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
        if(verbose) cat("\nCalculate the scores of",nrow(scoreMatrix),"times",nrow(scoreMatrix),"=",nrow(scoreMatrix)^2,"combinations \n")
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
  ## CHECK HERE STILL WITH THE SCORING; WHY ARE THOSE WEIGHTS THERE!?!?!?
      scoresHit <- pairwiseAlignment(toupper(origSequence), 
                                     toupper(hitSequence), 
                                     gapOpening=-1, 
                                     gapExtension=-1, 
                                     substitutionMatrix = mat,
                                     type="local")@score
      sch <- scoresHit / (window)
      
  #    sch <- (sch + 101)/2
      
    } else {
      stop("Currently only diagonal plots are supported!")
#       for(i in 1:length(scoresLeft)){
#         if(HITneg){
#           scoresLeft[i] <- getMatch(SSChr=SSChr, SSStart=SSLeft[length(SSLeft)+1-i], SSEnd=SSLeft[length(SSLeft)-i],
#                                     hitChr=hitChr, hitStart=hitRight[length(hitRight)+1-i], hitEnd=hitRight[length(hitRight)-i], hitOrga=hitOrga)@score      
#         } else {
#           scoresLeft[i] <- getMatch(SSChr=SSChr, SSStart=SSLeft[i], SSEnd=SSLeft[i+1],
#                                     hitChr=hitChr, hitStart=hitLeft[i], hitEnd=hitLeft[i+1], hitOrga=hitOrga)@score
#         }
#       }
#       for(i in 1:length(scoresRight)){
#         if(HITneg){
#           scoresRight[i] <- getMatch(SSChr=SSChr, SSStart=SSRight[i], SSEnd=SSRight[1+i],
#                                      hitChr=hitChr, hitStart=hitLeft[length(hitLeft)-i], hitEnd=hitLeft[length(hitLeft)-i+1], hitOrga=hitOrga)@score      
#         } else {
#           scoresRight[i] <- getMatch(SSChr=SSChr, SSStart=SSRight[i], SSEnd=SSRight[i+1],
#                                      hitChr=hitChr, hitStart=hitRight[i], hitEnd=hitRight[i+1], hitOrga=hitOrga)@score
#         }
#       }
#       scoresHit <- getMatch(SSChr=SSChr, SSStart=SSStart, SSEnd=SSEnd,
#                             hitChr=hitChr, hitStart=hitStart, hitEnd=hitEnd, hitOrga=hitOrga)@score      
#       
#       
#       scl <- round(scoresLeft / (2*window)*100)
#       scr <- round(scoresRight / (2*window)*100)
#       sch <- round(scoresHit / (2*window)*100)
#       
#       scl <- (scl + 101)/2
#       scr <- (scr + 101)/2
#       sch <- (sch + 101)/2
    }
    
    
    plot(-100,-100, ylim=yRange, xlim=c(xRange), xaxt="n", yaxt="n", ylab="", xlab="Chromosomal positions of original and target organism")
    
    
    redGreenFun <- colorRampPalette(c("red", "green"))
    redGreen <- redGreenFun(99)
    
    if(is.numeric(diagonal)){
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
          # DO HERE STILL THE ALIGNMENT TO GET THE CORRECT COLOR!!!
          polygon(c(otherHits$origStart[otherRun], otherHits$origEnd[otherRun], otherHits$hitEnd[otherRun], otherHits$hitStart[otherRun]),c(1,1,-1,-1),col="black", border=NA)
        }
      }
            
    } else{
      # Plot the hitbox
      rect(origStart,-1,origEnd,1,col=redGreen[sch], border=NA)
      # Plot the Left box
      for(i in 1:length(scoresLeft)){
        rect(origLeft[i],-1,origLeft[i+1],1,col=redGreen[scl[i]], border=NA)
      }
      # Plot the right Box
      for(i in 1:length(scoresRight)){
        rect(origRight[i],-1,origRight[i+1],1,col=redGreen[scr[i]], border=NA)
      }
      # Now add the scores
      if(labels){
        text(mean(origStart,origEnd)+window/2, 0, srt=90, adj = 0, labels = round(scoresHit), xpd = TRUE)
        for(i in 1:length(scoresLeft)){
          text(mean(origLeft[i],origLeft[i+1])-window/2, 0, srt=90, adj = 0, labels = round(scoresLeft[i]), xpd = TRUE )
        }
        for(i in 1:length(scoresRight)){
          text(mean(origRight[i],origRight[i+1])+window/2, 0, srt=90, adj = 0, labels = round(scoresRight[i]), xpd = TRUE )
        }  
      }    
    }
    
    
    lines(c(xRange[1],xRange[2]),c(1,1))
    origTicks <- seq(xRange[1], xRange[2],length.out = nTick)
    for(tick in 1:nTick){
      lines(c(origTicks[tick],origTicks[tick]),c(1,1.05))
    }
    origTicksNice <- prettyNum(as.character(round(origTicks,0)), big.mark=",")
    text(origTicks,1.1,paste(origChr,origTicksNice,sep=":"))
    
    lines(c(xRange[1],xRange[2]),c(-1,-1))
    hitTicks <- seq(hitStartFlank, hitEndFlank,length.out = nTick)
    for(tick in 1:nTick){
      lines(c(origTicks[tick],origTicks[tick]),c(-1,-1.05))
    }
    hitTicksNice <- prettyNum(as.character(round(hitTicks,0)), big.mark=",")
    text(origTicks,-1.1,paste(hitChr,hitTicksNice,sep=":"))
    
    lines(c(origStart, origEnd),c(1,1),lwd=5)
    lines(c(origStart, origEnd),c(-1,-1),lwd=5) 
  #  if(HITneg){
  #    lines(c(origStart, hitStart-hitEnd+origStart),c(-1,-1),lwd=5)     
  #  } else {
  #    lines(c(origStart, hitEnd-hitStart+origStart),c(-1,-1),lwd=5)
  #  }
    
    if(HITneg){
      axis(2,at=c(-1,1),label=c("Neg. strand","Pos. strand"))
    } else {
      axis(2,at=c(-1,1),label=c("Pos. strand","Pos. strand"))
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
  #      if(is.numeric(exonMatch)){
  #        SSseq <- getSequence(SSChr, xRange[1], xRange[2], "pig", FALSE) 
  #      }
        for(i in 1:nrow(hitAnnot)){
          lines(c(hitAnnot$V4[i]+ xRange[1] - hitRange[1], hitAnnot$V5[i] + xRange[1]-hitRange[1]), c(-1.5,-1.5), lwd=5)
          tempGeneName <-gsub(' \"', '', strsplit(strsplit(hitAnnot$V9[i],"gene_name")[[1]][2],'\";')[[1]][1])
        # If the previous version hasn't found a gene name, try another way
          if(is.na(tempGeneName)){
            tempGeneName <- strsplit(strsplit(hitAnnot$V9[i],"gene=")[[1]][[2]],";")[[1]][[1]]
          }
          text(mean(c(hitAnnot$V4[i]+ xRange[1] - hitRange[1], hitAnnot$V5[i] + xRange[1]-hitRange[1])), -1.6, tempGeneName) 
          
    #      temp <- as.character(HITexon$gene_id)
    #      if(substr(temp[1],1,1)==" ") temp <- gsub(" ","",temp) 
    #      tempExon <- HITexon[temp==as.character(HITgenes$geneID[i]),]
    #      if(nrow(tempExon)>0){
    #        # Here we plot the bars for the exon - that means, here we can get also the exon similarity values:
    #        for(j in 1:nrow(tempExon)){
    #          # Extract the temporary exon sequence in the hit:
    #          if(is.numeric(exonMatch)){
    #            if(tempExon$V3[j]=="exon"){
    #              exonSeq <- getSequence(tempExon$V1[j], tempExon$V4[j], tempExon$V5[j], hitOrga, HITneg)
    #              exonSeqMatch <- pairwiseAlignment(SSseq, exonSeq, gapOpening=-2, gapExtension=-2, type="local")
    #              exonSeqRatio <- exonSeqMatch@score/(nchar(exonSeq)*2)
    #              if(exonSeqRatio>exonMatch){
    #                xpoly1 <- tempExon[j,4]+ xRange[1]-hitRange[1] + exonSeqMatch@subject@range@start
    #                xpoly2 <- xRange[1] + exonSeqMatch@pattern@range@start
    #                xpoly3 <- xpoly2 + exonSeqMatch@pattern@range@width
    #                xpoly4 <- xpoly1 + exonSeqMatch@subject@range@width
    #                if(xpoly1>xRange[1] & xpoly4<xRange[2])polygon(c(xpoly1, xpoly2 , xpoly3 ,xpoly4),c(-1.5,-1,-1,-1.5),col=redGreen[min(99,round(exonSeqRatio*100))], border=NA)
    #              }
    #            }
    #          }
    #          lines(c(tempExon[j,4]+ xRange[1]-hitRange[1],tempExon[j,5]+ xRange[1]-hitRange[1]), c(-1.5,-1.5), lwd=5)             
    #        }  
    #     }
          
        }
      }
    } # if(annot)
    
    if(coverage){
      features <- GRanges( seqnames = rep(SSChr,length(SSLeft)+length(SSRight)),
                           ranges = IRanges(seq(xRange[1], xRange[2], window),
                                            width=window)
      )
      olap <- summarizeOverlaps(features, bfl)
      counts12 <- apply(assays(olap)$counts[,1:2],1,mean)
      counts34 <- apply(assays(olap)$counts[,3:4],1,mean)
      maxCounts <- max(c(counts12,counts34))
      counts12 <- counts12/maxCounts
      counts34 <- counts34/maxCounts
      counts12 <- counts12 + 2
      counts34 <- counts34 + 2
      
      centers <- as.matrix(features@ranges)
      centers <- centers[,1] - centers[,2]/2 
      
      if(is.numeric(smoothPara)){
        lo12 <- loess(counts12~centers, span=smoothPara)
        lo34 <- loess(counts34~centers, span=smoothPara)
        lines(centers,predict(lo12), col='red', lwd=2)
        lines(centers,predict(lo34), col='blue', lwd=2)
      } else {
        lines(centers,counts12, col='red', lwd=2)
        lines(centers,counts34, col='blue', lwd=2)    
      }
      
      lines(xRange,c(2,2))
      # lines(xRange,c(3,3))
      # Y-Axis Label
      if(HITneg){
        axis(2,at=c(2,3),label=c("0",maxCounts))
      } else {
        axis(2,at=c(2,3),label=c("0",maxCounts))
      }
    }
    
    if(output) list(SHit=scoresHit, SL=scoresLeft,SR=scoresRight, scoreMat=scoreMatrix)
    if(!is.null(figureFolder)) dev.off()
    if(verbose) cat("\nFigure created (",date(),")\n")
  }
}
