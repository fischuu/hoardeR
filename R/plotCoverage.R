plotCoverage <- function(x, use.sqrt=TRUE){

  maxXValue <- max(unlist(sapply(sapply(x,"[",1),"[",1)))
  maxYValue <- max(unlist(sapply(sapply(x,"[",1),"[",2)))
  
  if(use.sqrt) maxYValue <- sqrt(maxYValue)
  
  # If the input is not of class coverageDensity it is assumed to be a path that leeds to a set of bam files  
  if(class(x)!="coverageDensity") x <- coverageDensity(x)

  par(mfrow=c(ceiling(sqrt(length(x))), floor(sqrt(length(x))) ),
      oma = c(5,4,0,0) + 0.1,
      mar = c(0,0,1,1) + 0.1)
  # Get y-range
  for(chrRun in 1:length(x)){
    curRange <- c()
    for(sampleRun in 1:length(x[[chrRun]])){
      curRange <- c(curRange,range(x[[chrRun]][[sampleRun]]$y))
    }
    
    plotX <- x[[chrRun]][[1]]$x 
    plotY <- x[[chrRun]][[1]]$y

    if(use.sqrt) plotY <- sqrt(plotY)
    
        
    plot(x=plotX,
         y=plotY, 
         lwd=1,
         type="l",
         xaxt="n",
         yaxt="n",
         xlim=c(0,maxXValue),
         ylim=c(0,maxYValue))
    
    for(i in 2:length(x[[1]])){
      plotX <- x[[chrRun]][[i]]$x 
      plotY <- x[[chrRun]][[i]]$y
      
      if(use.sqrt) plotY <- sqrt(plotY)
      
      lines(plotX, plotY, lwd=1)
    }    
  }
}
  
