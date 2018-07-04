plotCoverage <- function(x){

  # If the input is not of class coverageDensity it is assumed to be a path that leeds to a set of bam files  
  if(class(x)!="coverageDensity") x <- coverageDensity(x)

  par(mfrow=c(length(x),1))
  # Get y-range
  for(chrRun in 1:length(x)){
    curRange <- c()
    for(sampleRun in 1:length(x[[chrRun]])){
      curRange <- c(curRange,range(x[[chrRun]][[sampleRun]]$y))
    }
    ylim <- range(curRange)
    
    plot(x[[chrRun]][[1]],
         ylim = ylim,
         lwd=1)
    
    for(i in 2:length(x[[1]])){
      lines(x[[chrRun]][[i]], lwd=1)
    }    
  }
}
  
