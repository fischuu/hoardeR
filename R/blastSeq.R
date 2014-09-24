# TODO:
# Keep the timing values and return them also as a result

blastSeq <- function(seq, n_blast=20, delay_req=3, delay_rid=60, email=NULL, xmlFolder=NULL, keepInMemory=TRUE, database="chromosome", verbose=TRUE){

# Polite system sleeps as requested from NCBI  
  if(delay_req<3) stop("Sending more requests than once every 3 seconds is considered to be rude from NCBI!")  
  if(delay_rid<60) stop("Polling RID results more often than once per minute is considered to be rude from NCBI!")  
  if(is.null(email)) stop("NCBI requires to provide an email address, please give one in the function call!")

# Getting some needed variables
  totalSeq <- length(seq)
  if(is.null(names(seq))) names(seq) <- 1:totalSeq

# Check about the xml folder settings
  writeXML <- FALSE
  if(!is.null(xmlFolder)){
    writeXML <- TRUE
    dir.create(xmlFolder, showWarnings = FALSE)
  }
  
# Store here the blast RIDs
  RID <- rep(0,totalSeq)
# Store here the results
  res <- list()
  sendThis <- 1
  curRunning <- 0
  ready <- 0
  active <- NULL
  times <- c()
  timeAvg <- c("00:00:00")
# This is very optimistic, maybe I should take also a time break, in case one result doesn't get ready
  while(ready < totalSeq){
    if((curRunning < n_blast) & (sendThis <= totalSeq)){
      Sys.sleep(delay_req)
      RID[sendThis] <- sendFA(seq[sendThis],email=email, database=database)
      curRunning <- curRunning + 1    
      active <- c(active,sendThis)
      sendThis <- sendThis + 1
    } else {
      Sys.sleep(delay_rid) 
      for(i in active){
        res[[i]] <- getBlastResult(RID[i])
        Sys.sleep(delay_req)
        if(res[[i]]$ready==TRUE){
           ready <- ready + 1
           curRunning <- curRunning - 1
           active <- active[-which(active==i)]
           times <- c(times,res$time)
           timeAvg <- timeStat(times)
        # Write here then the XML file to the folder
           if(writeXML){
             file.create(paste(xmlFolder,names(seq)[i],".xml",sep=""))
             fileConn <- file(paste(xmlFolder,names(seq)[i],".xml",sep=""))
             writeLines(res[[i]]$blastRes, fileConn)
             close(fileConn)         
           }
           if(!keepInMemory){
             res[[i]] <- NULL
           }
        }
      }
    }
    if(verbose){
      cat("Missing:",totalSeq-ready,"\n")
      cat("Running:",active,"\n")
      cat("Finished:",ready,"\n")
      cat("Avg. Blast Time:",timeAvg,"\n")
      cat("Total running time: ------ \n")
      cat("---------------------------------------------------------------\n")
    }
  }
 result <- list(RID=RID, res=res, database=database)
 class(result) <- "blastRes"
 result
}
