# TODO:
# Keep the timing values and return them also as a result

blastSeq <- function(seq, n_blast=20, delay_req=3, delay_rid=60, email=NULL, xmlFolder=NULL, logFolder=NULL, keepInMemory=FALSE, database="refseq_genomes", verbose=TRUE, createLog=TRUE){

# Developer version variable
  useLast <- FALSE
  
  startTime <- Sys.time()
  firstRun <- TRUE
  
# Polite system sleeps as requested from NCBI  
  if(delay_req<10) stop("Sending more requests than once every 10 seconds is considered to be rude from NCBI!")  
  if(delay_rid<60) stop("Polling RID results more often than once per minute is considered to be rude from NCBI!")  
  if(is.null(email)) stop("NCBI requires to provide an email address, please give one in the function call!")

# Getting some needed variables
  totalSeq <- length(seq)
  if(is.null(names(seq))) names(seq) <- 1:totalSeq

# Check about the xml folder settings
  writeXML <- TRUE
  if(is.null(xmlFolder)){
   xmlFolder <- paste("hoardeR-",format(Sys.time(), "%d.%m.%Y@%X"),sep="") 
   xmlFolder <- gsub(":","-", xmlFolder)
  }
  
  dir.create(xmlFolder, showWarnings = FALSE)  
# Set up the log folder
  if(createLog){
    if(is.null(logFolder)){
      if(is.null(xmlFolder)){
        if(useLast){
          if(verbose) cat("The useLast option is set and no log/xml is given. Search in the working directory for the last hoardeR log file...\n")
          allFiles <- list.files(recursive = TRUE)
          allLogs <- allFiles[grepl("\\.log",allFiles)]
          timeDiffs <- Sys.time() - file.info(allLogs)$mtime
          takeThisLog <- allLogs[which.min(timeDiffs)]
          xmlFolder <- file.path(getwd(),strsplit(takeThisLog,"/logs")[[1]][1])
          message("Use this folder as xml:", xmlFolder)
          message("Use this folder as log:", file.path(xmlFolder,"logs") )
        } else {
          stop("No log/xml path given.")
        }
      } else {
        #logFolder <- strsplit(xmlFolder,"/")[[1]]
        #logFolder <- logFolder[nchar(logFolder)>0]
        #logFolder <- logFolder[1:(length(logFolder)-1)]
        logFolder <- file.path(xmlFolder,"logs")
        message("Create/use log folder: ", logFolder)
      }
    }
    dir.create(logFolder, showWarnings = FALSE)    
  }

# Create the RID/sequence info table
  seqInfo <- data.frame(seqNames=names(seq),
                        seqRID=rep("0",length(seq)),
                        seqFinished=rep(FALSE,length(seq)),
                        seqRuntime=rep("00:00:00",length(seq)),
                        stringsAsFactors=FALSE)

# Read/Write the RID/sequence info table
  if(createLog){
    if(file.exists(file.path(logFolder,"seqRID-info.csv"))){
      seqInfoImported <- read.table(file.path(logFolder,"seqRID-info.csv"), header=TRUE, stringsAsFactors=FALSE)
      if(sum(seqInfo$seqNames==seqInfoImported$seqNames)!=nrow(seqInfo)) stop("Log-file mismatch! Please provide the right log file for the existing project or change the logFolder option.")
      seqInfo <- seqInfoImported
    } else {
      write.table(seqInfo,file=file.path(logFolder,"seqRID-info.csv"), quote=FALSE, row.names=FALSE)      
    }
  }

# Write the project Log
  if(createLog){
    if(file.exists(file.path(logFolder,"hoardeR.log"))){
      fileConn<-file(file.path(logFolder,"hoardeR.log"))
      logLines <- readLines(fileConn)
      close(fileConn)      
      startTime <- as.POSIXct(logLines[1])
    } else {  
      fileConn<-file(file.path(logFolder,"hoardeR.log"))
      writeLines(c(as.character(startTime),
                   "---------------------------------",
                   "Settings:",
                   paste("n_blast",n_blast)),fileConn)
      close(fileConn)      
      
    }

  }

# Store here the blast RIDs
  RID <- rep(0,totalSeq)

# Store here the results (In case no log file is created, otherwise get the last stored result)
  if(createLog){
  # Adjust this here still according to 'keepInMemory' settings
    res <- list()
    sendThis <- min(which(seqInfo$seqRID==0))
    curRunning <- sum((seqInfo$seqRID!=0) & seqInfo$seqFinished==FALSE)
    ready <- sum(seqInfo$seqFinished==TRUE)
    active <- which(((seqInfo$seqRID!=0) & seqInfo$seqFinished==FALSE)==TRUE)
    RID <- seqInfo$seqRID
    timeAvg <- c("00:00:00")  
  } else {
    res <- list()
    sendThis <- 1
    curRunning <- 0
    ready <- 0
    active <- NULL
    timeAvg <- c("00:00:00")  
  }
  
# This is very optimistic, maybe I should take also a time break, in case one result doesn't get ready
  while(ready < totalSeq){
    if((curRunning < n_blast) & (sendThis <= totalSeq)){
      Sys.sleep(delay_req)
      RID[sendThis] <- sendFA(seq[sendThis],email=email, database=database, logFolder=logFolder, verbose=verbose)
      curRunning <- curRunning + 1    
      active <- c(active,sendThis)
      seqInfo$seqRID[sendThis] <- RID[sendThis]
      sendThis <- sendThis + 1
      if(!keepInMemory) write.table(seqInfo,file=file.path(logFolder,"seqRID-info.csv"), quote=FALSE, row.names=FALSE) 
    } else {
      if(!firstRun){
        Sys.sleep(delay_rid)
      } else {
        firstRun <- FALSE 
      }
      for(i in active){
        res[[i]] <- getBlastResult(RID[i]) 
        Sys.sleep(delay_req)
        if(res[[i]]$ready==TRUE){
           ready <- ready + 1
           curRunning <- curRunning - 1
           active <- active[-which(active==i)]
           seqInfo$seqFinished[seqInfo$seqNames==names(seq)[i]] <- TRUE
           timings <- seqInfo$seqRuntime[(seqInfo$seqFinished==TRUE)]
           timeAvg <- timeStat(timings[timings!="0"])
        # Write here then the XML file to the folder
           if(writeXML){
             xmlFile <- paste(names(seq)[i],".xml",sep="")
             xmlFile <- gsub(">","", xmlFile)
             xmlFile <- gsub(":",".", xmlFile)
             file.create(file.path(xmlFolder, xmlFile  ))
             fileConn <- file(file.path(xmlFolder,xmlFile ))
             writeLines(res[[i]]$blastRes, fileConn)
             close(fileConn)  
             write.table(seqInfo,file=file.path(logFolder,"seqRID-info.csv"), quote=FALSE, row.names=FALSE)  
           }
           if(!keepInMemory){
             res[[i]] <- NULL
           }
        } else {
          seqInfo$seqRuntime[seqInfo$seqNames==names(seq)[i]] <- res[[i]]$time
          write.table(seqInfo,file=file.path(logFolder,"seqRID-info.csv"), quote=FALSE, row.names=FALSE)  
        }
      }
    }
    if(verbose){
      timeAvgPrint <- timeAvg
      if(timeAvgPrint=="00:00:00") timeAvgPrint <- "NA"
      cat("Missing:",totalSeq-ready,"\n")
      cat("Running:",length(active),"\n")
      cat("Finished:",ready,"\n")
      cat("Avg. Blast Time:",timeAvgPrint,"\n")
      cat("Total running time:",secToTime(as.numeric(Sys.time() - startTime, units="secs")),"\n")
      cat("---------------------------------------------------------------\n")
    }
  }
 result <- list(RID=RID, res=res, database=database)
 class(result) <- "blastRes"
 result
}
