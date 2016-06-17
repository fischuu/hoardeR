extractTSinfo <- function(x){
  ortholog <- strsplit(x,"target=new>")
  ortholog <- strsplit(ortholog[[1]][2],"</a></td>")[[1]][1]
  
  geneName <- strsplit(x,"<td nowrap>")
  geneName <- strsplit(geneName[[1]][2],"</td>")[[1]][1]
    
  totalSites <- strsplit(x,"<B>") 
  consSites <- strsplit(totalSites[[1]][2],"</B>")[[1]][1]
  poorlySites <- strsplit(totalSites[[1]][3],"</B>")[[1]][1]   
  
  res <- data.frame(Ortholog=ortholog,
                    geneName=geneName,
                    consSites=consSites,
                    poorlySites=poorlySites)
  res
}

targetScan <- function(mirna=NULL, species=NULL, release="7.1", maxOut=NULL){

 if(is.null(species)){
  # Guessing species from mirna
    if(substr(mirna,1,3)=="hsa") species <- "Human"
    if(substr(mirna,1,3)=="mmu") species <- "Mouse"
    if(substr(mirna,1,3)=="rno") species <- "Rat"
    if(substr(mirna,1,3)=="ptr") species <- "Chimpanzee"
    if(substr(mirna,1,3)=="mml") species <- "Rhesus"
    if(substr(mirna,1,3)=="bta") species <- "Cow"
    if(substr(mirna,1,3)=="cfa") species <- "Dog"
    if(substr(mirna,1,3)=="mdo") species <- "Opossum"
    if(substr(mirna,1,3)=="gga") species <- "Chicken"
    if(substr(mirna,1,3)=="xtr") species <- "Frog"
  }
  
# Input checks
  if(is.null(mirna)) stop("No mirne name given. Use e.g. 'miR-9-5p'.")
  species <- paste(toupper(substr(species, 1, 1)), tolower(substr(species, 2, nchar(species))), sep="")
  species <- match.arg(species, c("Human", "Mouse", "Rat", "Chimpanzee", "Rhesus", "Cow", "Dog", "Opossum", "Chicken", "Frog"))
  release <- gsub("\\.","",release)
 
# Retrieve the content 
  tsAddress <- paste("http://www.targetscan.org/cgi-bin/targetscan/vert_71/targetscan.cgi?species=",species,"&mirg=",mirna,sep="")
  tsOut <- scan(tsAddress, what = "", sep = "\n", quiet = TRUE)

# Check first if the miRNA is in the targetScan database  
  if(sum(grepl("is not in our miRNA database",tsOut[1:min(100,length(tsOut))]))>0){
       warning(mirna," is not in the targetScan database!\n")
       res <- data.frame(Ortholog=NULL,
                         geneName=NULL,
                         consSites=NULL,
                         poorlySites=NULL)
       res
    
  }else {
      # Check first, if the targetScen result is unique
        if(sum(grepl("matches multiple families in our miRNA database",tsOut[1:min(100,length(tsOut))]))>0){
          multFams <- tsOut[grepl("=miR",tsOut)]
          newMirnas <- character(length(multFams))
          for(i in 1:length(multFams)){
            temp <- strsplit(multFams[i], "</A")[[1]][1]
            newMirnas[i] <- strsplit(temp,'\">')[[1]][2]
          }
          warning("Multiple matches multiple families in the targetScan database for ",mirna,":\n",paste(newMirnas,collapse="; "),"\nOnly the first one is used!")
          mirna <- newMirnas[1]  
          tsAddress <- paste("http://www.targetscan.org/cgi-bin/targetscan/vert_71/targetscan.cgi?species=",species,"&mirg=",mirna,sep="")
          tsOut <- scan(tsAddress, what = "", sep = "\n", quiet = TRUE)
        }
      
      # Find the rows of interest (Assume it to be in the first 100 rows, if this isn't the case extent the search area)
        startRow <- grepl("<th>total</th>",tsOut[1:min(100,length(tsOut))])
        if(sum(startRow)!=1)  startRow <- grepl("<th>total</th>",tsOut)
        if(sum(startRow)!=1) stop("ERROR: No table provided by targetScan.org!")
        startRow <- which(startRow==1)
        ifelse(is.null(maxOut), maxOut <- length(tsOut)-1, maxOut <- startRow + maxOut - 1)
        
        
      # Now extract the information and put them into a dataframe
      
      # The first row is a bit different, as it is contained in the header row, all others are then standardized
        temp1 <- strsplit(tsOut[startRow],"<td>")
        firstEntry <- paste(temp1[[1]][2],"<td>",temp1[[1]][3],sep="")
        res <- extractTSinfo(firstEntry)
        for(i in (startRow+1):maxOut){
          res <- rbind(res,extractTSinfo(tsOut[i]))
        }
  }
  res
}