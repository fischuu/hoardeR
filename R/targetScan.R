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

targetScan <- function(mirna=NULL, species="Human", release="7.1", maxOut=NULL){
# Input checks
  if(is.null(mirna)) stop("No mirne name given. Use e.g. 'miR-9-5p'.")
  species <- paste(toupper(substr(species, 1, 1)), tolower(substr(species, 2, nchar(species))), sep="")
  species <- match.arg(species, c("Human", "Mouse", "Rat", "Chimpanzee", "Rhesus", "Cow", "Dog", "Opossum", "Chicken", "Frog"))
  release <- gsub("\\.","",release)
 
# Retrieve the content 
  tsAddress <- paste("http://www.targetscan.org/cgi-bin/targetscan/vert_71/targetscan.cgi?species=",species,"&mirg=",mirna,sep="")
  tsOut <- scan(temp, what = "", sep = "\n", quiet = TRUE)

# Find the rows of interest (Assume it to be in the first 100 rows, if this isn't the case extent the search area)
  startRow <- substr(tsOut[1:100],0,14)=="<th>total</th>"
  if(sum(startRow)!=1)  startRow <- substr(tsOut,0,14)=="<th>total</th>"
  if(sum(startRow)!=1) stop("ERROR: No table provided by targetScan.org!")
  startRow <- which(startRow==1)
  ifelse(is.null(maxOut), maxOut <- length(tsOut-1), maxOut <- startRow + maxOut - 1)
  
  
# Now extract the information and put them into a dataframe

# The first row is a bit different, as it is contained in the header row, all others are then standardized
  temp1 <- strsplit(tsOut[startRow],"<td>")
  firstEntry <- paste(temp1[[1]][2],"<td>",temp1[[1]][3],sep="")
  res <- extractTSinfo(firstEntry)
  for(i in (startRow+1):maxOut){
    res <- rbind(res,extractTSinfo(tsOut[i]))
  }
  res
}

res <- targetScan(mirna="miR-9-5p", species="Cow", maxOut=5)