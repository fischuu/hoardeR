# This function is not exported! It is used internally to populate the species table
# It requires in addion the 

# library(RCurl)  #getURL

# Help function to extract information from the README vector
getInfo <- function(x, entry){
  takeThis <- which(grepl(entry, x))
  out <- "NA"
  if(length(takeThis)>0)  out <- strsplit(x[takeThis], paste(entry,":\t",sep=""))[[1]][2]
  out
}

createSpecies <- function(baseDir="ftp://ftp.ncbi.nlm.nih.gov/genomes/", verbose=TRUE){
  
  curl <- getCurlHandle()
  
  out <- data.frame(Organism.Name = character(),
                    Organism.Common.Name = character(),
                    Taxid = character(),
                    Assembly.Name = character(),
                    Assembly.Accession = character(),
                    Assembly.Submitter = character(),
                    Assembly.Data = character(),
                    Assembly.Type = character(),
                    Annotation.Release.Name = character(),
                    Annotation.Release.Date = character(),
                    Annotation.Report = character(),
                    NCBI.Url = character(),
                    stringsAsFactors=FALSE)
  
# First get all folders
  filenames <- getURL(baseDir, ftp.use.epsv = TRUE, dirlistonly = TRUE)
  filenames <- sort(strsplit(filenames,"\n")[[1]])
  
# Remove already known, ununser folders
  filenames <- filenames[-which(is.element(filenames,c("all",
                                                       "archive",
                                                       "genbank",
                                                       "refseq",
                                                       "ASSEMBLY_REPORTS",
                                                       "GENOME_REPORTS",
                                                       "HUMAN_MICROBIOM",
                                                       "INFLUENZA",
                                                       "TARGET",
                                                       "TOOLS",
                                                       "Viruses")))]
  
  filenames <- filenames[which(!grepl("README",filenames))]
  
# Now check for each entry, if CHR_X, GFF, and README_CURRENT_RELEASE
  for(i in 1:length(filenames)){
  # Get the subfolders for the filename
    currentURL <- paste(baseDir,filenames[i],"/",sep="")
    tmp <- strsplit(getURL(currentURL, ftp.use.epsv = FALSE, dirlistonly = TRUE),"\n")[[1]]
    
  # Check if it contains the required README file
    if(sum(grepl("README_CURRENT_RELEASE", tmp))==1){
    # Get the README file
      tmpRM <- strsplit(getURL(paste(baseDir,filenames[i],"/README_CURRENT_RELEASE",sep=""), curl=curl),"\n")[[1]]
    # Extract the required information
      tmpEntry <- data.frame(Organism.Name = getInfo(tmpRM, "ORGANISM NAME"),
                             Organism.Common.Name = getInfo(tmpRM, "ORGANISM COMMON NAME"),
                             Taxid = getInfo(tmpRM, "TAXID"),
                             Assembly.Name = getInfo(tmpRM, "ASSEMBLY NAME"),
                             Assembly.Accession = getInfo(tmpRM, "ASSEMBLY ACCESSION"),
                             Assembly.Submitter = getInfo(tmpRM, "ASSEMBLY SUBMITTER"),
                             Assembly.Data = getInfo(tmpRM, "ASSEMBLY DATE"),
                             Assembly.Type = getInfo(tmpRM, "ASSEMBLY TYPE"),
                             Annotation.Release.Name = getInfo(tmpRM, "ANNOTATION RELEASE NAME"),
                             Annotation.Release.Date = getInfo(tmpRM, "ANNOTATION RELEASE DATE"),
                             Annotation.Report = getInfo(tmpRM, "ANNOTATION REPORT"),
                             NCBI.Url = currentURL,
                             stringsAsFactors=FALSE)
      out <- rbind(out,tmpEntry)
    }
    if(verbose) cat("Species", filenames[i],i,"/",length(filenames),"processed \n")
  }
  out
}