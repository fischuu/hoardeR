getFastaFromBed <- function(bed, species=NULL, release = "84", fastaFolder=NULL, version=NULL, verbose=TRUE, export=NULL, fileName=NULL){
  
  if(is.null(fileName)) fileName <- deparse(substitute(bed))
  exportFileName <- fileName
  
# Input checks
  if(is.null(bed)) stop("No bed data given!")
  if(is.null(nrow(bed))) bed <- t(as.matrix(bed))
  bed <- data.frame(Chr=bed[,1],
                    Start=bed[,2],
                    End=bed[,3],
                    Gene=bed[,4])
  if(is.null(species)) stop("No species given!")    
  if(is.null(fastaFolder)) stop("No directory with fasta files given!")   
  if(is.null(version)){
    temp <- hoardeR::species
    speciesVersion <- temp$Ensembl.Assembly[temp$Scientific.name==species]    
  }   
  if(verbose) cat("Using species version:", speciesVersion,"\n")
  
# First the required original fasta file
  dir.create(fastaFolder, showWarnings=FALSE)
  species.int <- gsub(" ", "_",tolower(species))
  ensemblURL <- paste("ftp://ftp.ensembl.org/pub/release-",release,"/fasta/",species.int,"/dna/",sep="")

# Now loop through all bed rows, download the required fasta files and extract the fasta information from there
  novelFA <- list()
  for(bedRun in 1:nrow(bed)){
    fileName <-  paste(cap(species.int),".",speciesVersion,".dna.chromosome.",bed$Chr[bedRun],".fa.gz",sep="")
    .file = file.path(fastaFolder, fileName)
    if(!file.exists(.file)) download.file(paste(ensemblURL,fileName,sep=""), .file)
    if(verbose) cat("Read original fasta file:\n   ", .file ,"\n")
    seqSpecies <- read.fasta(.file,seqtype="DNA")
    novelFA[bedRun] <- paste(seqSpecies[[1]][bed$Start[bedRun]:bed$End[bedRun]],collapse="")
  }
  names(novelFA) <- paste(">",bed[,1],":",bed[,2],"-",bed[,3],sep="")
  if(!is.null(export)){
    fileConn<-file(file.path(fastaFolder,paste(exportFileName,".fa",sep="")))
    writeLines( mixVectors(names(novelFA), toupper(unlist(novelFA))), fileConn)
    close(fileConn)
  }
  class(novelFA) <- "fa"
  novelFA
  
}