getFastaFromBed <- function(bed, species=NULL, assembly = NULL, fastaFolder=NULL, verbose=TRUE, export=NULL, fileName=NULL){
  
# Make the species object accessible here
  species.df <- hoardeR::species
  
# Check the species
  if(nrow(findSpecies(species))<1) error("Provided species could not be found in the 'species' table")
  
# If no fileName is given and the fasta should be exported, use the variable name of the bed object as fileName
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
  if(is.null(fastaFolder)){
    message("No directory with fasta files given! Use the working directory: \n", getwd())
    fastaFolder <- getwd()
  } else {
      if(substrRight(fastaFolder,1)=="/") fastaFolder <- removeRight(fastaFolder,1)
    } 
  if(is.null(assembly)){
    assemblyWasNULL <- TRUE
    speciesAssembly <- species.df$Assembly.Name[species.df$Organism.Name==species][1]    
  } else {
    assemblyWasNULL <- FALSE
    speciesAssembly <- assembly
  }  
  if(verbose) cat("Using species assembly:", speciesAssembly,"\n")
  
# First the required original fasta file
  dir.create(fastaFolder, showWarnings=FALSE)
  species.int <- gsub(" ", "_",tolower(species))
  #ensemblURL <- paste("ftp://ftp.ensembl.org/pub/release-",release,"/fasta/",species.int,"/dna/",sep="")
  NCBI.URL <- species.df$NCBI.Url[species.df$Organism.Name==species][1]
  
# Now loop through all bed rows, download the required fasta files and extract the fasta information from there
  novelFA <- list()
  for(bedRun in 1:nrow(bed)){
  # This is to the CHR_ folders:
  #  filePath <- paste(NCBI.URL,"CHR_",adjustCHRLabel(bed$Chr[bedRun]),"/", sep="") 
    filePath <- paste(NCBI.URL,"Assembled_chromosomes/seq/", sep="")
    CHRfile <- getURL(filePath, ftp.use.epsv = TRUE, dirlistonly = TRUE)
    CHRfile <- strsplit(CHRfile,"\n")[[1]]
    CHRfile <- CHRfile[grepl(paste(species.df$Assembly.Name[species.df$Organism.Name==species][1], "_chr",bed$Chr[bedRun],".fa.gz",sep=""), CHRfile)]
    
    if(!assemblyWasNULL){
      CHRfileRequested <- gsub(species.df$Assembly.Name[species.df$Organism.Name==species][1], speciesAssembly, CHRfile) 
      CHRfile <- CHRfileRequested
    }
    
    .file = file.path(fastaFolder, CHRfile)
    if(!file.exists(.file)){
      if(verbose) cat("Local file not found! Try to download fasta file: ", CHRfile ,"\n")
      download.file(paste(filePath,CHRfile,sep=""), .file, quiet = TRUE)
    }
    if(file.size(.file)<1){
      cat("File could not be downloaded. If you want to use the assembly", speciesAssembly, "for", species,"please provide the file: ", .file,"\n")       
    } else {
      if(verbose) cat("Read fasta file: ", .file ,"\n")
        tmp <- scanFa(.file)
        novelFA[bedRun] <- as.character(subseq(tmp,bed$Start[bedRun],bed$End[bedRun]))
    }
  }
  if(length(novelFA)>0) names(novelFA) <- paste(">",bed[,1],":",bed[,2],"-",bed[,3],sep="")
  if(!is.null(export)){
    fileConn<-file(file.path(fastaFolder,paste(exportFileName,".fa",sep="")))
    writeLines( mixVectors(names(novelFA), toupper(unlist(novelFA))), fileConn)
    close(fileConn)
  }
  class(novelFA) <- "fa"
  novelFA
  
}