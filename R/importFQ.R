# This function reads in a fasta file and prepares the vector from it
  importFQ <- function(file){
    
    .Deprecated("GenomicTools.fileHandler::importFQ", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
    
    res <- readLines(file)
  # Check if the Fasta file is alternating, one line label, the next line sequence
    sumAlternating <- sum(grepl("@",res)==c(TRUE,FALSE,FALSE,FALSE))
    if(sumAlternating!=length(res)) stop("Your Fastq file is malformed. Please ensure that name rows start with > and that names and sequence
                                         rows are alternating (plus the two additional rows for the quality).")
    seq <- res[seq(2,length(res),4)]
    qual <- res[seq(4,length(res),4)]
    names(seq) <- res[seq(1,length(res)-1,4)]
    names(qual) <- res[seq(1,length(res)-1,4)]
    
    structure(seq, qual=qual, class="fq")
} 
