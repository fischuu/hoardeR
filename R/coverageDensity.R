# Some help functions for the bam handling
# function for collapsing the list of lists into a single list
# as per the Rsamtools vignette
  .unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
# function for checking negative strand  
  check_neg <- function(x){
    if (intToBits(x)[5] == 1){ 
      return(T)
    } else {
      return(F)
    }
  }
  
# function for checking positive strand
  check_pos <- function(x){
    if (intToBits(x)[3] == 1){
      return(F)
    } else if (intToBits(x)[5] != 1){
      return(T)
    } else {
      return(F)
    }
  }

# And the main function
  coverageDensity <- function(folder, chr=c(1:22,"X","Y","MT"), posneg=FALSE, verbose=TRUE, bw="SJ"){
  
  # Get the filenames of the bams
    bamNames <- list.files(folder, pattern="*.bam$")

  # Create the bam list object
    bam_df_list <- list()

  # Now loop through the bam, import them and store them as dataFrame in R  
    for(i in 1:length(bamNames)){
    # Import the bam
      tmpBam <- scanBam(file.path(folder,bamNames[i]))  
    # store names of BAM fields
      bam_field <- names(tmpBam[[1]])
    # go through each BAM field and unlist
      bam_list <- lapply(bam_field, function(y) .unlist(lapply(tmpBam, "[[", y)))
    # store as data frame
      bam_df <- do.call("DataFrame", bam_list)
      names(bam_df) <- bam_field
    # Write the dataFrame into the list 
      bam_df_list[[i]] <- bam_df
    # OUtput message
      if(verbose) message("Imported sample",bamNames[i])
    }
    
    names(bam_df_list) <- bamNames
    chrIndex <- 1
    
  # CASE: Consider postive and negative strand reads  

    if(posneg){
    # Prepare the list objects for the densities
      chr_neg_density <- list()
      chr_pos_density <- list()

      negOut <- list()    
      posOut <- list()
    # Loop through all the Chromosomes   
      for(chrRun in chr){
      # Now loop again through all samples and calculate the densities
        for(i in 1:length(bam_df_list)){
        # store the mapped positions on the plus and minus strands
          chr_neg <- bam_df_list[[i]][bam_df_list[[i]]$rname == chrRun &
                                        apply(as.data.frame(bam_df_list[[i]]$flag), 1, check_neg),
                                        'pos'
                                        ]
        
          chr_pos <- bam_df_list[[i]][bam_df_list[[i]]$rname == chrRun &
                                      apply(as.data.frame(bam_df_list[[i]]$flag), 1, check_pos),
                                      'pos'
                                      ]
        # calculate the densities
          chr_neg_density[[i]] <- density(chr_neg, bw=bw)
          chr_pos_density[[i]] <- density(chr_pos, bw=bw)
        
        # display the negative strand with negative values
          chr_neg_density[[i]]$y <- chr_neg_density[[i]]$y * -1
        
        # Output   
          message("Processed sample",bamNames[i],"in Chromosome",chrRun)
        }
        negOut[[chrIndex]] <- chr_neg_density
        posOut[[chrIndex]] <- chr_pos_density     
        chrIndex <- chrIndex + 1
      }
      
      out <- list(pos=posOut, 
                  neg=negOut)
      
    # CASE: Consider merged reads, regardless which strand
    } else {
      
      chr_density <- list()
      out <- list()
      
      for(chrRun in chr){
        # Now loop again through all samples and calculate the densities
        for(i in 1:length(bam_df_list)){
          # store the mapped positions on the plus and minus strands
          chr <- bam_df_list[[i]][which(bam_df_list[[i]]$rname == chrRun),'pos']
          
          # calculate the density
          chr_density[[i]] <- density(chr, bw=bw)
          
          # Output   
          message("Processed sample",bamNames[i],"in Chromosome",chrRun)
        }
        out[[chrIndex]] <- chr_density
        chrIndex <- chrIndex + 1
      }
      
    }
    
    class(out) <- "coverageDensity"
    out
  }