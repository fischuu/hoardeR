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
  
# SlideWindow
  slideWindowSum <- function(x, from, to, step.size, window.size){
    out <- c()
    
    tmpStart <- from
    tmpEnd <- from+window.size
    
    index <- 1
    
    while(tmpEnd < to){
      out[index] <- sum(x>=tmpStart & x<=tmpEnd)
      
      tmpStart <- tmpStart + step.size
      tmpEnd <- tmpEnd + step.size
      index <- index + 1
    }
    
    out
  }
  
# And the main function
  coverageDensity <- function(folder, chr=c(1:22,"X","Y","MT"), chr.length=NULL, posneg=FALSE, verbose=TRUE, use.sqrt=FALSE, kernel.package="slideWindowSum", step.size=50000, window.size=100000){
  
  # Input checks
    if(is.null(chr.length)) warning("No length information for chromosomes provided, density estimation might go wrong without it!")
    kernel.package <- match.arg(kernel.package,c("slideWindowSum","stats","KernSmooth"))
    
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
      
          # calculate the density
          if(kernel.package=="stats"){
            tmp_neg_density <- density(chr_neg, bw=bw)
            tmp_pos_density <- density(chr_pos, bw=bw)
            chr_neg_density[[i]] <-  list(x=tmp_neg_density$x,
                                          y=tmp_neg_density$y)
            chr_pos_density[[i]] <-  list(x=tmp_pos_density$x,
                                          y=tmp_pos_density$y)
          } else if(kernel.package=="KernSmooth"){
            chr_neg_density[[i]] <- bkde(chr_neg, kernel = "normal", gridsize = (chr.length[chrIndex]+(step.size-chr.length[chrIndex]%%step.size))/step.size+1, range.x=c(0,chr.length[chrIndex]))            
            chr_pos_density[[i]] <- bkde(chr_pos, kernel = "normal", gridsize = (chr.length[chrIndex]+(step.size-chr.length[chrIndex]%%step.size))/step.size+1, range.x=c(0,chr.length[chrIndex]))
          
          } else if(kernel.package=="slideWindowSum")
          
        # Transfor the values, if requested
          if(use.sqrt){
            chr_neg_density[[i]]$y <- sqrt(chr_neg_density[[i]]$y)
            chr_pos_density[[i]]$y <- sqrt(chr_pos_density[[i]]$y)
          } 
          
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
          if(kernel.package=="stats"){
            tmp_density <- density(chr, bw=bw, n=5000, from=0, to=250000000)
            chr_density[[i]] <-  list(x=tmp_density$x,
                                      y=tmp_density$y)
          } else if(kernel.package=="KernSmooth"){
            chr_density[[i]] <- bkde(chr, kernel = "normal", gridsize = (chr.length[chrIndex]+(step.size-chr.length[chrIndex]%%step.size))/step.size+1, range.x=c(0,chr.length[chrIndex]+(step.size-chr.length[chrIndex]%%step.size)), bandwidth=floor(seqlengths(Hsapiens)[1:25]/bw.para) )            
          } else if(kernel.package=="slideWindowSum"){
            tmp_density <- slideWindowSum(chr, from=0, to=chr.length[chrIndex], step.size=50000, window.size = 100000)
            chr_density[[i]] <- list(x=1:length(tmp_density),
                                     y=tmp_density)
          }

          # Transfor the values, if requested
          if(use.sqrt){
            chr_density[[i]]$y <- sqrt(chr_density[[i]]$y)
          } 
          
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