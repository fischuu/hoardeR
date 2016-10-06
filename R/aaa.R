# Collection of help functions

# Get the last element of a vector
lastElement <- function(x){
  x[length(x)]
}
#--------------------------------------------------------------------------------------------------------

# This function mixes two vectors alternating
mixVectors <- function(x,y){
  unlist(c(rbind(x, y)) )
}

# Merge two vectors into a data.frame by using the names of them as key
mergeVectors <- function(x,y){
  names.x <- names(x)
  names.y <- names(y)
  x.df <- data.frame(x=x)
  y.df <- data.frame(y=y)
  x.df$ID <- names.x
  y.df$ID <- names(y)
  merge(x.df,y.df,by="ID")
}
#--------------------------------------------------------------------------------------------------------

# Get a substring from the right end side of a string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#---------------------------------------------------------------------------------------------------------
# Remove the last n elements from a string
removeRight <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

#---------------------------------------------------------------------------------------------------------

# Order function, if one matrix should be brought into the order of another one
getSameOrder <- function(x,y){
 newPos <- numeric(length(y))
 for(i in 1:length(y)){
   newPos[i] <- which(y[i]==x)
 }
 newPos
}

#----------------------------------------------------------------------------------------------------

# Format the output to fit the NCBI ft server

adjustCHRLabel <- function(x){

  tmp <- suppressWarnings(as.numeric(x))
  if(!is.na(tmp)) x <- sprintf("%02d",tmp)
  x
}

#--------------------------------------------------------------------------------------------------------

# Get gene names from the V9 column of the annotation 

getGeneName <- function(x){
  if(grepl("gene_name", x)){
    geneName <-gsub(' \"', '', strsplit(strsplit(x,"gene_name")[[1]][2],'\";')[[1]][1])
  } else if(grepl("gene=", x)){
    geneName <- strsplit(strsplit(x,"gene=")[[1]][[2]],";")[[1]][[1]]
  } else if(grepl("gene_id", x)){
    geneName <- gsub(' \"', '', strsplit(strsplit(x,"gene_id")[[1]][2],'\";')[[1]][1])
  } else if(grepl("Name=",x)){
    geneName <- gsub(' \"', '', strsplit(strsplit(x,"Name=")[[1]][2],'\";')[[1]][1])
  }
  geneName
}

