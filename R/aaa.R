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

# Order function, if one matrix should be brought into the order of another one
getSameOrder <- function(x,y){
 newPos <- numeric(length(y))
 for(i in 1:length(y)){
   newPos[i] <- which(y[i]==x)
 }
 newPos
}

#----------------------------------------------------------------------------------------------------