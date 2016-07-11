# Get the last element of a vector
lastElement <- function(x){
  x[length(x)]
}
#--------------------------------------------------------------------------------------------------------

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

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#---------------------------------------------------------------------------------------------------------

