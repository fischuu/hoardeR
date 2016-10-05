findSpecies <- function(string){
# First find all the cells that contain the string
  cellsOI <- apply(species, 2, function(x) grepl(string, x))
# Now determine all the rows that of interest
  rowsOI <- apply(cellsOI,1,sum)
# Display the rows of interest
  species[rowsOI>0,]
}