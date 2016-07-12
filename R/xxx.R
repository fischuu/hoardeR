# Make the species data available without loading it.

.onLoad <- function(libname, pkgname) {
  utils::data("species", package="hoardeR", envir=parent.env(environment()))
}