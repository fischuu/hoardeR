\docType{data}
\name{species}
\alias{species}
\title{Available species at NCBI}
\format{A data frame with 348 species.}
\source{
As downloaded on 05.10.2016 from

ftp://ftp.ncbi.nlm.nih.gov/genomes/
}
\description{
This is a list of all organisms/species that are provided by NCBI and hence could end up in the Blast run. Further, it defines the default versions of
the assuemblies that will be downloaded if no further version is specified in \code{plotHit}, \code{getAnnotation} or \code{getFastaFromBed}.
}
\examples{

data(species)
summary(species) 

}
\keyword{datasets}

