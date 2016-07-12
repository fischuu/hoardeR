\docType{data}
\name{species}
\alias{species}
\title{Available species at Ensembl}
\format{A data frame with XXX species.}
\source{
As downloaded on 20.6.2016 from

http://www.ensembl.org/info/about/species.html
}
\description{
This is a list of all organisms/species that are provided by Ensembl and hence could end up in the Blast run. Further, it defines the default versions of
the assuemblies that will be downloaded if no further version is specified in \code{plotHit}, \code{getAnnotation} or \code{getFastaFromBed}.
}
\examples{

data(species)
summary(species) 

}
\keyword{datasets}

