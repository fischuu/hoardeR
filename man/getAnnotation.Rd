\name{getAnnotation}
\alias{getAnnotation}

\title{
 Downloading or Importing of Annotation Data
}

\description{
  This function downloads (if needed) the annotation file from a given species from Ensembl and loads it into the namespace.
}
\usage{
  getAnnotation(species=NULL, release=84, version=NULL, annotationFolder, type="gtf")
}

\arguments{
  \item{species}{The scientific name of the species (\code{String}).}
  \item{release}{The Ensembl release version.}
  \item{version}{The annotation version.}
  \item{annotationFolder}{The folder where the file is stored.}
  \item{type}{The file extension/encoding of the annotation file.}
}

\details{ This function downloads for a given \code{species} the annotation file, as provided from Ensembl. The main parameters basically define
the URL, where the file is located. The file is then downloaded into the folder, provided in \code{annotationFolder} and then imported to the namespace.

If a file has been downloaded previously, it will be loaded directly from that folder. In case the user wants to use an annotation that is not
located at Ensembl, the corresponding files can also be placed into the folder provided in \code{annotationFolder} and the function will load
it from there.
}

\value{
  A \code{data.table} with the annotation information.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
susScrofa <- getAnnotation(species = "Sus scrofa", annotationFolder="/home/user/annotation")
homoSapiens <- getAnnotation(species = "Homo sapiens", annotationFolder="/home/user/annotation")
}
}

\keyword{methods}
 
