\name{getAnnotation}
\alias{getAnnotation}

\title{
 Downloading or Importing of Annotation Data
}

\description{
  This function downloads (if needed) the annotation file from a given species from NCBI and loads it into the namespace.
}
\usage{
  getAnnotation(species=NULL, assembly=NULL, annotationFolder=NULL, 
                type="gff3", verbose)
}

\arguments{
  \item{species}{The scientific name of the species (\code{String}).}
  \item{assembly}{The NCBI assembly version.}
  \item{annotationFolder}{The folder where the file will be stored.}
  \item{type}{The file extension/format of the annotation file.}
  \item{verbose}{Logical, if function gives feedback.}
}

\details{ This function downloads for a given \code{species} the annotation file, as provided from NCBI. The main parameters basically define the URL, where the file is located. The file is then downloaded into the folder, provided in \code{annotationFolder} and then imported to the namespace.

If a file has been downloaded previously, it will be loaded directly from that folder. In case the user wants to use an annotation that is not provided by NCBI, the corresponding files can also be placed into the same folder, following the naming scheme as suggested from the function and the function will load it from there.
}

\value{
  A \code{data.table} with the annotation information.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
susScrofa <- getAnnotation(species = "Sus scrofa", 
                           annotationFolder="/home/user/annotation")
                           
homoSapiens <- getAnnotation(species = "Homo sapiens", 
                             annotationFolder="/home/user/annotation")
}
}

\keyword{methods}
 
