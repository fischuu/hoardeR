\name{targetScan}
\alias{targetScan}

\title{
 Retrieving miRNA target information from targetscan.org
}

\description{
  This function requests from the webpage targetscan.org the stored information for mirnas.
}
\usage{
   targetScan(mirna=NULL, species=NULL, release="7.1", maxOut=NULL)
}

\arguments{
  \item{mirna}{The name of the mirna (\code{String}).}
  \item{species}{The species identifier, see details (\code{String}).}
  \item{release}{The release version of targetscan.org.}
  \item{maxOut}{The amount of target genes, default (\code{NULL}) is all.}
}

\details{
 This function sends a miRNA name to the targetscan.org webpage, retrieves the information and gives it back as a data.frame.
 Options for \code{species} are \code{"Human", "Mouse", "Rat", "Chimpanzee",} \code{"Rhesus", "Cow", "Dog", "Opossum", "Chicken", "Frog"}.
}

\value{
 A data.frame with the following columns

  \item{Ortholog}{The ortholog name of the target gene.}
  \item{geneName}{The long description of the target gene.}
  \item{consSites}{The total number of conserved sites.}
  \item{poorlySites}{The total number of poorly conserved sites.}
}

\references{
  \cite{V. Agarwal, G. Bell, J.Nam, et al. (2015): Predicting effective microRNA target sites in mammalian mRNAs. eLife, 4, pages 1-38, doi: 10.7554/eLife.05005}
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
targetScan(mirna="miR-9-5p", species="Cow", maxOut=5)
}
}

\keyword{methods}
 
