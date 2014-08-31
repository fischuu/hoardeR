\name{getEnsgInfo}
\alias{getEnsgInfo}

\title{
 Retrieve Gene Information From The NCBI Database.
}

\description{
  This function retrieves for a given Ensemble Number the corresponding information from the NCBI database.
}
\usage{
  getEnsgInfo(ensg)
}

\arguments{
  \item{ensg}{Ensemble ID (\code{String}).}
}

\details{
 This function retrieves for a given Ensemble Number the corresponding information from the NCBI database. The 
 object \code{ensg} can also be a vector of Ensemble IDs.
}
\value{
  A matrix with information.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
ensg <- c("ENSG00000174482", "ENSG00000113494")
getEnsgInfo(ensg)
}
}

\keyword{methods}
 