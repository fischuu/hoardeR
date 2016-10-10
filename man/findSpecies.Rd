\name{findSpecies}
\alias{findSpecies}

\title{
 Search in the \code{species}` Object.
}

\description{
  This function output rows from the \code{species} object that contain a certain string.
}
\usage{
  findSpecies(string)
}

\arguments{
  \item{string}{Search string.}
}

\details{
  This function output rows from the \code{species} object that contain a certain string. It uses the \code{grepl} function to 
  identify the corresponding rows.
}
\value{
 A \code{data.frame}.
}

\author{
  Daniel Fischer
}

\seealso{
  \code{\link{species}}, \code{\link{grepl}}
}

\examples{
findSpecies("cattle")
}

\keyword{methods}
 
