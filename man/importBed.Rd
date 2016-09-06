\name{importBed}
\alias{importBed}

\title{
 Importing a Bed File.
}

\description{
  This function imports a standard bed file.
}
\usage{
  importBed(file, header=FALSE)
}

\arguments{
  \item{file}{Specifies the filename/path.}
  \item{header}{Logical, is a header present.}
}

\details{
 This function imports a standard bed-file into a data.frame. It is basically a convenience wrapper around \code{read.table}.
}
\value{
 A \code{data.frame}.
}

\author{
  Daniel Fischer
}

\seealso{
  \code{\link{exportBed}}, \code{\link{read.table}}
}

\examples{
\dontrun{
novelBed <- data.frame(Chr=c(11,18,3),
                         Start=c(72554673, 62550696, 18148822),
                         End=c(72555273, 62551296, 18149422),
                         Gene=c("LOC1", "LOC2", "LOC3"))

exportBed(novelBed, file="myLocs.bed")

novelBed.imp <- importBed(file="myLocs.bed")
}
}

\keyword{methods}
 
