\name{exportBed}
\alias{exportBed}

\title{
 Exporting a Bed File.
}

\description{
  This function exports a standard bed file.
}
\usage{
  exportBed(x, file=NULL, header=FALSE)
}

\arguments{
  \item{x}{A data.frame}
  \item{file}{Specifies the filename/path.}
  \item{header}{Logical, shall a header be written.}
}

\details{
 This function exports a data.frame to a standard bed file. If no file name is given, the variable name will be used instead.
}
\value{
 A bed file.
}

\author{
  Daniel Fischer
}


\examples{
\dontrun{
novelBed <- data.frame(Chr=c(11,18,3),
                         Start=c(72554673, 62550696, 18148822),
                         End=c(72555273, 62551296, 18149422),
                         Gene=c("LOC1", "LOC2", "LOC3"))

exportBed(novelBed, file="myLocs.bed")
exportBed(novelBed, file="myLocs_wHeader.bed", header=TRUE)
}
}

\keyword{methods}
 
