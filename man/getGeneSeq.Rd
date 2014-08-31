\name{getGeneSeq}
\alias{getGeneSeq}

\title{
 Extracting a gene sequence from NCBI database.
}

\description{
  This function retrieves a gene sequence from the NCBI database.
}
\usage{
  getGeneSeq(chr, start, end, organism)
}

\arguments{
  \item{chr}{Chromosome number, numeric/string}
  \item{start}{Start position, numeric}
  \item{end}{End position, numeric}
  \item{organism}{Name of the organism, string}
}

\details{
 Extracting a gene sequence from NCBI database.
}
\value{
  A string that contains the genomic sequence.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
# Extracting for Sus Scrofa, build version 3:
getGeneSeq(1,2134,14532,"susScr3")
}
}

\keyword{methods}
 
