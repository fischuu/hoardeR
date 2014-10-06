\name{subPhased}
\alias{subPhased}

\title{
 Rewrite the Phased File from a Beagle Output
}

\description{
  This function takes a phased Beagle output and rewrites the output.
}
\usage{
  subPhased(file=NULL, vmmk = NULL, out=NULL, chunkSize=100000, verbose=TRUE,
            removeInsertions=TRUE)
}

\arguments{
  \item{file}{Location of the original Beagle file (\code{String}).}
  \item{vmmk}{Location of the Variant Map Master key (\code{String}).}
  \item{out}{Name and location of the output file (\code{String}).}
  \item{chunkSize}{For large Beagle files, the chunk size.}
  \item{verbose}{The function gives feedback.}
  \item{removeInsertions}{All Indels will be removed.}
}

\details{
 This function takes a Beagle phased file and rewrites the alleles from numerical to character, based
 on the information provided in a variant map master key. For larger files the function can process
 the rewriting in chunks in order to save memory.
}

\value{
  A rewritten beagle phased file.
}

\author{
  Daniel Fischer
}


\keyword{methods}
 
