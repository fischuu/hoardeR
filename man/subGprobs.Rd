\name{subGprobs}
\alias{subGprobs}

\title{
 Rewrite the Gprobs File from a Beagle Output
}

\description{
  This function takes a Gprobs Beagle output and rewrites the output.
}
\usage{
  subGprobs(file=NULL, vmmk=NULL, out=NULL, chunkSize=100000, removeInsertions=TRUE,
             verbose = TRUE)
}

\arguments{
  \item{file}{Location of the original Beagle file (\code{String}).}
  \item{vmmk}{Location of the Variant Map Master key (\code{String}).}
  \item{out}{Name and location of the output file (\code{String}).}
  \item{chunkSize}{For large Beagle files, the chunk size.}
  \item{removeInsertions}{All Indels will be removed.}
  \item{verbose}{The function gives feedback.}
}

\details{
 This function takes a Beagle Gprobs file and rewrites the alleles from numerical to character, based
 on the information provided in a variant map master key. 
 For larger files the function can process
 the rewriting in chunks in order to save memory.
}

\value{
  A rewritten beagle Gprobs file.
}

\author{
  Daniel Fischer
}


\keyword{methods}
 
