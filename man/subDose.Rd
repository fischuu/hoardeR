\name{subDose}
\alias{subDose}

\title{
 Rewrite the Dose File from a Beagle Output
}

\description{
  This function takes a Dose Beagle output and rewrites the output.
}
\usage{
  subDose(file=NULL, vmmk=NULL, out=NULL, removeInsertions=TRUE, verbose=TRUE)
}

\arguments{
  \item{file}{Location of the original Beagle file (\code{String}).}
  \item{vmmk}{Location of the Variant Map Master key (\code{String}).}
  \item{out}{Name and location of the output file (\code{String}).}
  \item{verbose}{The function gives feedback.}
  \item{removeInsertions}{All Indels will be removed..}
}

\details{
 This function takes a Beagle Dose file and rewrites the alleles from numerical to character, based
 on the information provided in a variant map master key. 
}

\value{
  A rewritten beagle phased file.
}

\author{
  Daniel Fischer
}


\keyword{methods}
 
