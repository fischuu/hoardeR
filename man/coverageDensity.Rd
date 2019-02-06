\name{coverageDensity}
\alias{coverageDensity}

\title{
Calculation of the coverage density
}

\description{
Calculates the coverage density.
}
\usage{
  coverageDensity(folder, chr=c(1:22,"X","Y","MT"), chr.length=NULL,
                  posneg=FALSE, verbose=TRUE, use.sqrt=FALSE,
                  kernel.package="slideWindowSum",
                  step.size=50000, window.size=100000, bw=100)
}

\arguments{
  \item{folder}{folder with bam files}
  \item{chr}{Chromosome names to be plotted.}
  \item{chr.length}{Length of chromosome}
  \item{posneg}{Logical, plot pos and neg strand}
  \item{verbose}{Logical, verbose output} 
  \item{use.sqrt}{Logical, apply sqrt transformation}
  \item{kernel.package}{Class of kernel smoother} 
  \item{step.size}{Step size in bases}
  \item{window.size}{Window size in bases}
  \item{bw}{Bandwidth parameter}
}

\details{
 This function calculates the coverage of bam-files
 
}

\author{
  Daniel Fischer
}

\keyword{methods}
 
