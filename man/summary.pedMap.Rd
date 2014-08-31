\name{summary.pedMap}
\docType{methods}
\alias{summary.pedMap}
\alias{summary,pedMap-method}

\title{Summarize an pedMap Object}
\description{
  Summarizes an \code{pedMap} object in an informative way.
}
\usage{
 \method{summary}{pedMap}(object, ...)
}

\arguments{
  \item{object}{Object of class \code{pedMap}.}
  \item{...}{Additional parameters.}
}

\details{
 Summary for a \code{pedMap} object, providing the dimensions of the \code{map}, the \code{fam} and the
 \code{geno} matrix as well as the total amount of Allele A/A, A/B, B/B as well as the amount of missing
 data and the monomorhic locations.
}

\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
