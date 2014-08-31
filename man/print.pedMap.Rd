\name{print.pedMap}
\docType{methods}
\alias{print,pedMap-method}
\alias{print.pedMap}

\title{Print an pedMap Object}
\description{
  Prints an \code{pedMap} object.
}
\usage{
 \method{print}{pedMap}(x, nrow=5, ncol=10, ...)
}

\arguments{
  \item{x}{Object of class \code{pedMap}.}
  \item{nrow}{Amount of rows to be displayed, numeric.}
  \item{ncol}{Amount of cols to be displayed, numeric.}
  \item{...}{Additional parameters.}
}

\details{
The print function displays an \code{pedMap} object. By default just the first 5 rows of each
list item and the first 10 columns of the \code{geno} matrix are displayed.
}

\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
