\name{print.fa}
\docType{methods}
\alias{print,fa-method}
\alias{print.fa}

\title{Print an fa Object}
\description{
  Prints an \code{fa} object.
}
\usage{
 \method{print}{fa}(x, n=2, seq.out=50, ...)
}

\arguments{
  \item{x}{Object of class \code{fa}.}
  \item{n}{Amount of elements to be displayed, numeric.}
  \item{seq.out}{Length of each element to be displayed, numeric..}  
  \item{...}{Additional parameters.}
}

\details{
The print function displays an \code{fa} object. By default just the first two elements with their first 50 bases are
displayed. To display the full sequence, set \code{seq.out=NULL}.
}

\author{Daniel Fischer}

\keyword{methods}
\keyword{print}
