\name{tableSpecies}
\alias{tableSpecies}

\title{
Tables the species in xml file
}

\description{
Tables the species in xml file
}
\usage{
  importFA(file)
}

\arguments{
  \item{file}{Specifies the filename/path.}
}

\details{
 This function imports a standard fasta file. It assumes that label and sequence lines are alternating,
 meaning in the odd lines is the sequence name given, starting with > and in the even rows are the 
 corresponding sequences.
}
\value{
 An object of class \code{fa} containing the sequences. The names coorespond to the sequence names given in the fasta file.
}

\author{
  Daniel Fischer
}

\seealso{
\code{\link{print.fa}}, \code{\link{summary.fa}}
}

\examples{
\dontrun{
importFA(file="myFasta.fa")
}
}

\keyword{methods}
 
