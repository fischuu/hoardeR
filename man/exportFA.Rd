\name{exportFA}
\alias{exportFA}

\title{
 Exporting a Fasta File.
}

\description{
  This function exports a standard fasta file.
}
\usage{
  exportFA(fa, file)
}

\arguments{
  \item{fa}{An \code{fa} object.}
  \item{file}{Specifies the filename/path.}
}

\details{
 This function exports a \code{fa} object to a standard fasta file. 
}
\value{
 A fasta file.
}

\author{
  Daniel Fischer
}

\seealso{
\code{\link{print.fa}}, \code{\link{summary.fa}}, \code{\link{importFA}}
}

\examples{
\dontrun{
novelBed <- data.frame(Chr=c(11,18,3),
                         Start=c(72554673, 62550696, 18148822),
                         End=c(72555273, 62551296, 18149422),
                         Gene=c("LOC1", "LOC2", "LOC3"))

myFasta <- getFastaFromBed(novelBed, species="Bos taurus", 
+                          fastaFolder="/home/daniel/fasta/")

exportFA(myFasta, file="myFasta.fa")
}
}

\keyword{methods}
 
