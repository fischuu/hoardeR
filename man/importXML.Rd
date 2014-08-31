\name{importXML}
\alias{importXML}

\title{
  Import a XML File
}

\description{
  This function imports a xml file produced from \code{blastSeq}.
}
\usage{
  importXML(seqNames, folder, which=c(1,2), idTH = 0.8, verbose=TRUE)
}

\arguments{
  \item{seqNames}{Sequence names.}
  \item{folder}{Fodler , where the xml files are stored.}
  \item{which}{Which of the provided sequence names should be imported.}
  \item{idTH}{Identity threshold, see details.}
  \item{verbose}{Logical, function give status messages.}
}

\details{
  This function imports a xml files produced from the \code{blastSeq} function.
  The \code{idTh} options sets the limit, what the minimum id threshold is until a hit will be
  taken into the result data frame.
}
\value{
  A data frame containing the results.
}

\author{
  Daniel Fischer
}

 
