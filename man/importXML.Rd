\name{importXML}
\alias{importXML}

\title{
  Import an XML File produced by the blastSeq function
}

\description{
  This function imports an xml file produced from \code{blastSeq}.
}
\usage{
  importXML(folder, seqNames=NULL, which=NULL, idTH = 0.8, verbose=TRUE)
}

\arguments{
  \item{folder}{Folder where the xml files are stored.}
  \item{seqNames}{Sequence names.}
  \item{which}{Which of the provided sequence names should be imported.}
  \item{idTH}{Identity threshold, see details.}
  \item{verbose}{Logical, function give status messages.}
}

\details{
  This function imports xml files produced from the \code{blastSeq} function.
  The \code{idTh} options sets the limit, what the minimum id threshold is until a hit will be
  taken into the result data frame.
  
  If \code{seqNames} is \code{NULL} (default) and code{which} is \code{NULL} (default) all sequences are imported. The list can be limited, either by
  specifying the sequence name (the file name of it), or its position in the order list of names, using \code{which}.
}
\value{
  A data frame containing the results.
}

\author{
  Daniel Fischer
}

 
