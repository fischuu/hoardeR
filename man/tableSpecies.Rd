\name{tableSpecies}
\alias{tableSpecies}

\title{
Tables the species in xml file
}

\description{
Tables the species in xml file
}
\usage{
  tableSpecies(xml, species=NULL, type="chr", minOutput=TRUE, exclude="",
               locations=FALSE)
}

\arguments{
  \item{xml}{The xml file.}
  \item{species}{Restrict species to a certain set.}
  \item{type}{Filter option.}
  \item{minOutput}{Logical, should the output be minimal.}
  \item{exclude}{Names of species to exclude.}
  \item{locations}{Logical, shall the hit locations be given as well.}
  
}

\details{
 Function provides a table of identified species. This table can e.g. be put into the \code{barplot} function to visualize the findings.
 
 Further, if the option \code{locations} is set to \code{TRUE} the function not only tables the species, but also the individual locations
 of the hits. This output is required for the further steps. Hence, this function plays a important role in the identification pipeline.
}
\value{
A \code{table} with the species from the XML file
}

\author{
  Daniel Fischer
}


\examples{
\dontrun{
tableSpecies(xmls)
pigHits <- tableSpecies(xmls, species="Sus scrofa", locations = TRUE)
}
}

\keyword{methods}
 
