\name{intersectXMLAnnot}
\alias{intersectXMLAnnot}

\title{
Intersect XML object with annotation object
}

\description{
For a annotation object this function intersects the loci of it with the output of the tableSpecies function.
}
\usage{
  intersectXMLAnnot(tabSpecies, annot, level="gene", flanking=NULL)
}

\arguments{
  \item{tabSpecies}{The table with locations from \code{tableSpecies}.}
  \item{annot}{The annotation object.}
  \item{level}{The level of intersection.}
  \item{flanking}{Allowed flanking space for intersection.}
}

\details{
 Function expects as an input table from \code{tableSpecies} with the option \code{locations=TRUE}. Further, it needs an annotation object,
 as provided by the \code{getAnnotation} function. With that it intersects then the loci on the level as specified in \code{level}. Currently
 only \code{"gene"} is supported. 
 
 The \code{flanking} option allows for flanking space up- and down-stream of the genes. This is especially then useful if the novel gene
 candidates are in the extension of known genes (e.g. responsible for regulation or if they are novel exons.)
}
\value{
A table with intersection loci.
}

\author{
  Daniel Fischer
}


\examples{
\dontrun{

pigHits <- tableSpecies(xmls, species="Sus scrofa", locations = TRUE)
ssannot <- getAnnotation(species = "Sus scrofa", annotationFolder="/home/user/annotation")
pigInter <- list()
for(i in 1:nrow(pigHits)){
   pigInter[[i]] <- intersectXMLAnnot(pigHits[i,], ssannot)
}
}

}

\keyword{methods}
 
