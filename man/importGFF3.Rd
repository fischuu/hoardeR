\name{importGFF3}
\alias{importGFF3}

\title{
  Import a GFF3 File
}

\description{
  This function imports a gff3 file.
}
\usage{
  importGFF3(gff)
}

\arguments{
  \item{gff}{File name of the gff3 file}
}

\details{
  This function imports a gff file and splits the last column which is usually tricky to handle as the order of the variables is not always
  the same.
}
\value{
  A data frame containing the columns of the gtf file, including the splitted last column.
}

\author{
  Daniel Fischer
}

 
