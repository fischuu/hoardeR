\name{importGTF}
\alias{importGTF}

\title{
  Import a GTF File
}

\description{
  This function imports a gtf file.
}
\usage{
  importGTF(gtf, skip = 0, nrow = -1)
}

\arguments{
  \item{gtf}{File name of the gtf file.}
  \item{skip}{Rows to skip from the top.}
  \item{nrow}{Total amount of rows read.}
}

\details{
  This function imports a gtf file and splits the column 9 which is usually tricky to handle as the order of the variables is not always
  the same.
}
\value{
  A data frame containing the columns of the gtf file, including the splitted last column.
}

\author{
  Daniel Fischer
}

 
