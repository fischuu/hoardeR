\name{getFastaFromBed}
\alias{getFastaFromBed}

\title{
Get fasta information based on locations in bed-format
}

\description{
For a given fasta and a bed file this function can extract the nucleotide sequences and stores them as fasta file.
}
\usage{
  getFastaFromBed(bed, species=NULL, release = "84", fastaFolder=NULL,
                  version=NULL, verbose=TRUE, export=NULL, filename=NULL)
}

\arguments{
  \item{bed}{The location in bed format, see details.}
  \item{species}{Define the species.}
  \item{release}{Ensembl release version.}
  \item{fastaFolder}{Location of the fasta files.}
  \item{version}{Version ID of the reference genome.}
  \item{verbose}{Logical, should informative status updates be given.}
  \item{export}{Foldername.}
  \item{filename}{Filename to store the FA object.}  
}

\details{
 Function expects as an input a \code{data.frame} in bed format. This means, the first column should contain the chromosome, the second
 the start-coordinates, the third the end-coordinates. The forth column contains the ID of the loci. 
 
 If a standard species, release and version is provided and the data is available from the Ensembl webpage, the functions automatically
 downloads the required files, takes the loci and extracts then the nucleotide sequences from it. 
 
 The export function, specifies then a folder to where the fasta file should be stored. If no filename is provided, the filename is then
 the object name passed to the \code{bed} function.
}
\value{
An \code{fa} object containing the nucleotide sequences in fasta format.
}

\author{
  Daniel Fischer
}


\examples{
\dontrun{

myBed <- data.frame(chr=c(1,2),
                    start=c(235265,12356742),
                    end=c(435265,12386742),
                    gene=c("LOC1", "LOC2"))

myFA <- getFastaFromBed(myBed, species="Homo sapiens", fastaFolder="/home/user/fasta/", export=TRUE)
}
}

\keyword{methods}
 
