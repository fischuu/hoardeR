\name{getFastaFromBed}
\alias{getFastaFromBed}

\title{
Get fasta information based on locations in bed-format
}

\description{
For a given fasta and a bed file this function can extract the nucleotide sequences and stores them as fasta file.
}
\usage{
  getFastaFromBed(bed, species=NULL, assembly = NULL, fastaFolder=NULL,
                  verbose=TRUE, export=NULL, fileName=NULL)
}

\arguments{
  \item{bed}{The location in bed format, see details.}
  \item{species}{Define the species.}
  \item{assembly}{Assembly identifier.}
  \item{fastaFolder}{Location of the fasta files.}
  \item{verbose}{Logical, should informative status updates be given.}
  \item{export}{Foldername.}
  \item{fileName}{Filename to store the FA object.}  
}

\details{
 Function expects as an input a \code{data.frame} in bed format. This means, the first column should contain the chromosome, the second
 the start-coordinates, the third the end-coordinates. The forth column contains the ID of the loci. 
 
 If a standard species is used (as defined in the \code{species} data frame), the function automatically downloads the required files
 from NCBI, takes the loci and extracts then the nucleotide sequences from it. If the corresponding assemly is not available from NCBI
 an own fasta file can be provided. For that the fa-file needs to be in the fastaFolder and follow the same naming system as the NCBI 
 files are labelled. In that case, the function suggests the correct filename for an unknown assembly.
 
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
 
