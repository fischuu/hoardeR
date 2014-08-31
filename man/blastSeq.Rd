\name{blastSeq}
\alias{blastSeq}

\title{
 Sending Genomic Sequences to NCBI Blast service
}

\description{
  This function sends genomic sequences to the NCBI Blast service.
}
\usage{
  blastSeq(seq, n_blast=20, delay_req=3, delay_rid=60, email=NULL, 
           xmlFolder=NULL, keepInMemory=TRUE, database="chromosome",
           verbose=TRUE)
}

\arguments{
  \item{seq}{The fasta sequence that should be blasted (\code{String}).}
  \item{n_blast}{Amount of parallel blast requests, in case \code{seq} is a vector.}
  \item{delay_req}{Seconds between the single Blast requests.}
  \item{delay_rid}{Seconds between the single result requests.}
  \item{email}{User email, required information from NCBI (\code{String}).}
  \item{xmlFolder}{Path to the result folder.}
  \item{keepInMemory}{Logical, shall the results be kept in the memory.}
  \item{database}{The NCBI database to use.}
  \item{verbose}{Shall the program give extensive feedback.}
}

\details{
 This function sends fasta sequences to the NCBI blast service. The defaults for the delays are required by NCBI and must not be
 smaller than the default values. Also, NCBI asks the user to provide an email address.
 
 The input \code{seq} can be a vector of strings. In that case the sequences are one after another processed. The option \code{n_blast}
 sets then the upper threshold of how many blast requests are send to the NCBI Blast service at a time and kept running there parallel.
 It is here in the users obligation not to misuse the service with too many parallel requests. 

 The \code{xmlFolder} parameter specifies the folder to where the XML results will be stored.  In case the folder does not exist, R will create it.

 In case the option \code{keepInMemory} is set to TRUE the Blast results will be kept in memory, otherwise they will
 be just written to the HDD. Especially if many sequences are used to Blast it is recommended to drop the result from the memory.
}

\value{
  An xml file that contains the the NCBI result.
}

\author{
  Daniel Fischer
}

\examples{
\dontrun{
blastSeq("ACGTGCATCGACTAGCTACGACTACGACTATC")
}
}

\keyword{methods}
 
