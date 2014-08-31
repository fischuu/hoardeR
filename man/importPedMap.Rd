\name{importPedMap}
\alias{importPedMap}

\title{
  Import a ped/map File Pair
}

\description{
  This function imports a ped/map file pair.
}
\usage{
  importPedMap(ped, map=NULL, pedSep="\t", pedHeader=FALSE, genoSep=" ",
               mapSep="\t", mapHeader=FALSE, na.value="0")
}

\arguments{
  \item{ped}{File name of the ped file.}
  \item{map}{File name of the map file, optional, see details.}
  \item{pedSep}{Column separator in the ped file.}
  \item{pedHeader}{Logical, ped file contains header.}
  \item{genoSep}{Separator for Genotype, see details.}
  \item{mapSep}{Column separator in the map file.}
  \item{mapHeader}{Logical, map file contains header.}
  \item{na.value}{Character, encoding of missing values.}
  
}

\details{
  This function imports a ped/map file pair. For that it is sufficient to provide the file name of the ped file, if the
  map file has the same name, but just the .map ending (e.g. myFile.ped and myFile.map). Also, the file suffix .ped can
  be ommited. 
  
  The \code{genoSep} option provides the separator between the Alleles within one Genotype, e.g 'A A' (\code{genoSep=" "}) or
  'A/A' (\code{genoSep="/"}).
}
\value{
  A list of type \code{pedMap} with the three list items:
  
  \item{map}{Matrix with the Genotype Map information.}
  \item{fam}{Matrix with the family information.}
  \item{geno}{Matrix with the genotype information.}
}

\author{
  Daniel Fischer
}

 
