\name{plotHit}
\alias{plotHit}

\title{
Visualization of a cross-species hit
}

\description{
For each cross-species hit the function plots the similarity within that area together with an optional annotation and coverage track.
}
\usage{
  plotHit(hits, flanking=1, window=NULL, annot=TRUE, coverage=FALSE,
          smoothPara=NULL, diagonal=0.25, verbose=TRUE, output=FALSE,
          hitSpecies=NULL, hitSpeciesVersion=NULL, origSpecies=NULL,
          origSpeciesVersion=NULL, fastaFolder=NULL, release=84,
          origAnnot=NULL, hitAnnot=NULL, nTick=5, which=NULL, figureFolder=NULL,
          figurePrefix=NULL, indexOffset=0, bamFolder=NULL, bamFiles=NULL,
          groupIndex=NULL, groupColor=NULL, countWindow=NULL)
}

\arguments{
  \item{hits}{The hit object to be plotted.}
  \item{flanking}{Allowed flanking site in Mb.}
  \item{window}{Moving window size of similarity measure.}
  \item{annot}{Logical, add annotation track}
  \item{coverage}{Logical, add coverage track} 
  \item{smoothPara}{Smoothing parameter for coverage}
  \item{diagonal}{Threshold for allowed diagonal similarity} 
  \item{verbose}{Logical, shall the function give status updates}
  \item{output}{Logical, shall numerical results be given}
  \item{hitSpecies}{Scientific identifier of the hit species.}
  \item{hitSpeciesVersion}{Version of the hit species assembly} 
  \item{origSpecies}{Scientific name of the original species}
  \item{origSpeciesVersion}{Version of the original species} 
  \item{fastaFolder}{Location of the fasta files}
  \item{release}{Ensembl release version}
  \item{origAnnot}{Annotation object of the original species}
  \item{hitAnnot}{Annotation object of the hit species} 
  \item{nTick}{Number of ticks on the annotation track}
  \item{which}{Which hits should be plotted}
  \item{figureFolder}{Folder where Figures should be stored}
  \item{figurePrefix}{Prefix of the figure filenames}
  \item{indexOffset}{Offset of the running index of the filenames}
  \item{bamFolder}{Folder with the bam-files}
  \item{bamFiles}{Filenames of the bam-files}
  \item{groupIndex}{Index of subgroups in the bamfiles}
  \item{groupColor}{Vector with colors, one for each subgroup}
  \item{countWindow}{Window size to count the reads from bam-files.}
}

\details{
 This function is the workhorse of hoardeR and visualizes the findings of the blast and intersection runs. It is really flexibel to handle the hits and  
 hence there are many different options. The required options are \code{hits}, \code{hitSpecies}, \code{origSpecies} and \code{fastaFolder}.
 
 The hit object is an object as provided by \code{intersectXMLAnnot} and contains all intersections of interest (=intersections that are in close
 proximity of a gene in the hit species). Naturally the hit and the original species have to be specified as well as the folder, where the required fasta
 files are stored, or to where they should  be downloaded. If the species are the default species from Ensembl (as can be seen in the data.frame
 \code{species}), the annotation and assembly will be automatically downloaded to the specified location on the harddrive. Changes from that
 version can be adjusted with the \code{release} and the \code{version} option, but the filenames have still to match the convention, as they
 are provided by Ensembl. For example the human assembly for chromosome 18 is labeled as
 
 \url{Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz}
 
 That means, fasta files need to be gzipped and they have to be available chromosome-wise. First there is the scientific species identifier, then 
 the assembly version, followed by a fixed string 'dna.chromosome'. With that schema any possible fasta file can be provided to \code{plotHit}. If in 
 addition to the similarity also a coverage track should be added, the option \code{coverage} has to be set to \code{TRUE}. The option 
 \code{smoothPara} sets then the level of smoothing of the coverage. By default no smoothing will be applied. 
 
 In case an annotation track is requested (\code{annot=TRUE}), the annotation objects need to be provided to the \code{origAnnot} and \code{hitAnnot} options.
 
 The option \code{diagonal} defines the minimum level of similarity so that a (diagonal) match will be plotted. The colors are then towards green for
 total similarity and towards red for total disagree, based on a nucleotide mismatch matrix.
 
 If the option \code{verbose=TRUE} is set, the function gives a verbose output while running. Further, if \code{output=TRUE} then, in addition to the
 figure also a data.frame with the numerical results is provided. 
 
 In case that \code{hits} contains more than one hit, the \code{plotHit} function plots for each hit a figure. In that case a folder should be
 provided to where the figures should be stored, this can be done with the \code{figureFolder} and \code{figurePrefix} options. In case only
 asserted hits of \code{hits} shall be plotted, they can be selected with the \code{which} option.
 
 The function can also plot a coverage track over the similarity. For that, the option \code{coverage=TRUE} has to be set and a folder that 
 contains the necessary bam-files has to be specified in \code{bamFolder}. By default all bam files in that folder are used, if only a subset
 is requested, the filenames can be specified in \code{bamFiles}. In case several bam-files are given, the average coverage at each loci is used.
 Further, if the data contains subgroups (e.g. case/control), the vector \code{groupIndex} gives the group labels. Naturally its length should be
 similar to \code{bamFiles} (or similar to the total amount of files in the bam-folder). In case that more than one group is plotted in the
 coverage track, their colors can be defined in \code{groupColor}. Of course, this vector has to be as long as the number of groups are defined. 
 The option \code{countWindow} controls the moving window length in which the number of counts is calculated. The default is the same length as the
 hit.
 
}
\value{
Optional, a table with intersection loci.
}

\author{
  Daniel Fischer
}


\examples{
\dontrun{
pigInter.flank <- list()
for(i in 1:nrow(pigHits)){
   pigInter.flank[[i]] <- intersectXMLAnnot(pigHits[i,], ssannot, flanking=100)
}
# Basic usage:
plotHit(hits=pigInter.flank,
        flanking=100,
        hitSpecies = "Sus scrofa",
        origSpecies = "Bos taurus",
        fastaFolder = "/home/user/fasta/",
        figureFolder = "/home/user/figures/") 

# Annotation tracks added:
plotHit(hits=pigInter.flank,
        flanking=100,
        hitSpecies = "Sus scrofa",
        origSpecies = "Bos taurus",
        fastaFolder = "/home/user/fasta/",
        figureFolder = "/home/user/figures/",
        origAnnot=btannot,
        hitAnnot=ssannot)
        
# Annotation and coverage added:
plotHit(hits=pigInter.flank,
        flanking=100,
        hitSpecies = "Sus scrofa",
        origSpecies = "Bos taurus",
        fastaFolder = "/home/daniel/fasta/",
        figureFolder = "/home/user/figures/",
        origAnnot=btannot,
        hitAnnot=ssannot
        coverage=TRUE,
        bamFolder = "/home/users/bams/") 
}
}

\keyword{methods}
 
