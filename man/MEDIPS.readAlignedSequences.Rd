\name{MEDIPS.readAlignedSequences}
\alias{MEDIPS.readAlignedSequences}
\title{
Creates a MEDIPS SET by reading a suitable input file
}
\description{
Reads the input file and creates a MEDIPS SET.
After reading the input file, the MEDIPS SET contains the information about the input regions, like the input file name, 
the dependent organism, the chromosomes included in the input file, the length of the included chromosomes (automatically loaded), the number of regions, 
and the start, stop and strand informations of the regions. All further slots, for example for the weighting parameters and normalized data 
are still empty and will be filled during the workflow.
}
\usage{
MEDIPS.readAlignedSequences(file = NULL, BSgenome = NULL, numrows = -1)
}
\arguments{
  \item{file}{
Path and file name of the input data
}
  \item{BSgenome}{
The reference genome name as defined by BSgenome 
}
  \item{numrows}{
The number of short reads (number of rows) within the input file
}
}
\value{
An object of class MEDIPSset is returned where the region dependent informations are stored in the according slots.
These are informations about the input file, the reference genome, the total number of provided regions, the chromosomes which are covered by the regions, the total chromosome lengths, and the start and stop positions and strand informations of the regions.
}
\author{
Lukas Chavez
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")

CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)

}
