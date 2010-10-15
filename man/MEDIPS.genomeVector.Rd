\name{MEDIPS.genomeVector}
\alias{MEDIPS.genomeVector}
\title{
Calculates the genome wide short read coverage on a user specified resolution 
}
\description{
Based on the regions included within a previously created MEDIPS SET (see MEDIPS.readAlignedSequiences),
the function calculates the genome wide coverage on a user specified resolution.
Each chromosome inside the MEDIPS SET will be divided into bins of size bin_size and the short read coverage will be calculated on this resolution.
The bin representation of the genome is the 'genome vector'.  
}
\usage{
MEDIPS.genomeVector(data = NULL, extend = 400, bin_size = 50)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{extend}{
defines the number of bases by which the region will be extended before the genome vector is calculated.
Regions will be extended along the plus or the minus strand as defined by their provided strand information. 
}
  \item{bin_size}{
defines the size of genome wide bins and therefore, the size of the genome vector. 
Read coverages and coupling factors will be calculated for bins separated by bin_size base pairs.
}
}
\value{
The slots of the stated MEDIPS SET object associated to the genome vector will be occupied afterwards. 
These are the informations about the bin_size, the extend value, the chromosome and position of the bins, and the number regions within the MEDIPS SET that overlap with the genomic bin.
}
\author{
Lukas Chavez
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)

CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, bin_size = 50, extend = 400)

}

