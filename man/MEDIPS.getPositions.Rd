\name{MEDIPS.getPositions}
\alias{MEDIPS.getPositions}
\title{
Identifies genomic sequence pattern positions within the reference genome. 
}
\description{
The function identifies the genomic positions of the stated sequence pattern (e.g. CpGs).
For sequence pattern that are reverse complementary, only the positions on the plus strand will be returned.
Otherwise, all genomic positions of the pattern on the plus and minus strand will be returned.
The reference genome is the genome (or only some chromosomes of a genome) that was specified by excecuting the MEDIPS.readAlignedSequences function. 
}
\usage{
MEDIPS.getPositions(data = NULL, pattern = NULL)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{pattern}{
defines the sequence pattern, e.g. CG for CpGs.
}
}
\value{
The slots of the stated MEDIPS SET object associated to the sequence pattern will be occupied afterwards.
These are informations about the pattern itself and their chromosome and genomic positions.
}
\author{
Joern Dietrich
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)

CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")

}

