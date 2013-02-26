\name{MEDIPS.seqCoverage}
\alias{MEDIPS.seqCoverage}
\title{
The function identifies the number of CpGs (or any other predefined sequence pattern) covered by the given short reads.
}
\description{
The main idea of the sequence pattern coverage analysis is to test the number of CpGs (or any other predefined sequence pattern) covered by the given short reads and to test the depth of coverage.
}
\usage{
MEDIPS.seqCoverage(file = NULL, BSgenome = NULL, pattern = "CG", extend = 0, shift = 0, uniq = TRUE, chr.select = NULL, paired = F)
}
\arguments{
 \item{file}{
Path and file name of the input data
}
 \item{BSgenome}{
The reference genome name as defined by BSgenome 
}
 \item{pattern}{
defines the sequence pattern, e.g. CG for CpGs.
} 
  \item{extend}{
defines the number of bases by which the region will be extended before the genome vector is calculated.
Regions will be extended along the plus or the minus strand as defined by their provided strand information.
Please note, the extend and shift parameter are mutual exclusive.
}
  \item{shift}{
defines the number of bases by which the region will be shifted before the genome vector is calculated.
Regions will be shifted along the plus or the minus strand as defined by their provided strand information.
Please note, the extend and shift parameter are mutual exclusive.
}
  \item{uniq}{
if uniq=TRUE (default), all reads with exactly the same start and end positions will be replaced by one representative before the saturation analysis is performed.
}
  \item{chr.select}{
specify a subset of chromosomes for which the saturation analysis is performed.
}
\item{paired}{
option for paired end reads
}
}
\author{
Lukas Chavez
}
\examples{

library(MEDIPSData)
library(BSgenome.Hsapiens.UCSC.hg19)
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")

cr = MEDIPS.seqCoverage(file=bam.file.hESCs.Rep1.MeDIP, BSgenome="BSgenome.Hsapiens.UCSC.hg19", pattern="CG", extend=250, shift=0, uniq=TRUE, chr.select="chr22")

}

