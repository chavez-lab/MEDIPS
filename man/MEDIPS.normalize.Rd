\name{MEDIPS.normalize}
\alias{MEDIPS.normalize}
\title{
Function that normalizes raw signals by local sequence pattern (e.g. CpG) densities.
}
\description{
The normalization function accesses the pre-calculated slope and intercept values derived from the MEDIPS.calibrationCurve function in order to weight the raw signals. 
The relative methlyation score (rms) for the genomic bins is then defined by rms = x\((y-intercept)/slope),
where x is the raw signal and y is the coupling factor of a genomic bin. 
Based on the total number of regions within the MEDIPS SET, the rms values will be transformed into a reads per million format and afterwards transformed into the log2 scale.
In order to make the rms values visualizable by common genome browsers, MEDIPS transforms its resulting data range into the consistent interval [0, 1000] before finally stored.
}
\usage{
MEDIPS.normalize(data = NULL)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
}
\value{
The slot of the stated MEDIPS SET object associated to the rms values will be occupied afterwards. 
}
\author{
Lukas Chavez
}

\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)
CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, bin_size = 50, extend = 400)
CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")
CONTROL.SET = MEDIPS.couplingVector(data = CONTROL.SET, fragmentLength = 700, func = "count")
CONTROL.SET = MEDIPS.calibrationCurve(data = CONTROL.SET)

CONTROL.SET = MEDIPS.normalize(data = CONTROL.SET)
}
