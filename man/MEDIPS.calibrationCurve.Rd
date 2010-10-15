\name{MEDIPS.calibrationCurve}
\alias{MEDIPS.calibrationCurve}
\title{
Function that calculates the calibration curve
}
\description{
Based on the calculated genome vector (MEDIPS.genomeVector) and on the coupling vector (MEDIPS.couplingVector) of a MEDIPS SET,
the function examins the dependency of local MeDIP-Seq signal intensities and local pattern (e.g. CpG) densities.
Calculation of the calibration curve is achieved by first dividing the total range of coupling factors into several levels.
Second, all genomic bins are partitioned into these levels by considering their associated coupling factors.   
Finally, for each level of coupling factors, MEDIPS calculates the mean signal and mean coupling factor of all genomic bins that fall into this level.
The calibration curve represents these averaged signals and coupling factors over the full range of coupling factors. 
It indicates the experiment specific dependency between signal intensity and CpG density.
Subsequently, the function performs a linear regression for small coupling factors of the calibration curve and records the slope and intercept of the resulting linear curve.
}
\usage{
MEDIPS.calibrationCurve(data = NULL)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
}
\value{
The slots of the stated MEDIPS SET object associated to the calibration curve will be occupied afterwards. 
These are the informations about the mean signals and mean coupling factors representing the calibration curve and the estimated normalization parameters intercept and slope.
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
}

