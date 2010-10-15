\name{MEDIPS.plotCalibrationPlot}
\alias{MEDIPS.plotCalibrationPlot}
\title{
Plots the results of the MEDIPS.calibrationCurve function.
}
\description{
Visualizes the dependency of raw MeDIP-Seq signals and CpG densities together with the results of the calcibration curve calculation.
}
\usage{
MEDIPS.plotCalibrationPlot(data=NULL, xrange=NULL, linearFit=FALSE, plot_chr="all", rpm=F)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{xrange}{
The mean signal range of the calibration curve typically falls into a low signal range. 
By setting the xrange parameter to e.g. 50, the calibration plot will only plot genomic bins associated with signals <=50. 
Therefore, the effect of an increased CpG density to an increased signal can be better visualized, especially if the data contains genomic bins with high signals.
}
  \item{rpm}{
can be either TRUE or FALSE. If set to TRUE, the signals will be transformed into reads per million (rpm) before plotted.
The coupling values remain untouched.
}
  \item{linearFit}{
When the parameter linearFit is set to TRUE, the plot contains the calculated linear curve that represents the dependency between signals and CpG densities.
}
  \item{plot_chr}{
default="all". Please don't forget to call a e.g. png("file.png") function before calling the plot command using "all" because R might not be able to plot the full amount of data in reasonable time.
Alternatively, you can specify a selected chromosome (e.g. chr1). 
Here, the plot_chr parameter only affects the plot and does not affect the MEDIPS SET.
}
}
\value{
The calibration plot will be visualized.
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

MEDIPS.plotCalibrationPlot(data = CONTROL.SET, linearFit = TRUE, plot_chr = "chr22")
}

