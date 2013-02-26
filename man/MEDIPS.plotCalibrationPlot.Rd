\name{MEDIPS.plotCalibrationPlot}
\alias{MEDIPS.plotCalibrationPlot}
\title{
Creates the calibration plot.
}
\description{
Visualizes the dependency between MeDIP-seq signals and CpG densities together with the calcibration curve and the results of the linear modelling.
}
\usage{
MEDIPS.plotCalibrationPlot(MSet=NULL, CSet=NULL, plot_chr="chr1", rpkm=F, main="Calibration plot", xrange=T, input=FALSE)
}
\arguments{
  \item{MSet}{
a MEDIPS SET object
}
\item{CSet}{
an according COUPLING SET object
}
\item{plot_chr}{
default="chr1". It is recommended to call a graphics device (e.g. png("calibrationPlot.png")) before calling the plot command, because R might not be able to plot the full amount of data in reasonable time.
}
\item{rpkm}{
can be either TRUE or FALSE. If set to TRUE, the methylation signals will be transformed into rpkm before plotted.
The coupling values remain untouched.
}
\item{main}{
The title of the calibration plot.
}
\item{xrange}{
The signal range of the calibration curve typically falls into a low signal range. 
By setting the xrange parameter to TRUE, the calibration plot will visualize the low signal range only. 
}
\item{input}{
Set input to TRUE when you are processing Input data (without MeDIP enrichment). 
The function will not calculate the linear model which represents the dependency between MeDIP signals and CpG densities.
}

}
\value{
The calibration plot will be visualized.
}
\author{
Lukas Chavez
}
\examples{
library(MEDIPSData)
data(hESCs_MeDIP)
data(CS)

MEDIPS.plotCalibrationPlot(CSet=CS, input=FALSE, main="Calibration Plot", MSet=hESCs_MeDIP[[1]], plot_chr="chr22",  rpkm=TRUE,  xrange=TRUE)

}

