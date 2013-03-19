\name{MEDIPS.plotCalibrationPlot}
\alias{MEDIPS.plotCalibrationPlot}
\title{
Creates the calibration plot.
}
\description{
Visualizes the dependency between MeDIP-seq signals and CpG densities together with the calcibration curve and the results of the linear modelling.
}
\usage{
MEDIPS.plotCalibrationPlot(MSet=NULL, ISet=NULL, CSet=NULL, plot_chr="chr1", rpkm=T, main="Calibration plot", xrange=T)
}
\arguments{
  \item{MSet}{
a MEDIPS SET object
}
\item{ISet}{
an INPUT SET (i.e. a MEDIPS SET created from Input sequening data)
}
\item{CSet}{
an according COUPLING SET object
}
\item{plot_chr}{
default="chr1". It is recommended to call a graphics device (e.g. png("calibrationPlot.png")) before calling the plot command, because R might not be able to plot the full amount of data in reasonable time.
}
\item{rpkm}{
can be either TRUE or FALSE. If set to TRUE, the methylation signals will be transformed into rpkm before plotted.
The coupling values remain untouched. It is necessary to set rpkm=T, if both, a MSet and an ISet are given and plotted at the same time.
}
\item{main}{
The title of the calibration plot.
}
\item{xrange}{
The signal range of the calibration curve typically falls into a low signal range. 
By setting the xrange parameter to TRUE, the calibration plot will visualize the low signal range only. 
}

}
\value{
The calibration plot will be visualized.
}
\author{
Lukas Chavez, Matthias Lienhard
}
\examples{
library(MEDIPSData)
data(hESCs_MeDIP)
data(CS)

MEDIPS.plotCalibrationPlot(CSet=CS, main="Calibration Plot", MSet=hESCs_MeDIP[[1]], plot_chr="chr22",  rpkm=TRUE,  xrange=TRUE)

}
