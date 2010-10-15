\name{MEDIPS.plotSaturation}
\alias{MEDIPS.plotSaturation}
\title{
Function plots the results of the MEDIPS.saturationAnalysis function.
}
\description{
The results of the saturation analysis will be visualized by the function.
}
\usage{
MEDIPS.plotSaturation(saturationObj = NULL)
}
\arguments{
  \item{saturationObj}{
The saturation results object returned by the MEDIPS.saturationAnalysis function
}
}
\value{
The coverage plot will be visualized.
}
\author{
Lukas Chavez
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)
sr.control = MEDIPS.saturationAnalysis(data = CONTROL.SET, bin_size = 50, extend = 400, no_iterations = 10, no_random_iterations = 1)

MEDIPS.plotSaturation(sr.control)

}
