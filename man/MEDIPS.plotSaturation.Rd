\name{MEDIPS.plotSaturation}
\alias{MEDIPS.plotSaturation}
\title{
Function plots the results of the MEDIPS.saturationAnalysis function.
}
\description{
The results of the saturation analysis will be visualized by the function.
}
\usage{
MEDIPS.plotSaturation(saturationObj = NULL, main="Saturation analysis")
}
\arguments{
  \item{saturationObj}{
The saturation results object returned by the MEDIPS.saturationAnalysis function
}
\item{main}{
The title of the coverage plot.
}
}
\value{
The coverage plot will be visualized.
}
\author{
Lukas Chavez
}
\examples{

library(MEDIPSData)
library(BSgenome.Hsapiens.UCSC.hg19)
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")

sr=MEDIPS.saturation(file=bam.file.hESCs.Rep1.MeDIP, BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=TRUE, extend=250, shift=0, window_size=100, chr.select="chr22", nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
MEDIPS.plotSaturation(saturationObj = sr, main="Saturation analysis")

}
