\name{MEDIPS.plotCoverage}
\alias{MEDIPS.plotCoverage}
\title{
Function plots the results of the MEDIPS.coverageAnalysis function.
}
\description{
The results of the coverage analysis will be visualized by the function.
}
\usage{
MEDIPS.plotCoverage(coverageObj = NULL)
}
\arguments{
  \item{coverageObj}{
The coverage results object returned by the MEDIPS.coverageAnalysis function
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
CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")
cr.control = MEDIPS.coverageAnalysis(data = CONTROL.SET, extend = 400, no_iterations = 10)

MEDIPS.plotCoverage(cr.control)

}
