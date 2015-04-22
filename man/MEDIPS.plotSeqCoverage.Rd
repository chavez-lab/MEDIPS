\name{MEDIPS.plotSeqCoverage}
\alias{MEDIPS.plotSeqCoverage}
\title{
Function plots the results of the MEDIPS.seqCoverage function.
}
\description{
The results of the sequence pattern coverage analysis will be visualized in two possible ways.
}
\usage{
MEDIPS.plotSeqCoverage(seqCoverageObj=NULL, main=NULL, type="pie", cov.level = c(0,1,2,3,4,5), t="Inf")
}
\arguments{
  \item{seqCoverageObj}{
The coverage results object returned by the MEDIPS.seqCoverage function.
}
\item{main}{
The title of the coverage plot.
}
\item{type}{
there are two types of visualization. The pie chart (default) illustrates the fraction of CpGs covered by the given reads at different coverage level (see also the parameter cov.level).
As an alternative, a histogram over all coverage level can be ploted ("hist").
}
\item{cov.level}{
The pie chart illustrates the fraction of CpGs covered by the given reads according to their coverage level.
The visualized coverage levels can be adjusted by the cov.level parameter.
}
\item{t}{
specifies the maximal coverage depth to be plotted, if type="hist"
}
}
\value{
The sequence pattern coverage plot will be visualized.
}
\author{
Lukas Chavez
}
\examples{

library(MEDIPSData)
library(BSgenome.Hsapiens.UCSC.hg19)
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")

cr=MEDIPS.seqCoverage(file=bam.file.hESCs.Rep1.MeDIP, pattern="CG", BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select="chr22", extend=250, shift=0, uniq=1e-3, paired=FALSE)

MEDIPS.plotSeqCoverage(seqCoverageObj=cr, main="Sequence pattern coverage", type="pie", cov.level = c(0,1,2,3,4,5))


}

