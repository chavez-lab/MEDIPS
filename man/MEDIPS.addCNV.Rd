\name{MEDIPS.addCNV}
\alias{MEDIPS.addCNV}
\title{
Function to run a copy number variation analysis.
}
\description{
Function calculates a CNV analysis based on two INPUT SETs by employing the DNAcopy package.
The results are attached to a provided result table.
}
\usage{
MEDIPS.addCNV(ISet1, ISet2, results, cnv.Frame=1000)
}
\arguments{
  \item{ISet1}{
First group of INPUT SETs
}
  \item{ISet2}{
Second group of INPUT SETs
}
\item{results}{
result table as returned by the MEDIPS.meth function
}
\item{cnv.Frame}{
window size used for calculating CNVs. Can be of different size than the result table.
}
}
\value{
The result table with an additional column containing DNAcopy's log-ratio.
}
\author{
Joern Dietrich
}

\examples{
library(MEDIPSData)
library("BSgenome.Hsapiens.UCSC.hg19")

bam.file.hESCs.Input = system.file("extdata", "hESCs.Input.chr22.bam", package="MEDIPSData")
bam.file.DE.Input = system.file("extdata", "DE.Input.chr22.bam", package="MEDIPSData")

hESCs.Input = MEDIPS.createSet(file=bam.file.hESCs.Input, BSgenome="BSgenome.Hsapiens.UCSC.hg19", extend=250, shift=0, uniq=TRUE, window_size=100, chr.select="chr22")
DE.Input = MEDIPS.createSet(file=bam.file.DE.Input, BSgenome="BSgenome.Hsapiens.UCSC.hg19", extend=250, shift=0, uniq=TRUE, window_size=100, chr.select="chr22")

data(resultTable)

resultTable = MEDIPS.addCNV(cnv.Frame=10000, ISet1=hESCs.Input, ISet2=DE.Input, results=resultTable)
}
