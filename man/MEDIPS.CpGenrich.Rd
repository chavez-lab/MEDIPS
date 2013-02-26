\name{MEDIPS.CpGenrich}
\alias{MEDIPS.CpGenrich}
\title{
Calculates CpG enrichment of provided short reads compared to the reference genome.
}
\description{
As a quality check for the enrichment of CpG rich DNA fragments obtained by the immunoprecipitation step of a MeDIP experiment, this function provides the functionality to calculate CpG enrichment values. 
The main idea is to check, how strong the regions are enriched for CpGs compared to the reference genome.
For this, the function counts the number of Cs, the number of Gs, the number CpGs, and the total number of bases within the stated reference genome. 
Subsequently, the function calculates the relative frequency of CpGs and the observed/expected ratio of CpGs present in the reference genome.
Additionally, the function calculates the same for the DNA sequences underlying the given regions.
The final enrichment values result by dividing the relative frequency of CpGs (or the observed/expected value, respectively) of the regions by 
the relative frequency of CpGs (or the observed/expected value, respectively) of the reference genome.
}
\usage{
MEDIPS.CpGenrich(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=TRUE, chr.select=NULL, paired=F)
}
\arguments{
   \item{file}{
Path and file name of the input data
}
  \item{BSgenome}{
The reference genome name as defined by BSgenome 
}
  \item{extend}{
defines the number of bases by which the region will be extended before the genome vector is calculated.
Regions will be extended along the plus or the minus strand as defined by their provided strand information.
}
\item{shift}{
As an alternative to the extend parameter, the shift parameter can be specified. Here, the reads are not extended but shifted by the specified number of nucleotides with respect to the given strand infomation.
One of the two parameters extend or shift has to be 0. 
}
\item{uniq}{
MEDIPS will replace all reads which map to exactly the same start and end positions by only one representative, if uniq=TRUE
}
\item{chr.select}{
only data at the specified chromosomes will be processed.
}
\item{paired}{
option for paired end reads
}
}
\value{
\item{regions.CG}{the numbe of CpGs within the regions}
\item{regions.C}{the number of Cs within the regions}
\item{regions.G}{the number of Gs within the regions}
\item{regions.relH}{the relative frequency of CpGs within the regions}
\item{regions.GoGe}{the observed/expected ratio of CpGs within the regions}
\item{genome.CG}{the numbe of CpGs within the reference genome}
\item{genome.C}{the number of Cs within the reference genome}
\item{genome.G}{the number of Gs within the reference genome}
\item{genome.relH}{the relative frequency of CpGs within the reference genome}
\item{genome.GoGe}{the observed/expected ratio of CpGs within the reference genome}
\item{enrichment.score.relH}{regions.relH/genome.relH}
\item{enrichment.score.GoGe}{regions.GoGe/genome.GoGe}
}
\author{
Joern Dietrich and Matthias Lienhard
}
\examples{

library(MEDIPSData)
library("BSgenome.Hsapiens.UCSC.hg19")
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")

#er=MEDIPS.CpGenrich(file=bam.file.hESCs.Rep1.MeDIP, BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select="chr22", extend=0, shift=0, uniq=TRUE)
}
