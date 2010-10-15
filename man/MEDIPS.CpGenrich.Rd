\name{MEDIPS.CpGenrich}
\alias{MEDIPS.CpGenrich}
\title{
Calculates the enrichment of provided CpG rich regions compared to the reference genome.
}
\description{
As a quality check for the enrichment of CpG rich DNA fragments obtained by the immunoprecipitation step of a MeDIP experiment, this function provides the functionality to calculate CpG enrichment values. 
The main idea is to check, how strong the regions are enriched for CpGs compared to the reference genome.
For this, the function counts the number of Cs, the number of Gs, the number CpGs, and the total number of bases within the reference genome of the stated MEDIPS SET. 
Subsequently, the function calculates the relative frequency of CpGs and the observed/expected ratio of CpGs present in the reference genome.
Additionally, the function calculates the same for the DNA sequences underlying the given regions.
The final enrichment values result by dividing the relative frequency of CpGs (or the observed/expected value, respectively) of the regions by 
the relative frequency of CpGs (or the observed/expected value, respectively) of the reference genome.
}
\usage{
MEDIPS.CpGenrich(data = data, extend = NULL)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{extend}{
defines the number of bases by which the region will be extended before the genome vector is calculated.
Regions will be extended along the plus or the minus strand as defined by their provided strand information.
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
Joern Dietrich
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET=MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)

enrich=MEDIPS.CpGenrich(data=CONTROL.SET)

enrich

}
