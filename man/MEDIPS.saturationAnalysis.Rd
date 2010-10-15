\name{MEDIPS.saturationAnalysis}
\alias{MEDIPS.saturationAnalysis}
\title{
Function calculates the saturation/reproducibility of the provided MeDIP-Seq data.
}
\description{
The saturation analysis addresses the question, whether the number of input regions is sufficient to generate a saturated and reproducible methylation profile of the reference genome.
The main idea is that an insufficent number of short reads will not result in a saturated methylation profile. 
Only if there is a sufficient number of short reads, the resulting genome wide methylation profile will be reproducible by another independent set of a similar number of short reads.
}
\usage{
MEDIPS.saturationAnalysis(data = NULL, no_iterations = 10, no_random_iterations = 1, empty_bins = TRUE, rank = FALSE, extend = 400, bin_size = NULL)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{no_iterations}{
defines the number of subsets created from the full sets of available regions (default=10)
}
  \item{no_random_iterations}{
approaches that randomly select data entries may be processed several times in order to obtain more stable results. 
By specifying the no_random_iterations parameter (default=1) it is possible to run the saturation analysis several times. 
The final results returned to the saturation results object are the averaged results of each random iteration step. 
}
  \item{empty_bins}{
can be either TRUE or FALSE (default TRUE). This parameter effects the way of calculating correlations between the resulting genome vectors. 
A genome vector consists of concatenated vectors for each included chromosome. The size of the vectors is defined by the bin_size parameter.
If there occur genomic bins which contain no overlapping regions, neither from the subsets of A nor from the subsets of B, 
these bins will be neglected when the paramter is set to FALSE.
}
  \item{rank}{
can be either TRUE or FALSE (default FALSE). This parameter also effects the way of calculating correlations between the resulting genome vectors.
If rank is set to TRUE, the correlation will be calculated for the ranks of the bins instead of considering the counts.
Setting this parameter to TRUE is a more robust approach that reduces the effect of possible occuring outliers (these are bins with a very high number of overlapping regions) to the correlation.
}
  \item{extend}{
defines the number of bases by which the region will be extended before the genome vector is calculated.
Regions will be extended along the plus or the minus strand as defined by their provided strand information.
}
  \item{bin_size}{
defines the size of genome wide bins and therefore, the size of the genome vector. 
Read coverages will be calculated for bins separated by bin_size base pairs.
}
}
\value{
\item{distinctSets}{Contains the results of each iteration step (row-wise) of the saturation analysis.
The first column is the number of considered regions in each set, the second column is the resulting pearson correlation coefficient when comparing the two independent genome vectors.}
\item{estimation}{Contains the results of each iteration step (row-wise) of the estimated saturation analysis.
The first column is the number of considered regions in each set, the second column is the resulting pearson correlation coefficient when comparing the two independent genome vectors.}
\item{distinctSets}{the total number of available regions}
\item{maxEstCor}{contains the best pearson correlation (second column) obtained by considering the artifically doubled set of reads (first column)}
\item{distinctSets}{contains the best pearson correlation (second column) obtained by considering the total set of reads (first column)}
}
\author{
Lukas Chavez
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)

sr.control = MEDIPS.saturationAnalysis(data = CONTROL.SET, bin_size = 50, extend = 400, no_iterations = 10, no_random_iterations = 1)

sr.control

}
