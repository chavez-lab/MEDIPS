\name{MEDIPS.coverageAnalysis}
\alias{MEDIPS.coverageAnalysis}
\title{
The function identifies the number of CpGs (or any other predefined sequence pattern) covered by the given short reads.
}
\description{
The main idea of the coverage analysis is to test the number of CpGs (or any other predefined sequence pattern) covered by the given short reads and to have a look at the depth of coverage.
Before you can start the coverage analysis, it is necessary that you have created a MEDIPS SET and excecuted the MEDIPS.getPositions function.
For the coverage analysis, the total set of available regions is divided into random subsets of equal size where the number of subsets is determined by the parameter no_iterations. 
The coverage analysis iteratively selects an increasing number of subsets and and tests how many CpGs are covered by the available regions.
Moreover, it is tested how many CpGs are covered at least 1x, 2x, 3x, 4x, 5x, and 10x. 
These levels of coverage depths can be adjusted by setting the coverages parameter (see below).
As the regions are typically of short length (e.g. 36bp), it is recommended to extend the region length by an extend value.
}
\usage{
MEDIPS.coverageAnalysis(data = NULL, coverages = c(1, 2, 3, 4, 5, 10), no_iterations = 10, no_random_iterations = 1, extend = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{coverages}{
default is c(1, 2, 3, 4, 5, 10). The coverages define the depth levels for testing how often a CpG was covered by the given regions.
Just specify any other vector of coverage depths you would like to test.
}
  \item{no_iterations}{
defines the number of subsets created from the full set of available regions (default=10). 
}
  \item{no_random_iterations}{
approaches that randomly select data entries may be processed several times in order to obtain more stable results. 
By specifying the no_random_iterations parameter (default=1) it is possible to run the coverage analysis several times. 
The final results returned to the coverage results object are the averaged results of each random iteration step.  
}
  \item{extend}{
extends the region lengths before the coverage analysis is performed.
}
}
\value{
\item{matrix}{Contains the number of covered CpGs in each iteration (rows) and for different levels of coverages (columns)}
\item{maxPos}{is the total number of sequence patterns (e.g. CpGs) within the refernce genome}
\item{pattern}{is the defined sequence pattern}
\item{coveredPos}{shows the number of covered sequence pattern (e.g. CpGs) using the total set of available regions for several depths of coverages (columns).
The last row shows the percentage of covered sequence pattern relative to the total number of available sequence patterns within the reference genome.}
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

cr.control

}

