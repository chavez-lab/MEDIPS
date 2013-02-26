\name{MEDIPS.coverageAnalysis}
\alias{MEDIPS.coverageAnalysis}
\title{
The function identifies the number of CpGs (or any other predefined sequence pattern) covered by the given short reads.
}
\description{
This function has been deprecated. 
Please see MEDIPS.seqCoverage instead.
}
\usage{
MEDIPS.coverageAnalysis(data=NULL, coverages=c(1,2,3,4,5,10), no_iterations=10, no_random_iterations=1, extend=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
MEDIPS SET
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
print("The function has been deprecated. Please see MEDIPS.seqCoverage.")

}

