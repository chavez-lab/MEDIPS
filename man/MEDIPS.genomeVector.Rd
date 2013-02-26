\name{MEDIPS.genomeVector}
\alias{MEDIPS.genomeVector}
\title{
Calculates the genome wide short read coverage on a user specified resolution 
}
\description{
This function has been deprecated. 
Please see MEDIPS.createSet instead.

Based on the regions included within a previously created MEDIPS SET (see MEDIPS.readAlignedSequiences),
the function calculates the genome wide coverage on a user specified resolution.
Each chromosome inside the MEDIPS SET will be divided into bins of size bin_size and the short read coverage will be calculated on this resolution.
The bin representation of the genome is the 'genome vector'.  
}
\usage{
MEDIPS.genomeVector(data = NULL, extend = 400, bin_size = 50)
}
\arguments{
  \item{data}{
MEDIPS SET
}
  \item{extend}{
Regions will be extended w.r.t. the extend parameter along the plus or the minus strand (as defined by their provided strand information).
After extending the regions, their length will be 'extend' (i.e. the extend parameter is NOT added to the given read lengths but all regions will be of size 'extend' afterwards).  
}
  \item{bin_size}{
defines the size of genome wide bins and therefore, the size of the genome vector. 
Read coverages will be calculated for bins separated by bin_size base pairs.
}
}
\value{
The slots of the stated MEDIPS SET object associated to the genome vector will be occupied afterwards. 
These are the informations about the bin_size, the extend value, the chromosome and position of the bins, and the number regions within the MEDIPS SET that overlap with the genomic bins.
}
\author{
Lukas Chavez
}
\examples{

print("This function has been deprecated. Please see MEDIPS.createSet instead.")

}

