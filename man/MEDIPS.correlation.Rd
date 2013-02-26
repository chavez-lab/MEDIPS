\name{MEDIPS.correlation}
\alias{MEDIPS.correlation}
\title{
Calculates pairwise Pearson correlations between provided MEDIPS SETs
}
\description{
The function calculates genome wide Pearson correlations between all pairs of provided MEDIPS SETs.
}
\usage{
MEDIPS.correlation(MSets=NULL, plot = T)
}
\arguments{
  \item{MSets}{
a concatenated set of MEDIPS SETs
}
  \item{plot}{
if specified, the correlation will be depicted as a scatter plot
}
}
\value{
a correlation matrix
}
\author{
Lukas Chavez
}
\examples{

library(MEDIPSData)
data(hESCs_MeDIP)
data(DE_MeDIP)

correlation = MEDIPS.correlation(MSets=c(hESCs_MeDIP[[1]], DE_MeDIP[[1]]), plot = FALSE)

}

