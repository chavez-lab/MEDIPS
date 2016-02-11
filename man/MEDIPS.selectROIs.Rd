\name{MEDIPS.selectROIs}
\alias{MEDIPS.selectROIs}
\title{
Selects row-wise subsets of a result table as returned by the MEDIPS.meth function.
}
\description{
MEDIPS provides the functionality to select subsets of the result matrix returned by the MEDIPS.meth function according to any given set of regions of interest (ROIs).
}
\usage{
MEDIPS.selectROIs(results=NULL, rois=NULL, columns=NULL, summarize=NULL)
}
\arguments{
  \item{results}{
a result table as returned by the function MEDIPS.meth
}
  \item{rois}{
A matrix containing genomic coordinates of any regions of interest.
}
  \item{columns}{
Only selected columns will be returned as determined by the columns parameter.
It is possible to specify one or more concrete column names, please see an example below.
}
\item{summarize}{
By setting summarize=NULL (default) all windows included in the genomic ranges of the given ROIs will be returned.
As an alternative, it is possible to calculate mean (summarize = "avg") or sum (summarize = "sum") values over the individual windows included in each ROI,
or to select only the most significant window within each given ROI (summarize="minP").
}
}
\author{
Lukas Chavez, Matthias Lienhard
}
\examples{

library(MEDIPSData)
data(resultTable)

rois=data.frame(chr=c("chr22","chr22"), start=c(19136001, 19753401), stop=c(19136200, 19753500), ID=c("ID_1", "ID_2"))
columns=names(resultTable)[grep("counts|rpkm|logFC",names(resultTable))]
s = MEDIPS.selectROIs(results=resultTable, rois=rois, columns=columns, summarize=NULL)


}


