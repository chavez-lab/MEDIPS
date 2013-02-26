\name{MEDIPS.selectROIs}
\alias{MEDIPS.selectROIs}
\title{
Selects row-wise subsets of a result table as returned by the MEDIPS.meth function.
}
\description{
MEDIPS provides the functionality to select subsets of the result matrix returned by the MEDIPS.meth function according to any given set of regions of interest (ROIs).
}
\usage{
MEDIPS.selectROIs(results=NULL, rois=NULL, columns=NULL, summarize=F)
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
It is possible to specify concrete column names or to provide several identifiers of interest like e.g. c("counts", "rpkm").
}
\item{summarize}{
As an alternative for extracting all windows included in the genomic ranges of the given ROIs,
it is possible to calculate mean values over the individual windows for each ROI- a behaviour that is controlled by the parameter summarize.
For each ROI and for each specified column, there will be one mean value returned, if summarize=TRUE.
}
}
\author{
Lukas Chavez
}
\examples{

library(MEDIPSData)
data(resultTable)

rois=data.frame(chr=c("chr22","chr22"), start=c(19136001, 19753401), stop=c(19136200, 19753500), ID=c("ID_1", "ID_2"))
columns=names(resultTable)[grep("counts|rpkm|logFC",names(resultTable))]
s = MEDIPS.selectROIs(results=resultTable, rois=rois, columns=columns, summarize=FALSE)


}


