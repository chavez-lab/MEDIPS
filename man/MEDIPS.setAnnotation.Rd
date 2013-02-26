\name{MEDIPS.setAnnotation}
\alias{MEDIPS.setAnnotation}
\title{
Funtion to annotate a matrix of genomic coordinates (i.e. a result table) by a given annotation object.
}
\description{
The function appends any annotation IDs included in the given annotation object to the given regions object.
An annotation object can be retrived by the MEDIPS.getAnnotation function and the regions object is typically a (filtered) result table as returned by the MEDIPS.meth function.
An annotation ID is appended to a genomic region if their genomic coordinates overlap by at least one base.
There will be as many columns added to the regions object as overlapping annotations exist in the annotation object.
}
\usage{
MEDIPS.setAnnotation(regions, annotation, cnv=F)
}
\arguments{
  \item{regions}{
a matrix that contains row-wise genomic regions, e.g. as a result of the MEDIPS.meth function.
}
  \item{annotation}{
the annotation data object contains the genomic coordinates of annotations.
An annotation object can be e.g. retrived by the MEDIPS.getAnnotation function.
}
  \item{cnv}{
the MEDIPS.setAnnotation function is also internally used by the MEDIPS.addCNV function which automatically sets this parameter to TRUE. 
Otherwise cnv should be set to FALSE.
}
}
\value{
The provided result object with added columns containing overlapping annotations.
}
\author{
Joern Dietrich, Matthias Lienhard
}

\examples{
library(MEDIPSData)
data(resultTable)

sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)
sig = MEDIPS.mergeFrames(frames=sig, distance=1)
ens_gene = MEDIPS.getAnnotation( annotation="GENE",chr="chr22")
sig = MEDIPS.setAnnotation(regions=sig, annotation=ens_gene)

}
