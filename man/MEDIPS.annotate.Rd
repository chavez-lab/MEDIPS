\name{MEDIPS.annotate}
\alias{MEDIPS.annotate}
\title{
Funtion to annotate given genomic coordinates.
}
\description{
This function has been deprecated. 
Please see MEDIPS.getAnnotation and MEDIPS.setAnnotation instead.
}
\usage{
MEDIPS.annotate(region, anno)
}
\arguments{
  \item{region}{
a matrix that contains row-wise genomic regions, e.g. DMRs.
The columns are: chromosome, start,  stop.
}
  \item{anno}{
the annotation data object contains row-wise the genomic coordinates of annotations.
The columns are: chromosome, start, stop, ID
}
}
\value{
The annotation function returns a matrix where the rows contain the regions from the given frames object (here DMRs) and the columns are:
\item{chr}{the chromosome name of the DMR}
\item{start}{the start position of the DMR}
\item{stop}{the stop position of the DMR}
\item{annotation}{the name of the annotation}
}
\author{
Joern Dietrich, Matthias Lienhard
}

\examples{
print("The function has been deprecated. Please see MEDIPS.getAnnotation and MEDIPS.setAnnotation instead.")
}
