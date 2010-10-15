\name{MEDIPS.annotate}
\alias{MEDIPS.annotate}
\title{
Funtion to annotate given genomic coordinates.
}
\description{
The function annotates any matrix containing genomic coordinates (region) by a given annotation file (anno) containing genomic regions of interest.
During a typical MEDIPS workflow, this is of interest for anotating identified differentially methylated regions (DMRs)
derived after having xcecuted the MEDIPS.selectSignificants() or MEDIPS.mergeFrames() function. 
For annotating DMRs, you have to provide an annotation file that contains pre-defined ROIs.
For each provided region, the function returns all annotations from the provided annotation file.
In case there are several overlapping annotations, the region is returned several times in sparated rows, each row associated to one annotation.
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
Joern Dietrich
}

\examples{

region = list(chr="chr22", start=25170186, stop=25170687)
anno = system.file("extdata", "hg19.chr22.txt", package="MEDIPS")

annotated = MEDIPS.annotate(region=region, anno=anno)

annotated

}
