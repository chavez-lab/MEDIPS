\name{MEDIPS.mergeFrames}
\alias{MEDIPS.mergeFrames}
\title{
Merges genomic coordinates of neighboring windows into one supersized window
}
\description{
After having filtered the result table returned by the MEDIPS.meth function using the MEDIPS.selectSig function, there might be neighboring significant frames.
For these cases it is worthwhile to merge neighbouring regions into one supersized frame.
}
\usage{
MEDIPS.mergeFrames(frames=NULL, distance=1)
}
\arguments{
  \item{frames}{
  is a filtered result table received by the MEDIPS.selectSig function.
}
\item{distance}{
  allows an according number of bases as a gap between neighboring significant windows to be merged.
 The default value is 1 in order to merge adjacent windows.
}
}
\value{
The remaining distinct frames are represented only by their genomic coordinates within the returned results table
\item{chromosome}{the chromosome of the merged frame}
\item{start}{the start position of the merged frame}
\item{stop}{the stop position of the merged frame}
\item{ID}{a numbered ID of the merged frame}
The result table does not contain any merged significant values.
}
\author{
Lukas Chavez
}
\examples{

regions=as.data.frame(list(chr=c("chr22", "chr22"), start=c(1001, 1501), stop=c(1500,1750)))
regions.merged=MEDIPS.mergeFrames(regions)

regions.merged

}
