\name{MEDIPS.mergeFrames}
\alias{MEDIPS.mergeFrames}
\title{
Merges genomic coordinates of overlapping frames into one supersized frame
}
\description{
In case, the MEDIPS.diffMethyl function was excecuted by setting the value of the step parameter < the value of frame_size the parameter,
one may end up with overlapping significant frames.
For these cases it is worthwhile to merge overlapping regions into one supersized frame.
}
\usage{
MEDIPS.mergeFrames(frames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frames}{
  is a matrix received by the MEDIPS.selectSignificants() function.
}
}
\value{
The remaining distinct frames are represented only by their genomic coordinates within the returned results table
\item{chromosome}{the chromosome of the merged frame}
\item{start}{the start position of the merged frame}
\item{stop}{the stop position of the merged frame}
The results table does not contain any merged rpm, rms, variance, p.value, etc. values.
}
\author{
Lukas Chavez
}
\examples{

regions=as.data.frame(list(chr=c("chr22", "chr22"), start=c(1000, 1250), stop=c(1500,1750)))
regions.merged=MEDIPS.mergeFrames(regions)

regions.merged

}
