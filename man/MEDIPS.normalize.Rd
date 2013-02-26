\name{MEDIPS.normalize}
\alias{MEDIPS.normalize}
\title{
Function that normalizes raw signals by local sequence pattern (e.g. CpG) densities.
}
\description{
This function has been deprecated. 
Please see MEDIPS.rms instead.

The normalization function accesses the pre-calculated slope and intercept values derived from the MEDIPS.calibrationCurve function in order to weight the raw signals. 
The relative methlyation score (rms) for the genomic bins is then defined by rms = x\((y-intercept)/slope),
where x is the raw signal and y is the coupling factor of a genomic bin. 
Based on the total number of regions within the MEDIPS SET, the rms values will be transformed into a reads per million format.
}
\usage{
MEDIPS.normalize(data=NULL)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
}
\value{
RMS signals.
}
\author{
Lukas Chavez
}

\examples{

print("This function has been deprecated. Please see MEDIPS.rms instead.")

}
