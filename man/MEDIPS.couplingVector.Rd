\name{MEDIPS.couplingVector}
\alias{MEDIPS.couplingVector}
\title{
Calculates the sequence pattern densities at genome wide windows.
}
\description{
The function calculates the local densities of a defined sequence pattern (e.g. CpGs)
and returns a COUPLING SET object which is necessary for normalizing MeDIP data.
}
\usage{
MEDIPS.couplingVector(pattern="CG", refObj=NULL)
}
\arguments{
  \item{pattern}{
defines the sequence pattern, e.g. CG for CpGs.
}
  \item{refObj}{
a MEDIPS Set or MEDIPS ROI Set that serves as reference for the genome and window parameters. 
}
}
\value{
A COUPLING SET object.
}
\author{
Lukas Chavez
}
\examples{

library("MEDIPSData")
library("BSgenome.Hsapiens.UCSC.hg19")

data(hESCs_MeDIP)
CS = MEDIPS.couplingVector(pattern="CG", refObj=hESCs_MeDIP)

}

