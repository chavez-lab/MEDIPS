\name{MEDIPS.couplingVector}
\alias{MEDIPS.couplingVector}
\title{
Calculates the sequence pattern densities for genome wide bins.
}
\description{
Based on the coordinates of the bins of the genome vector included in the stated MEDIPS SET object,
the function calculates the local density of a pre-defined sequence pattern (e.g. CpGs).
}
\usage{
MEDIPS.couplingVector(data = NULL, distFile = "empty", fragmentLength = 700, func = "count")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{distFile}{
only accessed, if the parameter func=custom.
By setting the parameter func to custom, it is required to specify any custom distance weights file using the parameter distFile.
}
  \item{fragmentLength}{
Only sequence pattern (e.G. CpGs) within the range of (bin_position-fragmentLength), bin_position+fragmentLength] will contribute to the final local coupling factor.
The optimized value for the fragmentLength will reflect the estimated size of your sonicated DNA fragments.
}
  \item{func}{
There are several possible weghting function.
MEDIPS supports setting the weighting function parameter func to
count: simply count the number of CpGs within the predifined maximal distance to the current bin;
linear: the weights for CpGs decreases in a linear way and end at 0 at the predifined maximal distance to the current bin;
exp: the weights for CpGs decreases in an exponential way and end at 0 at the predifined maximal distance to the current bin;
log: the weights for CpGs decreases in a logarithmic way and end at 0 at the predifined maximal distance to the current bin;
custom: by setting the parameter to custom, it is required to specify a custom distance weights file using the parameter distFile.
You can create any of such a distance file by your own and specify it here. 
Here, the fragmentLength parameter will be neglected and the maximal distance within your provided distance file will be the limit.
}
}
\value{
The slots of the stated MEDIPS SET object associated to the coupling vector will be occupied afterwards. 
These are the informations about the selected distance function, possibly about the provided distance weights file, the fragment length and the calculated coupling factors for the genomic bins.
}
\author{
Lukas Chavez and Joern Dietrich
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)
CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, bin_size = 50, extend = 400)
CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")

CONTROL.SET = MEDIPS.couplingVector(data = CONTROL.SET, fragmentLength = 700, func = "count")

}

