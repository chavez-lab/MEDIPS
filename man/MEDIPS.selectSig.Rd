\name{MEDIPS.selectSig}
\alias{MEDIPS.selectSig}
\title{
Selects windows which show significant differential coverage between two MEDIPS SETs from a resultTable (as returned by the function MEDIPS.meth).
}
\description{
Based on the results table returned by the MEDIPS.meth function,
the function selects windows which show significant differential coverage between the two groups of MEDIPS SETs.
Selection of significant windows follows according to the specification of the available parameters.
}
\usage{
MEDIPS.selectSig(results=NULL, p.value=0.01, adj=T, ratio=NULL, bg.counts=NULL, CNV=F)
}
\arguments{
  \item{results}{
specifies the result table derived from the MEDIPS.meth function.
}
  \item{p.value}{
this is the p.value threshold as calculated either by the ttest or edgeR method
}
  \item{adj}{
this parameter specifies whether the p.value or the adjusted p.values is considered
}
\item{ratio}{
this parameter sets an additional thresold for the ratio where the ratio is either score.log2.ratio or edgeR.logFC depending on the previously selected method.
Please note, the specified value will be transformed into log2 internally.
}
  \item{bg.counts}{
as an additional filter parameter, it is possible to require a minimal number of reads per window in at least one of the MEDIPS SET groups.
For this, the mean of counts per group is considered. The parameter bg.counts can either be a concrete integer or an appropriate column name of the result table.
By specifying a column name, the 0.95 quantile of the according genome wide count distribution is determined and used as a minimal background threshold 
(please note, only count columns are reasonable).  
}
  \item{CNV}{
The information on CNVs present in the samples of interest can be used for correcting differential coverage observed in the corresponding IP data (e.g. MeDIP or ChIP data).
In case Input data has been provided for both conditions, MEDIPS is capable of calculating genome wide CNV ratios by employing the package DNAcopy.
In case the parameter CNV is set to TRUE, MEDIPS will subtract the CNV ratio from the IP ratio. 
Subsequently, only genomic windows having a CNV corrected IP ratio higher than the specified ratio threshold (specification of the ratio parameter is required in this case) 
will be considered as windows with sufficient differential IP coverage.
}
}
\author{
Lukas Chavez, Matthias Lienhard
}
\examples{

library(MEDIPSData)
data(resultTable)

sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)

}


