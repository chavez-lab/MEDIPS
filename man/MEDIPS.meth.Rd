\name{MEDIPS.meth}
\alias{MEDIPS.meth}
\alias{MEDIPS.diffMeth}

\title{
Funtion summarizes coverage profiles for given MEDIPS SETs and allows to calculate differental coverage and copy number vartiation, if applicable. 
}
\description{
The function summarizes coverage profiles (counts, rpkm) for given MEDIPS SETs at the slots MSet1, MSet2, ISet1, and ISet1.
In case the parameter MeDIP is set to TRUE and a COUPLING SET was provided at the slot CS, 
the function will calculate CpG density normalized methylation profiles (relative methylation score, rms) for the MEDIPS SETs at the slots MSet1 and MSet2.
In case two groups of MEDIPS SETs have been provided at MSet1 and MSet2, the function will calculate differential coverage.
In case two groups of MEDIPS SETs have been provided at ISet1 and ISet2 and the parameter CNV was set to TRUE, 
the function will calculate copy number variation.
The function allows for processing a variable number of provided MEDIPS SETs and therefore, the returned matrix is of variable length.
}
\usage{
MEDIPS.meth(MSet1 = NULL, MSet2 = NULL, CSet = NULL, ISet1 = NULL, ISet2 = NULL, chr = NULL, p.adj="bonferroni", diff.method="edgeR", CNV=FALSE, MeDIP=FALSE, type="rpkm",minRowSum=10, quantile=FALSE)
}
\arguments{
  \item{MSet1}{
has to be one or a concatenated list of MEDIPS SET objects (the control replicates)
}
  \item{MSet2}{
has to be one or a concatenated list of MEDIPS SET objects (the treatment data) or empty
}
  \item{CSet}{
has to be a COUPLING SET object (must fit the given MEDIPS SET objects with respect to reference genome and represented chromosomes)
}
  \item{ISet1}{
has to be one or a concatenated list of Input derived MEDIPS SET objects (general Input data or Inputs from the control replicates) or empty
}
   \item{ISet2}{
has to be one or a concatenated list of Input derived MEDIPS SET objects (Inputs from the treatment replicates) or empty
}
   \item{chr}{
specify one or several chromosomes (e.g. c("chr1", "chr2")), if only a subset of available chromosomes have to be processed. 
}
  \item{p.adj}{
in order to correct p.values derived from the differential coverage analysis for multiple testing, MEDIPS uses Rs' p.adjust function.
Therefore, the following methods are available: holm, hochberg, hommel, bonferroni (default) , BH, BY, fdr, none.
}
  \item{diff.method}{
method for calculating differential coverage. Available methods: ttest (default) and edgeR.
}
  \item{CNV}{
In case there are INPUT SETs provided at both Input slots (i.e. ISet1 and ISet2), 
copy number variation will be tested by applying the package DNAcopy to the window-wise log-ratios calculated based on the the means per group. 
By setting CNV=F this function will be disabled (default: CNV=TRUE). 
Please note, there is the function MEDIPS.addCNV which allows to run the CNV analysis on two groups of INPUT SETs using another (typically increased) window size.
}
  \item{MeDIP}{
This parameter determines, whether for the MEDIPS SETs given at the slots MSet1 and MSet2, CpG density dependent normalization values (rms and prob) will be calculated (default: MeDIP=TRUE).  
}
  \item{type}{
In case diff.method has been set to ttest, this parameter specifies, if differential coverage is calculated based on the rpkm (default) or rms values.
This parameter is ignored in case the edgeR method is selected as the underlying model requires count data.
}
  \item{minRowSum}{
threshold for the sum of counts in a window for the staistical test (default=10). 
}
   \item{quantile}{
This parameter enables quantile normalization on the count matrix prior to calculating differential coverage (default=FALSE).
}
}

\value{
\item{Chr}{the chromosome of the ROI}
\item{Start}{the start position of the ROI}
\item{Stop}{the stop position of the ROI}
\item{CF}{the number of CpGs in the window}
\item{*counts}{a variable number of columns (according to the number of provided MEDIPS SETs) containing for each set the number of (extended/shifted) reads that overlap with the window.}
\item{*rpkm}{a variable number of columns (according to the number of provided MEDIPS SETs) containing for each set the rpkm value of the window.}
\item{*rms}{optional (if MeDIP=TRUE): a variable number of columns (according to the number of provided MEDIPS SETs) containing for each set the rms value of the window.}
\item{*counts}{optional (if INPUT SETs given): a variable number of columns (according to the number of provided INPUT SETs) containing for each set the number of (extended/shifted) reads that overlap with the window.}
\item{*rpkm}{optional (if INPUT SETs given): a variable number of columns (according to the number of provided INPUT SETs) containing for each set the rpkm value of the window.}
\item{MSets1.counts.mean}{optional (if more than one MEDIPS SET given): the mean count over all MEDIPS SETs at MSet1.}
\item{MSets1.rpkm.mean}{optional (if more than one MEDIPS SET given): the mean rpkm value over all MEDIPS SETs at MSet1.}
\item{MSets1.rms.mean}{optional (if more than one MEDIPS SET given): the mean rms value over all MEDIPS SETs at MSet1.}
\item{MSets2.counts.mean}{optional (if more than one MEDIPS SET given): the mean count over all MEDIPS SETs at MSet2.}
\item{MSets2.rpkm.mean}{optional (if more than one MEDIPS SET given): the mean rpkm value over all MEDIPS SETs at MSet2.}
\item{MSets2.rms.mean}{optional (if more than one MEDIPS SET given): the mean rms value over all MEDIPS SETs at MSet2.}
\item{ISets1.counts.mean}{optional (if more than one INPUT SET given): the mean count over all INPUT SETs at ISet1.}
\item{ISets1.rpkm.mean}{optional (if more than one INPUT SET given): the mean rpkm value over all INPUT SETs at ISet1.}
\item{ISets2.counts.mean}{optional (if more than one INPUT SET given): the mean count over all INPUT SETs at ISet2.}
\item{ISets2.rpkm.mean}{optional (if more than one INPUT SET given): the mean rpkm value over all INPUT SETs at ISet2.}
\item{edgeR.logFC}{optional (if diff.method=edgeR): log fold change between MSet1 and MSet2 as returned by edgeR.}
\item{edgeR.logCPM}{optional (if diff.method=edgeR): logCPM between MSet1 and MSet2 as returned by edgeR.}
\item{edgeR.p.value}{optional (if diff.method=edgeR): p.value as returned by edgeR.}
\item{edgeR.adj.p.value}{optional (if diff.method=edgeR): adjusted p.value as calculated by the p.adjust function using edgeR's p.values as input.}
\item{score.log2.ratio}{optional (if diff.method=ttest): log2 fold change between the means of the groups MSet1 and MSet2.}
\item{score.p.value}{optional (if diff.method=ttest): p.value as returned by the t.test function.}
\item{score.adj.p.value}{optional (if diff.method=ttest): adjusted p.value as calculated by the p.adjust function using the ttest p.values as input.}
\item{score}{optional (if diff.method=ttest): score = (-log10(score.p.value)*10)*log(score.log2.ratio)}
\item{CNV.log2.ratio}{optional (if two INPUT SETs given and CNV=TRUE): the log2 ratio for segments as calculated by the DNAcopy package.}
}
\author{
Lukas Chavez, Matthias Lienhard, Joern Dietrich
}
\examples{
library(MEDIPSData)
data(hESCs_MeDIP)
data(DE_MeDIP)
data(hESCs_Input)
data(DE_Input)
data(CS)

resultTable = MEDIPS.meth(MSet1 = hESCs_MeDIP, MSet2 = DE_MeDIP, CSet = CS, ISet1 = hESCs_Input, ISet2 = DE_Input, chr = "chr22", p.adj="bonferroni", diff.method="ttest", CNV=FALSE, MeDIP=FALSE)
}


