\name{MEDIPS.methylProfiling}
\alias{MEDIPS.methylProfiling}
\title{
Funtion calculates mean methylation values (rpm, rms) and ams values, ratios, variances, and pvalues comparing two MEDIPS SETs for user supplied regions of interests (ROIs) or genome wide frames.
}
\description{
In order to compare two different conditions, first you have to create and process two MEDIPS SETs.
For the identification of DMRs, MEDIPS provides two alternative approaches. 
First, you can specify pre-defined regions of interest (ROIs).
Second, MEDIPS offers the possibility to calculate differential methylation for genome wide frames.
The function calculates summarized methylation values for the defined ROIs. 
Here, these are the mean values for both provided MEDIPS SETs as well as the ratio of means. 
Moreover, for each ROI, MEDIPS calculates p-values by comparing the set of rpm values (or rms values, respectively) within the ROI of the one MEDIPS SET 
against the set of rpm values (or rms values, respectively) within the ROI of the second MEDIPS SET using R's wilcox.test and t.test functions.
Additionally, it is recommended (but not necessary) to provide background data from an INPUT experiment (that is sequencing of none-enriched DNA fragments).
By providing an INPUT data set, MEDIPS additionally returns mean INPUT rpm values for the specified ROIs.
Please note, the function takes a long processing time when called for genome wide short windows (up to days).
}
\usage{
MEDIPS.methylProfiling(data1 = NULL, data2 = NULL, input = NULL, ROI_file = NULL, frame_size = NULL, math = "mean", step = NULL, select = 2, chr = NULL, transf=T)
}
\arguments{
  \item{data1}{
has to be a MEDIPS SET object (the control data)
}
  \item{data2}{
has to be a MEDIPS SET object (the treatment data)
}
  \item{input}{
has to be a MEDIPS SET object (the input data)
}
  \item{ROI_file}{
instead of processing genome wide frames using the parameters frame_size and step, 
here you can provide a file containing predefined ROIs.
}
  \item{frame_size}{
Besides summarizing methylation values for pre-defined ROIs, MEDIPS allows for calculating mean methylation values along the full chromosomes.
For this, you have to specify a desired frame size here.
}
  \item{math}{
default=mean; Here, you can specify other functions available in R for sumarizing values like median or sum. 
}
  \item{step}{
The step parameter defines the number of bases by which the frames are shifted along the chromosome. 
If you e.g. set the frame_size parameter to 500 and the step parameter to 250, then MEDIPS calculates mean methylation values for overlapping 500bp windows, 
where the size of the overlap will be 250bp for all neighbouring windows. 
}
  \item{select}{
can be either 1 or 2. If set to 1, the variance, ratio, and p-values will be calculated based on the rpm values; if set to 2, the rms values will be considered instead.
}
  \item{chr}{
only the specified chromosome will be evaluated (e.g. chr1)
}
  \item{transf}{
If set to TRUE, MEDIPS transforms the mean rms and ams values into log2 scale and subsequently transforms their resulting data range into the consistent interval $[0,1000]$ before finally stored.
}
}
\value{
\item{chr}{the chromosome of the ROI}
\item{start}{the start position of the ROI}
\item{stop}{the stop position of the ROI}
\item{length}{the number of genomic bins included in the ROI}
\item{coupling}{the mean coupling factor of the ROI}
\item{input}{the mean reads per million value of the INPUT MEDIPS SET at input (if provided)}
\item{rpm_A}{the mean reads per million value for the MEDIPS SET at data1}
\item{rpm_B}{the mean reads per million value for the MEDIPS SET at data2}
\item{rms_A}{the mean relative mathylation score for the MEDIPS SET at data1}
\item{rms_B}{the mean relative methylation score for the MEDIPS SET at data2}
\item{ams_A}{the mean absolute mathylation score for the MEDIPS SET at data1.
The ams scores are derived by dividing the mean rms value of the ROI by the mean coupling factor of the ROI before the log2 and interval transformations are performed.}
\item{ams_B}{the mean absolute mathylation score for the MEDIPS SET at data2.
The ams scores are derived by dividing the mean rms value of the ROI by the mean coupling factor of the ROI before the log2 and interval transformations are performed.}
\item{var_A}{the variance of the rpm or rms values (please see the parameter select) of the MEDIPS SET at data1}
\item{var_B}{the variance of the rpm or rms values (please see the parameter select) of the MEDIPS SET at data2}
\item{var_co_A}{the variance coefficient of the rpm or rms values (please see the parameter select) of the MEDIPS SET at data1}
\item{var_co_B}{the variance coefficient of the rpm or rms values (please see the parameter select) of the MEDIPS SET at data2}
\item{ratio}{rpm_A/rpm_B or rms_A/rms_B, respectively (please see the parameter select)}
\item{pvalue.wilcox}{the p.value returned by R's wilcox.test function for comparing the rpm values (or rms values, respectively; please see the parameter select) 
of the MEDIPS SET at data1 and of the MEDIPS SET at data2}
\item{pvalue.ttest}{the p.value returned by R's t.test function for comparing the rpm values (or rms values, respectively; please see the parameter select) 
of the MEDIPS SET at data1 and of the MEDIPS SET at data2}
}
\author{
Joern Dietrich
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)
CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, bin_size = 50, extend = 400)
CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")
CONTROL.SET = MEDIPS.couplingVector(data = CONTROL.SET, fragmentLength = 700, func = "count")
CONTROL.SET = MEDIPS.calibrationCurve(data = CONTROL.SET)
CONTROL.SET = MEDIPS.normalize(data = CONTROL.SET)
ROI_file = system.file("extdata", "hg19.chr22.txt", package="MEDIPS")

promoter = MEDIPS.methylProfiling(data1 = CONTROL.SET, ROI_file = ROI_file, math = mean, select = 2)
}


