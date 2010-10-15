\name{MEDIPS.selectSignificants}
\alias{MEDIPS.selectSignificants}
\title{
Selects candidate ROIs that show significant differential methylation between two MEDIPS SETs.
}
\description{
Based on the results matrix returned from the MEDIPS.diffMethyl function,
the function selects candidate ROIs that show significant differential methylation between the CONTROL.SET and the TREAT.SET in consideration of the background data included in the INPUT.SET.
Filtering for significant frames proceeds in the following order:
ROIs that do not contain any data either in the CONTROL.SET nor in the TREAT.SET are neglected first;
ROIs associated to p-values > p.value are neglected;
ROIs with a CONTROL/TREATMENT ratio < up (or > down, respectively) are neglected;
From the INPUT mean rpm distribution, a mean rpm threshold was defined by the quant parameter and
all ROIs that have a mean rpm value within the CONTROL.SET (or TREAT.SET, respectively) smaller than the estimated background rpm threshold are discarded;
The last filter is again based on the INPUT data. While the latter filter estimates a minimum rpm signal for the CONTROL.SET (or TREAT.SET, respectively) from the total background distribution,
we now define that the rpm value from the CONTROL SET (or TREAT.SET, respectively) of a ROI exceeds the local background data of the INPUT.SET by the parameter up.
This is, because MeDIP-Seq background data varies along the chromosomes due to varying DNA availability. 
}
\usage{
MEDIPS.selectSignificants(frames = NULL, input = T, control = T, up = 1.333333, down = 0.75, p.value = 0.01,quant = 0.9)
}
\arguments{
  \item{frames}{
specifies the results table derived from the MEDIPS.diffMethyl
}
  \item{input}{
default=T; Setting the parameter to TRUE requires that the results table includes a column for summarized rpm values of an INPUT SET.
In case, there is no INPUT data available, the input parameter has to be set to a rpm value that will be used as threshold during the subsequent analysis.
How to estimate such a threshold without background data is not yet solved by MEDIPS.
}
  \item{control}{
can be either TRUE or FALSE;
MEDIPS allows for selecting frames that are higher methylated in the CONTROL SET compared to the TREAT SET and vice versa but both approaches have to be perfomed in two independent runs.
By setting control=T, MEDIPS selects genomic regions, where the CONTROL SET is higher methylated.
By setting control=F, MEDIPS selects genomic regions, where the TREAT SET is higher methylated.
}
\item{up}{
default=1.333333; defines the lower threshold for the ratio CONTROL/TREAT as well as for the lower ratio for CONTROL/INPUT (if control=T) or TREATMENT/INPUT (if control=F), respectively.
}
  \item{down}{
default=0.75; defines the upper threshold for the ratio: CONTROL/TREATMENT (only if control=F).
}
  \item{p.value}{
default=0.01; defines the threshold for the p-values. 
One of the p-values derived from the wilcox.test or t.test function has to be <= p.value. 
}
  \item{quant}{
default=0.9; from the distribution of all summarized INPUT rpm values, 
MEDIPS calculates the rpm value that represents the quant quantile of the whole INPUT distribution.}
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
Lukas Chavez
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

file=system.file("extdata", "MeDIP_DE_chr22.txt", package="MEDIPS")
TREAT.SET = MEDIPS.readAlignedSequences(BSgenome = "BSgenome.Hsapiens.UCSC.hg19", file = file)
TREAT.SET = MEDIPS.genomeVector(data = TREAT.SET, bin_size = 50, extend = 400)
TREAT.SET = MEDIPS.getPositions(data = TREAT.SET, pattern = "CG")
TREAT.SET = MEDIPS.couplingVector(data = TREAT.SET, fragmentLength = 700, func = "count")
TREAT.SET = MEDIPS.calibrationCurve(data = TREAT.SET)
TREAT.SET = MEDIPS.normalize(data = TREAT.SET)

file=system.file("extdata", "Input_StemCells_chr22.txt", package="MEDIPS")
INPUT.SET = MEDIPS.readAlignedSequences(BSgenome = "BSgenome.Hsapiens.UCSC.hg19", file = file)
INPUT.SET = MEDIPS.genomeVector(data = INPUT.SET, bin_size = 50, extend = 400)

diff.methyl = MEDIPS.methylProfiling(data1 = CONTROL.SET, data2= TREAT.SET, input=INPUT.SET, chr="chr22", frame_size=1000, select=1)

diff.methyl.sig=MEDIPS.selectSignificants(diff.methyl)


}


