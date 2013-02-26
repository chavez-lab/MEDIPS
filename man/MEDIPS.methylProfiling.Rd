\name{MEDIPS.methylProfiling}
\alias{MEDIPS.methylProfiling}
\title{
Funtion calculates mean methylation values (rpm, rms, ams), ratios, variances, and p-values comparing two MEDIPS SETs for user supplied regions of interests (ROIs) or for genome wide frames.
}
\description{
This function has been deprecated. 
Please see MEDIPS.meth instead.

In order to compare two different conditions, first you have to create and process two sets of MEDIPS SETs.
For the identification of DMRs, MEDIPS provides two alternative approaches. 
First, you can specify pre-defined regions of interest (ROIs).
Second, MEDIPS offers the possibility to calculate differential methylation for genome wide frames.
The function calculates summarized methylation values for the defined ROIs. 
Here, these are the mean values for the provided MEDIPS SETs as well as the ratio of means. 
Moreover, for each ROI, MEDIPS calculates p-values by comparing the set of rpm values (or rms values, respectively) within the ROI of the one MEDIPS SET 
against the set of rpm values (or rms values, respectively) within the ROI of the second MEDIPS SET using R's wilcox.test and t.test functions.
Additionally, it is recommended (but not necessary) to provide background data from an INPUT experiment (that is sequencing of none-enriched DNA fragments).
By providing an INPUT data set, MEDIPS additionally returns mean INPUT rpm values for the specified ROIs.
Becuase the function allows for processing a variable number of provided MEDIPS SETs, the returned matrix is of variable length.
}
\usage{
MEDIPS.methylProfiling(data1 = NULL, data2 = NULL, input = NULL, ROI_file = NULL, frame_size = NULL, math = mean, step = NULL, select = 2, chr = NULL, transf = T)
}
\arguments{
  \item{data1}{
MEDIPS SET control
}
  \item{data2}{
MEDIPS SET treatment
}
  \item{input}{
INPUT SET
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
transforms the resulting data range into the interval [0:1000]
}
}
\value{
\item{Chr}{the chromosome of the ROI}
\item{Start}{the start position of the ROI}
\item{Stop}{the stop position of the ROI}
\item{Items}{the number of genomic bins included in the ROI}
\item{CF}{the mean coupling factor of the ROI}
\item{RPM_MSet1.*}{the mean reads per million value for the MEDIPS SET at position * of MSet1}
\item{RMS_MSet1.*}{the mean relative methylation score for the MEDIPS SET at position * of MSet1}
\item{AMS_MSet1.*}{the mean absolute mathylation score for the MEDIPS SET at position * of MSet1}
\item{RPM_MSet2.*}{the mean reads per million value for the MEDIPS SET at position * of MSet2 (if provided)}
\item{RMS_MSet2.*}{the mean relative methylation score for the MEDIPS SET at position * of MSet2 (if provided)}
\item{AMS_MSet2.*}{the mean absolute mathylation score for the MEDIPS SET at position * of MSet2 (if provided)}
\item{RPM_ISet1.*}{the mean reads per million value of the Input MEDIPS SET at position * of ISet1 (if provided)}
\item{RPM_ISet2.*}{the mean reads per million value of the Input MEDIPS SET at position * of ISet2 (if provided)}
\item{RPM_MSets1}{the mean reads per million value over all MEDIPS SETs at MSet1}
\item{RMS_MSets1}{the mean relative methylation score over all MEDIPS SETs at MSet1}
\item{AMS_MSets1}{the mean absolute methylation score over all MEDIPS SETs at MSet1}
\item{RPM_MSets2}{the mean reads per million value over all MEDIPS SETs at MSet2}
\item{RMS_MSets2}{the mean relative methylation score over all MEDIPS SETs at MSet2}
\item{AMS_MSets2}{the mean absolute methylation score over all MEDIPS SETs at MSet2}
\item{RPM_ISets1}{the mean reads per million value over all MEDIPS SETs at ISet1}
\item{RPM_ISets2}{the mean reads per million value over all MEDIPS SETs at ISet2}
\item{V_MSet1}{the variance of the rpm or rms values (please see the parameter select) over all MEDIPS SETs at MSet1}
\item{CV_MSet1}{the coefficient of variance of the rpm or rms values (please see the parameter select) over all MEDIPS SETs at MSet1}
\item{V_MSet2}{the variance of the rpm or rms values (please see the parameter select) of all MEDIPS SETs at MSet2 (if provided)}
\item{CV_MSet2}{the coefficient of variance of the rpm or rms values (please see the parameter select) over all MEDIPS SETs at MSet2 (if provided)}
\item{Ratio}{RPM_MSets1/RPM_MSets2 or RMS_MSets1/RMS_MSets2, respectively (please see the parameter select)}
\item{Wilcox}{the p.value returned by R's wilcox.test function for comparing the rpm values (or rms values, respectively; please see the parameter select) 
of the MEDIPS SETs at MSet1 and of the MEDIPS SETs at MSet2}
\item{T.test}{the p.value returned by R's t.test function for comparing the rpm values (or rms values, respectively; please see the parameter select) 
of the MEDIPS SETs at MSet1 and of the MEDIPS SETs at MSet2}
\item{Wilcox.adj}{multiple testing corrected p.value from slot Wilcox; calculated by R's p.adjust function}
\item{T.test.adj}{multiple testing corrected p.value from slot T.Test; calculated by R's p.adjust function}
}
\author{
Joern Dietrich
}
\examples{

print("This function has been deprecated. Please see MEDIPS.meth instead.")

}


