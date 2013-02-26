\name{MEDIPSroiSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MEDIPSroiSet-class}
\alias{MEDIPSroiSet}
\alias{rois}
\alias{bin_number}
\alias{ROI,MEDIPSroiSet-method}
\alias{bin_number,MEDIPSroiSet-method}
\alias{genome_name,MEDIPSroiSet-method}
\alias{chr_names,MEDIPSroiSet-method}
\alias{chr_lengths,MEDIPSroiSet-method}
\alias{sample_name,MEDIPSroiSet-method}
\alias{path_name,MEDIPSroiSet-method}
\alias{number_regions,MEDIPSroiSet-method}
\alias{genome_count,MEDIPSroiSet-method}
\alias{extend,MEDIPSroiSet-method}
\alias{show,MEDIPSroiSet-method}
\alias{shifted,MEDIPSroiSet-method}
\alias{uniq,MEDIPSroiSet-method}
\alias{rois,MEDIPSroiSet-method}
\alias{rois}

\title{MEDIPSroiSet class and internal functions}
\description{MEDIPSroiSet class is used in the MEDIPS library in order to store and extract objects and information of the specified regions of interest (ROI) from the input file as well as parameter settings specified during the workflow.}
\section{Objects from the Class}{
Objects of the classes contain information about the provided short reads, MeDIP raw/count signals, and further parameter settings.
A MEDIPS ROI SET object is created by the MEDIPS.createROIset() function.
According slots will be filled during the workflow.
}
\section{Slots}{
  \describe{
	\item{\code{genome_name}:}{Object of class \code{"character"} : the refernce genome }
	\item{\code{chr_names}:}{Object of class \code{"character"} : the names of the chromosomes included within the MEDIPS ROI SET }
	\item{\code{chr_lengths}:}{Object of class \code{"numeric"} : the lengths of the chromosomes included within the MEDIPS ROI SET}
	\item{\code{sample_name}:}{Object of class \code{"character"} : the name of the input file }
	\item{\code{path_name}:}{Object of class \code{"character"} : the path to the input file }
	\item{\code{number_regions}:}{Object of class \code{"numeric"} : the total number of included regions }
	\item{\code{genome_count}:}{Object of class \code{"numeric"} : the raw MeDIP-seq signals at the bins }    
	\item{\code{extend}:}{Object of class \code{"numeric"} : the length of the reads after extension }
	\item{\code{shifted}:}{Object of class \code{"numeric"} : the number of bases by which the reads are shifted along the sequencing direction }
	\item{\code{uniq}:}{Object of class \code{"logical"} : determines if reads mapping to exactly the same genomic position should be replaced by only on representative }
	\item{\code{ROI}:}{Object of class \code{"GRanges"}: the genomic positions of the regions of interest}                
	\item{\code{bin_number}:}{Object of class \code{"numeric"}: the number of bins per region}                

  }
}
\section{Methods}{
  \describe{
    \item{genome_name}{\code{signature(object = "MEDIPSroiSet")}: extracts the reference genome of the MEDIPS ROI SET}    
    \item{bin_number}{\code{signature(object = "MEDIPSroiSet")}: extracts the number of bins per ROI the bin_number slot of the MEDIPS ROI SET}
    \item{chr_names}{\code{signature(object = "MEDIPSroiSet")}: extracts the names of the chromosomes included within the MEDIPS ROI SET }
    \item{chr_lengths}{\code{signature(object = "MEDIPSroiSet")}: extracts the length of the chromosomes included within the MEDIPS ROI SET }
    \item{sample_name}{\code{signature(object = "MEDIPSroiSet")}: extracts the name of the input file}
    \item{path_name}{\code{signature(object = "MEDIPSroiSet")}: extracts the path to the input file}
    \item{number_regions}{\code{signature(object = "MEDIPSroiSet")}: extracts the total number of included regions}
    \item{genome_count}{\code{signature(object = "MEDIPSroiSet")}: extracts the raw MeDIP-Seq signals at the genomic bins }    
    \item{extend}{\code{signature(object = "MEDIPSroiSet")}: extracts the number of bases by which the regions are extended }
    \item{show}{\code{signature(object = "MEDIPSroiSet")}: prints a summary of the MEDIPS SET object content }
    \item{shifted}{\code{signature(object = "MEDIPSroiSet")}: extracts the number of bases by which the regions are shifted }
    \item{uniq}{\code{signature(object = "MEDIPSroiSet")}: extracts the specified value for the uniq parameter }
    \item{rois}{\code{signature(object = "MEDIPSroiSet")}: extracts the GRange object containing the regions of interest }       
    \item{MEDIPS.calibrationCurve}{\code{signature(MSet = "MEDIPSroiSet", CSet="COUPLINGset")}: internal function for calculating the calibration curve}
    \item{MEDIPS.negBin}{\code{signature(MSet="MEDIPSroiSet", CSet="COUPLINGset")}: internal function for calculating methylatiopn probabilities with respect to CpG density dependent negative binomial distributions}
    \item{MEDIPS.pois}{\code{signature(MSet="MEDIPSroiSet", CSet="COUPLINGset")}: internal function for calculating methylatiopn probabilities with respect to CpG density dependent poisson distributions}
    \item{MEDIPS.rms}{\code{signature(MSet="MEDIPSroiSet", CSet="COUPLINGset")}: internal function for calculating relative methylation scores}
    
    }
}
\author{Lukas Chavez, Joern Dietrich }
\examples{
showClass("MEDIPSroiSet")
}
\keyword{classes}
