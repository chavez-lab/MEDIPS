\name{MEDIPSset-class}
\Rdversion{1.1}
\docType{class}
\alias{MEDIPSset-class}
\alias{bin_size,MEDIPSset-method}
\alias{calcurve_mean_coupling,MEDIPSset-method}
\alias{calcurve_mean_signals,MEDIPSset-method}
\alias{calcurve_var,MEDIPSset-method}
\alias{cali_chr,MEDIPSset-method}
\alias{chr_lengths,MEDIPSset-method}
\alias{chr_names,MEDIPSset-method}
\alias{distFile,MEDIPSset-method}
\alias{distFunction,MEDIPSset-method}
\alias{extend,MEDIPSset-method}
\alias{fragmentLength,MEDIPSset-method}
\alias{genome_CF,MEDIPSset-method}
\alias{genome_chr,MEDIPSset-method}
\alias{genome_name,MEDIPSset-method}
\alias{genome_norm,MEDIPSset-method}
\alias{genome_pos,MEDIPSset-method}
\alias{genome_raw,MEDIPSset-method}
\alias{intercept,MEDIPSset-method}
\alias{number_pattern,MEDIPSset-method}
\alias{number_regions,MEDIPSset-method}
\alias{pattern_chr,MEDIPSset-method}
\alias{pattern_pos,MEDIPSset-method}
\alias{regions_chr,MEDIPSset-method}
\alias{regions_start,MEDIPSset-method}
\alias{regions_stop,MEDIPSset-method}
\alias{regions_strand,MEDIPSset-method}
\alias{sample_name,MEDIPSset-method}
\alias{seq_pattern,MEDIPSset-method}
\alias{show,MEDIPSset-method}
\alias{slope,MEDIPSset-method}
\alias{bin_size}
\alias{calcurve_mean_coupling}
\alias{calcurve_mean_signals}
\alias{calcurve_var}
\alias{cali_chr}
\alias{chr_lengths}
\alias{chr_names}
\alias{distFile}
\alias{distFunction}
\alias{extend}
\alias{fragmentLength}
\alias{genome_CF}
\alias{genome_chr}
\alias{genome_name}
\alias{genome_norm}
\alias{genome_pos}
\alias{genome_raw}
\alias{intercept}
\alias{number_pattern}
\alias{number_regions}
\alias{pattern_chr}
\alias{pattern_pos}
\alias{regions_chr}
\alias{regions_start}
\alias{regions_stop}
\alias{regions_strand}
\alias{sample_name}
\alias{seq_pattern}
\alias{show}
\alias{slope}
\alias{MEDIPS.distributeReads}
\alias{MEDIPS.transform}

\title{Class "MEDIPSset" }
\description{ This class is used in the MEDIPS library to store and extract all the MEDIPS objects and inormations generated during the workflow}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MEDIPSset", ...)}.
	 Objects of the class contain the information about the provided regions, their chromosome, start, stop and strand informations, 
	 the raw and normalized signals, the sequence pattern positions and informations about all specified parameters.
	 The MEDIPS SET will be created by reading the input regions file using the MEDIPS.readAlignedSequences function and all further slots will be filled during the workflow.
}
\section{Slots}{
  \describe{
    \item{\code{regions_chr}:}{Object of class \code{"character"} : the chromosomes of the regions }
    \item{\code{regions_start}:}{Object of class \code{"numeric"} : the start positions of the regions }
    \item{\code{regions_stop}:}{Object of class \code{"numeric"} : the stop positions of the regions }
    \item{\code{regions_strand}:}{Object of class \code{"character"} : the strand informations of the regions }
    \item{\code{number_regions}:}{Object of class \code{"numeric"} : the total number of included regions }
    \item{\code{pattern_chr}:}{Object of class \code{"character"} : the chromosomes of the sequence pattern }
    \item{\code{pattern_pos}:}{Object of class \code{"numeric"} : the positions of the sequence pattern }
    \item{\code{number_pattern}:}{Object of class \code{"numeric"} : the total number of sequence pattern }
    \item{\code{genome_chr}:}{Object of class \code{"character"} : the chromosome of the genomic bins }
    \item{\code{genome_pos}:}{Object of class \code{"numeric"} : the positions of the genomic bins }
    \item{\code{genome_CF}:}{Object of class \code{"numeric"} : the coupling factor at the genomic bins }
    \item{\code{genome_raw}:}{Object of class \code{"numeric"} : the raw MeDIP-Seq signals at the genomic bins }
    \item{\code{genome_norm}:}{Object of class \code{"numeric"} : the rms values at the genomic bins }
    \item{\code{genome_name}:}{Object of class \code{"character"} : the refernce genome }
    \item{\code{bin_size}:}{Object of class \code{"numeric"} : the bin size for the genome vector }
    \item{\code{extend}:}{Object of class \code{"numeric"} : the number of bases by which the regions are extended }
    \item{\code{fragmentLength}:}{Object of class \code{"numeric"} : the estimated fragment length of the DNA fragments }
    \item{\code{sample_name}:}{Object of class \code{"character"} : the name of the input file }
    \item{\code{chr_lengths}:}{Object of class \code{"numeric"} : the lengths of the chromosomes included within the MEDIPS SET }
    \item{\code{chr_names}:}{Object of class \code{"character"} : the names of the chromosomes included within the MEDIPS SET }
    \item{\code{seq_pattern}:}{Object of class \code{"character"} : the sequence pattern (e.g. CpG) }
    \item{\code{distFunction}:}{Object of class \code{"character"} : the distance function specified for calculating the coupling factors }
    \item{\code{distFile}:}{Object of class \code{"character"} : the file containing the custom distance weights, if provided }
    \item{\code{calcurve_mean_signals}:}{Object of class \code{"numeric"} : the mean signals of the calibration curve }
    \item{\code{calcurve_mean_coupling}:}{Object of class \code{"numeric"} : the mean coupling factors of the calibration curve }
    \item{\code{calcurve_var}:}{Object of class \code{"numeric"} : the signal variance of the levels of the calibration curve }
    \item{\code{intercept}:}{Object of class \code{"numeric"} : the intercept value calculated from the calibration curve }
    \item{\code{slope}:}{Object of class \code{"numeric"} : the slope value calculated from the calibration curve }
    \item{\code{cali_chr}:}{Object of class \code{"character"} : the chromosomes used for calibration }
  }
}
\section{Methods}{
  \describe{
    \item{bin_size}{\code{signature(object = "MEDIPSset")}: extracts the bins size from the bin_size slot }
    \item{calcurve_mean_coupling}{\code{signature(object = "MEDIPSset")}: extracts the mean coupling factors of the calibration curve }
    \item{calcurve_mean_signals}{\code{signature(object = "MEDIPSset")}: extracts the mean signals of the calibration curve }
    \item{calcurve_var}{\code{signature(object = "MEDIPSset")}: extracts the variances of the calibration curve }
    \item{cali_chr}{\code{signature(object = "MEDIPSset")}: extracts the chromosome used for calculating the calibration curve }
    \item{chr_lengths}{\code{signature(object = "MEDIPSset")}: extracts the length of the chromosomes included within the MEDIPS SET }
    \item{chr_names}{\code{signature(object = "MEDIPSset")}: extracts the names of the chromosomes included within the MEDIPS SET }
    \item{distFile}{\code{signature(object = "MEDIPSset")}: extracts the name of file containing the custom distance weights, if provided }
    \item{distFunction}{\code{signature(object = "MEDIPSset")}: extracts the distance function specified for calculating the coupling factors }
    \item{extend}{\code{signature(object = "MEDIPSset")}: extracts the number of bases by which the regions are extended }
    \item{fragmentLength}{\code{signature(object = "MEDIPSset")}: extracts the estimated fragment length of the DNA fragments }
    \item{genome_CF}{\code{signature(object = "MEDIPSset")}: extracts the coupling factor at the genomic bins }
    \item{genome_chr}{\code{signature(object = "MEDIPSset")}: extracts the chromosome of the genomic bins }
    \item{genome_name}{\code{signature(object = "MEDIPSset")}: extracts the refernce genome }
    \item{genome_norm}{\code{signature(object = "MEDIPSset")}: extracts the rms values at the genomic bins }
    \item{genome_pos}{\code{signature(object = "MEDIPSset")}: extracts the positions of the genomic bins }
    \item{genome_raw}{\code{signature(object = "MEDIPSset")}: extracts the raw MeDIP-Seq signals at the genomic bins }
    \item{intercept}{\code{signature(object = "MEDIPSset")}: extracts the intercept value calculated from the calibration curve }
    \item{number_pattern}{\code{signature(object = "MEDIPSset")}: extracts the total number of sequence pattern}
    \item{number_regions}{\code{signature(object = "MEDIPSset")}: extracts the total number of included regions}
    \item{pattern_chr}{\code{signature(object = "MEDIPSset")}: extracts the chromosomes of the sequence pattern }
    \item{pattern_pos}{\code{signature(object = "MEDIPSset")}: extracts the positions of the sequence pattern}
    \item{regions_chr}{\code{signature(object = "MEDIPSset")}: extracts the chromosomes of the regions}
    \item{regions_start}{\code{signature(object = "MEDIPSset")}: extracts the start positions of the regions}
    \item{regions_stop}{\code{signature(object = "MEDIPSset")}: extracts the stop positions of the regions}
    \item{regions_strand}{\code{signature(object = "MEDIPSset")}: extracts the strand informations of the regions}
    \item{sample_name}{\code{signature(object = "MEDIPSset")}: extracts the name of the input file}
    \item{seq_pattern}{\code{signature(object = "MEDIPSset")}: extracts the sequence pattern (e.g. CpG)}
    \item{show}{\code{signature(object = "MEDIPSset")}: prints a summary of the object content }
    \item{slope}{\code{signature(object = "MEDIPSset")}: extracts the slope value calculated from the calibration curve}
    \item{MEDIPS.distributeReads}{\code{signature(object = "MEDIPSset")}: help function for distributing the reads over the genome vector }
    \item{MEDIPS.transform}{\code{signature(object = "MEDIPSset")}: help function for transforming ams values }
    	 }
}
\author{Lukas Chavez, Joern Dietrich }
\examples{
showClass("MEDIPSset")
}
\keyword{classes}
