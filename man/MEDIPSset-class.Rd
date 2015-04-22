\name{MEDIPSset-class}
\Rdversion{1.1}
\docType{class}
\alias{MEDIPSset-class}
\alias{MEDIPSset}
\alias{genome_name,MEDIPSset-method}
\alias{bin_size,MEDIPSset-method}
\alias{chr_names,MEDIPSset-method}
\alias{chr_lengths,MEDIPSset-method}
\alias{sample_name,MEDIPSset-method}
\alias{path_name,MEDIPSset-method}
\alias{number_regions,MEDIPSset-method}
\alias{genome_count,MEDIPSset-method}
\alias{extend,MEDIPSset-method}
\alias{show,MEDIPSset-method}
\alias{shifted,MEDIPSset-method}
\alias{window_size,MEDIPSset-method}
\alias{uniq,MEDIPSset-method}

\alias{setTypes}
\alias{path_name}
\alias{getTypes}
\alias{getGRange}
\alias{adjustReads}
\alias{MEDIPS.Bed2Granges}
\alias{MEDIPS.Bam2GRanges}

\alias{genome_name}
\alias{bin_size}
\alias{chr_names}
\alias{chr_lengths}
\alias{sample_name}
\alias{number_regions}
\alias{genome_count}
\alias{extend}
\alias{show}
\alias{shifted}
\alias{uniq}
\alias{window_size}

\alias{MEDIPS.GenomicCoordinates}
\alias{MEDIPS.calibrationCurve}
\alias{MEDIPS.cnv}
\alias{MEDIPS.getPositions}
\alias{MEDIPS.rms}

\alias{matNnotNA}
\alias{matMin}
\alias{matDiff}
\alias{matMax}
\alias{matMean}
\alias{matSd}
\alias{matTtest}



\title{MEDIPSset class and internal functions}
\description{MEDIPSset class is used in the MEDIPS library in order to store and extract objects and information from the input file as well as parameter settings specified during the workflow.}
\section{Objects from the Class}{
Objects of the classes contain information about the provided short reads, MeDIP raw/count signals, and further parameter settings.
A MEDIPS SET object is created by the MEDIPS.genomeVector() function.
According slots will be filled during the workflow.
}
\section{Slots}{
  \describe{
	\item{\code{genome_name}:}{Object of class \code{"character"} : the refernce genome }
	\item{\code{window_size}:}{Object of class \code{"numeric"} : the window size for the genome vector }
	\item{\code{chr_names}:}{Object of class \code{"character"} : the names of the chromosomes included within the MEDIPS/COUPLING SET }
	\item{\code{chr_lengths}:}{Object of class \code{"numeric"} : the lengths of the chromosomes included within the MEDIPS/COUPLING SET}
	\item{\code{sample_name}:}{Object of class \code{"character"} : the name of the input file }
	\item{\code{path_name}:}{Object of class \code{"character"} : the path to the input file }
	\item{\code{number_regions}:}{Object of class \code{"numeric"} : the total number of included regions }
	\item{\code{genome_count}:}{Object of class \code{"numeric"} : the raw MeDIP-seq signals at the genomic bins }    
	\item{\code{extend}:}{Object of class \code{"numeric"} : the length of the regions after extension }
	\item{\code{shifted}:}{Object of class \code{"numeric"} : the number of bases by which the reads are shifted along the sequencing direction }
	\item{\code{uniq}:}{Object of class \code{"logical"} : determines if all reads mapping to exactly the same genomic position should be kept (uniq = 0), replaced by only one representative (uniq = 1), or if the number of stacked reads should be capped by a maximal number of stacked reads per genomic position determined by a poisson distribution of stacked reads genome wide and by a given p-value (1 > uniq > 0).}                
  }
}
\section{Methods}{
  \describe{
    \item{genome_name}{\code{signature(object = "MEDIPSset")}: extracts the reference genome of the MEDIPS SET}    
    \item{window_size}{\code{signature(object = "MEDIPSset")}: extracts the window size from the bin_size slot of the MEDIPS SET}
    \item{chr_names}{\code{signature(object = "MEDIPSset")}: extracts the names of the chromosomes included within the MEDIPS SET }
    \item{chr_lengths}{\code{signature(object = "MEDIPSset")}: extracts the length of the chromosomes included within the MEDIPS SET }
    \item{fragmentLength}{\code{signature(object = "MEDIPSset")}: extracts the estimated fragment length of the DNA fragments }
    \item{sample_name}{\code{signature(object = "MEDIPSset")}: extracts the name of the input file}
    \item{path_name}{\code{signature(object = "MEDIPSset")}: extracts the path to the input file}
    \item{number_regions}{\code{signature(object = "MEDIPSset")}: extracts the total number of included regions}
    \item{genome_count}{\code{signature(object = "MEDIPSset")}: extracts the raw MeDIP-Seq signals at the genomic bins }    
    \item{extend}{\code{signature(object = "MEDIPSset")}: extracts the number of bases by which the regions are extended }
    \item{show}{\code{signature(object = "MEDIPSset")}: prints a summary of the MEDIPS SET object content }
    \item{shifted}{\code{signature(object = "MEDIPSset")}: extracts the number of bases by which the regions are shifted }
    \item{uniq}{\code{signature(object = "MEDIPSset")}: extracts the specified value for the uniq parameter }
        
    \item{MEDIPS.GenomicCoordinates}{\code{signature(object = "MEDIPSset")}: internal function for calculating coordinates for the genomic bins }    
    \item{MEDIPS.readRegionsFile}{\code{signature(object = "MEDIPSset")}: internal function for reading short read information }
    \item{MEDIPS.calibrationCurve}{\code{signature(object = "MEDIPSset")}: internal function for calculating the calibration curve}
    \item{MEDIPS.cnv}{\code{signature(object = "MEDIPSset")}: internal function for calculating CNVs in case two groups of INPUT SETs have been provided to the MEDIPS.meth function}
    \item{MEDIPS.diffMeth}{\code{signature(object = "MEDIPSset")}: internal function for calculating differential coverage in case two groups of MEDIPS SETs have been provided to the MEDIPS.meth function}
    \item{MEDIPS.getPositions}{\code{signature(object = "MEDIPSset")}: internal function for receiving genomic coordinates of a given sequence pattern (e.g. CG)}
    \item{MEDIPS.rms}{\code{signature(object = "MEDIPSset")}: internal function for calculating relative methylation scores}

    \item{matNnotNA}{\code{signature(object = "MEDIPSset")}: internal function for vectorized calculation of the t.test }
    \item{matMin}{\code{signature(object = "MEDIPSset")}:   internal function for vectorized calculation of the t.test }
    \item{matDiff}{\code{signature(object = "MEDIPSset")}:  internal function for vectorized calculation of the t.test }
    \item{matMax}{\code{signature(object = "MEDIPSset")}:  internal function for vectorized calculation of the t.test }
    \item{matMean}{\code{signature(object = "MEDIPSset")}:  internal function for vectorized calculation of the t.test }
    \item{matSd}{\code{signature(object = "MEDIPSset")}:  internal function for vectorized calculation of the t.test }
    \item{matTtest}{\code{signature(object = "MEDIPSset")}:  internal function for vectorized calculation of the t.test  }
   
    
    }
}
\author{Lukas Chavez, Joern Dietrich }
\examples{
showClass("MEDIPSset")
}
\keyword{classes}
