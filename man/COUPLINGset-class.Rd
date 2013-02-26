\name{COUPLINGset-class}
\Rdversion{1.1}
\docType{class}
\alias{COUPLINGset}
\alias{COUPLINGset-class}
\alias{genome_name,COUPLINGset-method}
\alias{window_size,COUPLINGset-method}
\alias{chr_names,COUPLINGset-method}
\alias{chr_lengths,COUPLINGset-method}
\alias{seq_pattern,COUPLINGset-method}
\alias{genome_CF,COUPLINGset-method}
\alias{number_pattern,COUPLINGset-method}
\alias{show,COUPLINGset-method}

\alias{seq_pattern}
\alias{genome_CF}
\alias{number_pattern}

\title{COUPLINGset class and internal functions}
\description{COUPLINGset class is used in the MEDIPS library to store and extract information generated during the creation of a coupling vector.}
\section{Objects from the Class}{
Objects of the classes contain information about sequence pattern information, included chromosomes, and further parameter settings.
A COUPLING SET object is created by the MEDIPS.couplingVector() function.
According slots will be filled during the workflow.
}
\section{Slots}{
  \describe{
	\item{\code{genome_name}:}{Object of class \code{"character"} : the refernce genome }
	\item{\code{window_size}:}{Object of class \code{"numeric"} : the window size for the genome vector }
	\item{\code{chr_names}:}{Object of class \code{"character"} : the names of the chromosomes included within the MEDIPS/COUPLING SET }
	\item{\code{chr_lengths}:}{Object of class \code{"numeric"} : the lengths of the chromosomes included within the MEDIPS/COUPLING SET}
	\item{\code{seq_pattern}:}{Object of class \code{"character"} : the sequence pattern (e.g. CG) }  
	\item{\code{genome_CF}:}{Object of class \code{"numeric"} : the coupling factor at the genomic bins }
	\item{\code{number_pattern}:}{Object of class \code{"numeric"} : the total number of sequence pattern }
  }
}
\section{Methods}{
  \describe{	 
    \item{genome_name}{\code{signature(object = "COUPLINGset")}: extracts the reference genome of the COUPLING SET}
    \item{window_size}{\code{signature(object = "COUPLINGset")}: extracts the window size from the window_size slot COUPLING SET}
    \item{chr_names}{\code{signature(object = "COUPLINGset")}: extracts the names of the chromosomes included within the COUPLING SET }
    \item{chr_lengths}{\code{signature(object = "COUPLINGset")}: extracts the length of the chromosomes included within the COUPLING SET }
    \item{seq_pattern}{\code{signature(object = "COUPLINGset")}: extracts the sequence pattern (e.g. CpG)}
    \item{genome_CF}{\code{signature(object = "COUPLINGset")}: extracts the coupling factor at the genomic bins }
    \item{number_pattern}{\code{signature(object = "COUPLINGset")}: extracts the total number of sequence pattern}
    \item{show}{\code{signature(object = "COUPLINGset")}: prints a summary of the COUPLING SET object content }       
    }
}
\author{Lukas Chavez, Matthias Lienhard, Joern Dietrich}
\examples{
showClass("COUPLINGset")
}
\keyword{classes}
