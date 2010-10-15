\name{MEDIPS.exportWIG}
\alias{MEDIPS.exportWIG}
\title{
Exports reads per million, relative methylation score or sequence pattern density into a wiggle file.
}
\description{
The funtion allows for exporting the calculated methylation values (rpm or rms) or sequence pattern denisties from a MEDIPS SET into a wiggle (WIG) file.
The wiggle file contains values for all genomic bins of the genome vector and can be used for data visualization using appropiate genome browsers. 
}
\usage{
MEDIPS.exportWIG(data = NULL, file = NULL, raw = FALSE, descr = "", pattern.density = FALSE)
}
\arguments{
  \item{data}{
has to be a MEDIPS SET object
}
  \item{file}{
degines the name of the exported file
}
  \item{raw}{
if set to TRUE, the reads per million values will be exported. 
If set to FALSE, the rms values will be exported (default=FALSE).
}
  \item{descr}{
the exported wiggle file will include a track name and description that will be visualized by the utilized genome browser.
Both, track name and description will be the same as defined here.
}
  \item{pattern.density}{
if set to TRUE, the wiggle file will contain the sequence pattern densities of the coupling vector instead of the methylation values (default=FALSE).
}
}
\value{
the funtion exports the specified data from the MEDIPS SET into the stated file
}
\author{
Lukas Chavez
}
\examples{

library(BSgenome.Hsapiens.UCSC.hg19)
file=system.file("extdata", "MeDIP_hESCs_chr22.txt", package="MEDIPS")
CONTROL.SET = MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file=file)
CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, bin_size = 50, extend = 400)

MEDIPS.exportWIG(file = "example.output.WIG", data = CONTROL.SET, raw = TRUE, descr = "example")

}

