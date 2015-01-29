\name{MEDIPS.exportWIG}
\alias{MEDIPS.exportWIG}
\title{
Exports count, rpkm, or sequence pattern densities into a wiggle file.
}
\description{
The function allows for exporting the calculated methylation values (counts or rpkm) or sequence pattern densities from a MEDIPS or COUPLING SET into a wiggle (WIG) file.
The wiggle file will contain values for all genomic windows of the genome/coupling vector and can be used for data visualization using appropriate genome browsers.
Either a MEDIPS SET (parameter MSet) or a COUPLING SET (parameter CSet) has to be given. 
}
\usage{
MEDIPS.exportWIG(Set=NULL, CSet=NULL, file=NULL, format="rpkm", descr="")
}
\arguments{
  \item{Set}{
has to be a MEDIPS SET object. Required when the parameter 'format' is 'count', 'rpkm', or 'rms'.
}
  \item{CSet}{
has to be a COUPLING SET object. Required when the parameter 'format' is 'pdensity' or 'rms'.
}
  \item{file}{
defines the name of the exported file
}
  \item{format}{
if set to "count", there must be a MEDIPS SET at 'Set'. The number of overlapping (extended) short reads per window will be exported.
if set to "rpkm", there must be a MEDIPS SET at 'Set'. The rpkm values will be exported (default). 
If set to "pdensity", there must be a COUPLING SET at 'CSet'. The pattern densities (counts per window) will be exported (parameter Set will be ignored).
If set to 'rms', there must be a MEDIPS SET at 'Set' and a corresponding COUPLING SET at 'CSet'. The CpG density normalized methylatin estimates will be exported.
}
  \item{descr}{
the exported wiggle file will include a track name and description that will be visualized by the utilized genome browser.
Both, track name and description will be the same as defined here.
}  
}
\value{
the funtion exports the specified data from the MEDIPS or COUPLING SET into the stated file
}
\author{
Lukas Chavez
}
\examples{

library("BSgenome.Hsapiens.UCSC.hg19")

bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")
MSet=MEDIPS.createSet(file=bam.file.hESCs.Rep1.MeDIP, BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select="chr22", extend=250, shift=0, uniq=TRUE)

MEDIPS.exportWIG(Set=MSet, file="hESCs.Rep1.wig", format="rpkm", descr="hESCs.Rep1")

}

