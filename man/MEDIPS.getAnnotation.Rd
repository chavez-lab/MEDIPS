\name{MEDIPS.getAnnotation}
\alias{MEDIPS.getAnnotation}
\title{
Funtion to fetch annotations from biomaRt.
}
\description{
The function receives predifined annotations from ensembl biomaRt for subsequent annotation of a result table.
}
\usage{
MEDIPS.getAnnotation(host="www.ensembl.org",dataset=c("hsapiens_gene_ensembl","mmusculus_gene_ensembl")[1],annotation=c("TSS","EXON","GENE")[1],tssSz=c(-1000,500),chr=NULL)
}
\arguments{
  \item{host}{
BioMart database host you want to connect to.
For current ensembl version, use "www.ensembl.org". For other versions, set to the respective archive host, e.g. "may2012.archive.ensembl.org" for Ensembl 67
Please ensure that the ensembl version is compatible to the used genome version. 
}
  \item{dataset}{
The dataset you want to use.  To see the different datasets available within a biomaRt you can e.g. do: mart = useMart('ensembl'), followed by listDatasets(mart).
}
\item{annotation}{
Type of annotation you want to retrieve.  You can select "EXON" for exonic regions, "GENE" for gene regions including introns and "TSS" for regions at the transcirption start site. 
}
\item{tssSz}{
Defines the TSS region: start and end position relative to the strand and position of the transcript. 
}
\item{chr}{
Chromosome names for which the annotations should be filtered.
}
}
\value{
The MEDIPS.getAnnotation function returns a list of annotation tables where each table consists of
\item{id}{the id of the annotation}
\item{chr}{the chromosome of the annotation}
\item{start}{the start position (5') of the annotation}
\item{end}{the end position (3') of the annotation}
}
\author{
Joern Dietrich and Matthias Lienhard
}

\examples{

#homo sapiens, current ensembl version
#annotation_hs = MEDIPS.getAnnotation(dataset="hsapiens_gene_ensembl", annotation="TSS", chr=c("chr22"),tssSz=c(-1000,500))

#mus musculus, ensembl version 58 (mm9)
annotation_mm9 = MEDIPS.getAnnotation(host="may2010.archive.ensembl.org",dataset="mmusculus_gene_ensembl", annotation="TSS", chr=c("chr19"))

}
