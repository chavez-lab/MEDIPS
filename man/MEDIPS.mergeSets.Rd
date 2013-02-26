\name{MEDIPS.mergeSets}
\alias{MEDIPS.mergeSets}
\title{
Creates one merged MEDIPS SET out of two.
}
\description{
A MEDIPS SET contains a genome vector which is the count coverage at genome wide windows.
Moreover, the MEDIPS SET stores the total number of reads given for calculating the genome vector.
Two MEDIPS SETs can be merged whenever they have been construcetd based on the same reference genome, the same set of chromosomes and for the same window size.
The returned MEDIPS SET will contain a genome vector where at each window the counts of both given MEDIPS SETs are added. 
In addtion, the total number of reads will be the sum of both given MEDIPS SETs.
Please note, several other attributes like the extend or shift value can be different in both of the given MEDIPS SETs and will be empty in the merged MEDIPS SET.
The merged MEDIPS SET will not contain any path to a concrete input file anymore and therefore, cannot be used for the MEDIPS.addCNV function anymore. 
}
\usage{
MEDIPS.mergeSets(MSet1=NULL, MSet2=NULL, name="Merged Set")
}
\arguments{
  \item{MSet1}{
  A MEDIPS SET object as created by the MEDIPS.createSet function
}
\item{MSet2}{
  A MEDIPS SET object as created by the MEDIPS.createSet function
}
\item{name}{
  The new sample name of the merged MEDIPS SET
}
}
\author{
Lukas Chavez
}
\examples{

library(MEDIPSData)
data(hESCs_Input)
data(DE_Input)

merged_Set = MEDIPS.mergeSets(hESCs_Input, DE_Input, name="Merged_input")

merged_Set

}
