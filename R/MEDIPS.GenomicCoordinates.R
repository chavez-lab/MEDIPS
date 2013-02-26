##########
#Function calculates genome wide coordinates based on reference genome and w.r.t specified window size.
#A GRange object for the genome wide windows is returned.
##########
##Input:	Parameters that specify the reference genome and targeted window size.
##Param:	supersize_chr, no_chr_windows_, chromosomes, chr_lengths, window_size
##Output:	GRange object
##Modified:	11/10/2011
##Author:	Lukas Chavez

MEDIPS.GenomicCoordinates <- function(supersize_chr=NULL, no_chr_windows=NULL, chromosomes=NULL, chr_lengths=NULL, window_size=NULL){
	
	cat("Calculating genomic coordinates...")
	
	genomeVec_chr=vector(length=supersize_chr[length(chromosomes)], mode="character")
	genomeVec_pos=vector(length=supersize_chr[length(chromosomes)], mode="numeric")		
	
	total=length(chromosomes)
 	for(i in 1:length(chromosomes)){		
		if(i==1){			
			genomeVec_chr[1:no_chr_windows[i]]=chromosomes[i]
			genomeVec_pos[1:no_chr_windows[i]]=seq(1, chr_lengths[i], window_size)
		}
		if(i>1){
			genomeVec_chr[(supersize_chr[i-1]+1):(supersize_chr[i-1]+no_chr_windows[i])]=chromosomes[i]
			genomeVec_pos[(supersize_chr[i-1]+1):(supersize_chr[i-1]+no_chr_windows[i])]=seq(1, chr_lengths[i], window_size)
		}
        }
	
	cat("\nCreating Granges object for genome wide windows...\n")
	ROIs = GRanges(seqnames=genomeVec_chr, ranges=IRanges(start=genomeVec_pos, end=genomeVec_pos+window_size-1))
	
	gc()
	return(ROIs)
	
}
