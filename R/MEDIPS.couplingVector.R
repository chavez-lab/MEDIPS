##########################################################################
##The function counts how many of the specified sequence pattern exist in each window.
##########################################################################
##Input:	Parameters that specify the sequence pattern and the size of the windows w.r.t the reference genome
##		Parameters can be specified w.r.t. a given MEDIPS SET object.
##Param:	pattern, refObj
##Output:	COUPLING SET
##Requires:	gtools, BSgenome
##Modified:	08/31/2012
##Author:	Lukas Chavez

MEDIPS.couplingVector <- function(pattern="CG", refObj=NULL){	
	if(is.list(refObj)){
		refObj=refObj[[1]]
	}
	if(class(refObj)!="MEDIPSset" & class(refObj)!="MEDIPSroiSet")	{
		stop("You must provide an MEDIPSset or MEDIPSroiSet as reference object\n")
	}
	chr_lengths=chr_lengths(refObj)
	chromosomes = chr_names(refObj)
	BSgenome = genome_name(refObj)
	if(is.list(refObj)){
		refObj=refObj[[1]]
	}
	## Get the genomic positions of the sequence pattern
	cat("Get genomic sequence pattern positions...\n")
	GRanges.pattern = MEDIPS.getPositions(BSgenome, pattern, chromosomes)
	
	## Create the genome vector Granges object	
	if(class(refObj)=="MEDIPSset"){
		window_size = window_size(refObj)
		no_chr_windows = ceiling(chr_lengths/window_size)
		supersize_chr = cumsum(no_chr_windows)	
		Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)	
	}else{	
	#class(refObj)=="MEDIPSroiSet"
		Granges.genomeVec = rois(refObj)
		window_size=-1
	}
			
	##Count the number of sequence pattern in each window.
	cat(paste("Counting the number of ",pattern, "'s in each window...\n", sep=""))
	genomeCoup = countOverlaps(Granges.genomeVec, GRanges.pattern) 

	COUPLINGsetObj = new('COUPLINGset', seq_pattern=pattern, genome_name=BSgenome, genome_CF=genomeCoup, number_pattern=length(GRanges.pattern), window_size=window_size, chr_names=chromosomes, chr_lengths=chr_lengths)
	gc()
	
	return(COUPLINGsetObj)	
}
