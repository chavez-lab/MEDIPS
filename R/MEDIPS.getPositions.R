##########################################################################
##Function that identifies genomic positions of a given sequence pattern within the specified genome.
##########################################################################
##Input:	A BSgenome object, chromosomes and a sequence pattern (e.g. CG) to search for.
##Param:	BSgenome, pattern (character), chromosomes
##Output:	GRange object
##Requires:	BSgenome, Biostrings
##Modified:	11/10/2011
##Author:	Lukas Chavez, Joern Dietrich

MEDIPS.getPositions <-function(BSgenome=NULL, pattern=NULL, chromosomes=NULL){

	if(is.null(pattern)){stop("Must specify a sequence pattern!")}

	currentchr=NULL
	currentstart=NULL

	organism=ls(paste("package:", BSgenome, sep=""))
	genomedata=get(organism)

	##Determine pattern positions by accessing Biostrings
	for(chromosome in chromosomes){
		message("Searching in ", chromosome, " ...", appendLF=T)
	        {subject<-genomedata[[chromosome]]}

		plus_matches<-matchPattern(pattern,subject)
		start=start(plus_matches)
		currentchr=c(currentchr,rep(chromosome,length(start)))
		currentstart=c(currentstart,start)

		pattern<-DNAString(pattern)
		rcpattern<-reverseComplement(pattern)
		if(pattern!=rcpattern){
			minus_matches<-matchPattern(rcpattern,subject)
			start=start(minus_matches)
			currentchr=c(currentchr,rep(chromosome,length(start)))
			currentstart=c(currentstart,start)
		}
		message("[", length(start), "] found.", appendLF=T)
	}

	message("Number of identified ", pattern, " pattern: ", length(currentchr), appendLF=T)
	gc()
	return(GRanges(seqnames=currentchr, ranges=IRanges(start=currentstart, end=currentstart)))
}
