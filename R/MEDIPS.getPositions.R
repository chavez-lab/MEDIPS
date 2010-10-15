#####
#Function that identifies genomic positions of a given sequence pattern within the specified genome.
#####

MEDIPS.getPositions <-function(data=NULL, pattern=NULL){
	
	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")
	if(is.null(pattern)){stop("Must specify a sequence pattern like CG.")}	
	
	currentchr=NULL
	currentstart=NULL
	
	organism=ls(paste("package:", genome_name(data),sep=""))
	genomedata=get(organism)
	chromosomes=chr_names(data)	
	
	##Function which creates position files by accessing Biostrings
	#Function which creates the matrix  
	total=length(chromosomes)
        pb <- txtProgressBar(min = 0, max = total, style = 3)
	for(chromosome in chromosomes){
	        {subject<-genomedata[[chromosome]]}
		setTxtProgressBar(pb, which(chromosomes==chromosome))
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
	}		
 	MEDIPSsetObj = new('MEDIPSset', seq_pattern=as.character(pattern), pattern_chr=currentchr, pattern_pos=currentstart, number_pattern=length(currentstart), genome_chr=genome_chr(data), genome_pos=genome_pos(data), genome_raw=genome_raw(data), extend=extend(data), bin_size=bin_size(data), sample_name=sample_name(data), genome_name=genome_name(data), regions_chr=regions_chr(data), regions_start=regions_start(data), regions_stop=regions_stop(data), regions_strand=regions_strand(data), number_regions=number_regions(data), chr_names=chr_names(data), chr_lengths=chr_lengths(data))
	cat("\n")

	return(MEDIPSsetObj)	
}
