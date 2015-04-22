##########################################################################
##Function to calculate the number of covered genomic sequence patterns (e.g. CpGs) by given genomic regions.
##########################################################################
##Input:	bam file or tab (|) separated txt file "chr | start | stop  | strand"
##Param:	file, BSSgenome, file, pattern, extend, shift, uniq (T|F), chr.select
##Output:	seqCoverage result object 
##Requires:	gtools, BSgenome, MEDIPS.getPositions, MEDIPS.Bam2GRanges, MEDIPS.txt2Granges
##Modified:	11/10/2011
##Author:	Lukas Chavez

MEDIPS.seqCoverage<-                                  
function(file=NULL, BSgenome=NULL, pattern="CG", extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){
	
	## Proof of correctness....
	if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
	if(is.null(file)){stop("Must specify a bam or txt file.")}
		
	## Read file		
	fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
	path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/") 
	if(path==""){path=getwd()}		
	if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}
	
	dataset=get(ls(paste("package:", BSgenome, sep="")))
	if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}
	else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq, bwa=bwa)}
		
	## Sort chromosomes
	if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
	if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}
	
	## Get chromosome lengths for all chromosomes within data set.
	cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))		
	#chr_lengths=as.numeric(sapply(chromosomes, function(x){as.numeric(length(dataset[[x]]))}))
	chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])	
	
	## Get the genomic positions of the sequence pattern
	cat("Get genomic sequence pattern positions...\n")
	GRanges.pattern = MEDIPS.getPositions(BSgenome, pattern, chromosomes)	
		
	cat("Calculating sequence pattern coverage...\n")
	overlap = countOverlaps(GRanges.pattern, GRange.Reads) 
	overlap_2 = countOverlaps(GRange.Reads, GRanges.pattern)
	overlap_2 = length(overlap_2[overlap_2==0])
	
	
	seqCoverageObject=list(cov.res=overlap, pattern=pattern, numberReads=length(GRange.Reads), numberReadsWO=overlap_2)
	gc()
	return(seqCoverageObject)
}
