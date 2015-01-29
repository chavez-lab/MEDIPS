##########################################################################
##Creates MEDIPS SET of type S4 from input data (i.e. aligned short reads)
##########################################################################
##Input:	bam file or tab (|) separated txt file "chr | start | stop  | strand"
##Param:	[path]+file name, extend, shift, window_size, BSgenome, uniq (T | F), chr.select
##Output:	MEDIPS SET
##Requires:	gtools, BSgenome
##Modified:	11/1/2015
##Author:	Lukas Chavez

MEDIPS.createSet <- function(file=NULL, extend=0, shift=0, window_size=300, BSgenome=NULL, uniq=TRUE, chr.select=NULL, paired=F, sample_name=NULL, bwa=FALSE){
			
	## Proof of correctness....
	if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
	if(is.null(file)){stop("Must specify a bam or txt file.")}
		
	## Read region file		
	fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
	path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/") 
	if(path==""){path=getwd()}

	dataset=get(ls(paste("package:", BSgenome, sep="")))

	#default chromosome setting
	if(is.null(chr.select)){
		cat("All chromosomes in the reference BSgenome will be processed:\n")
		chr.select=seqnames(dataset)
		print(chr.select)
	}else{
		if(sum(!chr.select%in%seqnames(dataset))!=0){
			cat("The requested chromosome(s)", chr.select[!chr.select%in%seqnames(dataset)], "are not available in the BSgenome reference. Chromosomes available in the BSgenome reference:", seqnames(dataset), "\n")
			stop("Please check chromosome names.")
		}
	}	
	if(length(chr.select)>1){chr.select=mixedsort(unique(chr.select))}

	if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}
		ext=strsplit(fileName, ".", fixed=T)[[1]]	
	if (ext[length(ext)] %in% c("gz","zip","bz2"))
		ext=ext[-length(ext)]

	if(ext[length(ext)] %in% c("wig","bw","bigwig")){
		#read object from wig (or bigwig) file
		MEDIPSsetObj=getMObjectFromWIG(fileName, path, chr.select, BSgenome)
	}else{
		#read bam or bed file
		if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, uniq)}
		else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, uniq, bwa)}
				
		## Get chromosome lengths for all chromosomes within data set.
		chr_lengths=as.numeric(seqlengths(dataset)[chr.select])

		## Create the genome vector Granges object	
		no_chr_windows = ceiling(chr_lengths/window_size)	
		supersize_chr = cumsum(no_chr_windows)	
		Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr.select, chr_lengths, window_size)	
			
		##Distribute reads over genome
		cat("Calculating short read coverage at genome wide windows...\n")
		overlap = countOverlaps(Granges.genomeVec, GRange.Reads)
		
		datachr = unique(seqlevels(GRange.Reads))
		if(sum(!chr.select%in%datachr)!=0){
			cat("Please note, no data in the alignment file for chromosome(s):", chr.select[!chr.select%in%datachr], "\n")	
		}

		##Set sample name
		if(is.null(sample_name)){
			sample_name = fileName
		}

		## creating MEDIPSset object
    		MEDIPSsetObj = new('MEDIPSset', sample_name=sample_name,
				                        path_name=path,
				                        genome_name=BSgenome, 
							number_regions=length(GRange.Reads), 
							chr_names=chr.select, 
							chr_lengths=chr_lengths, 
							genome_count=overlap, 
							extend=extend,
							shifted=shift, 
							window_size=window_size,
							uniq=uniq)
	}
	return(MEDIPSsetObj) 	
}
