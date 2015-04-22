##########################################################################
##Estimates the reproducibility of genome wide short read coverages 
##by sampling subsets of the given data and comparing accordance.
##########################################################################
##Input:	bam file or tab (|) separated txt file "chr | start | stop  | strand"
##Param:	file, BSSgenome, nit, nrit, empty_bins, rank, extend, shift, window_size, uniq (T|F), chr.select
##Output:	Saturation result object
##Requires:	gtools, BSgenome, MEDIPS.Bam2GRanges, MEDIPS.txt2Granges, MEDIPS.GenomicCoordinates
##Modified:	11/10/2012
##Author:	Lukas Chavez


MEDIPS.saturation<-
function(file=NULL, BSgenome=NULL, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE, extend=0, shift=0,  window_size=500, uniq=1e-3, chr.select=NULL, paired=F){
	
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
	
	## Create the genome vector Granges object	
	no_chr_windows = ceiling(chr_lengths/window_size)
	supersize_chr = cumsum(no_chr_windows)	
	Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)	
				
	####################################
	#The first run (a1) is to calculate the saturation based on two distinct subsets of given regions.
	#The second run (a2) is to estimate a saturation based on an artifically doubled set of given regions.
	####################################
	for(a in 1:2){		
		if(a==1){
			distinctSets_size=floor(length(GRange.Reads)/2)
			subset_size=floor(distinctSets_size/nit)
			cat("Saturation analysis...\n")
		}		
		if(a==2){
			cat("Estimated saturation analysis...\n")
			distinctSets_size=length(GRange.Reads)
			GRange.Reads = c(GRange.Reads, GRange.Reads)			
			nit=nit*2
		}			
						
		####################################
		#The loops: First one is  for repeating the random process and results are averaged at the end.
		####################################
		for(r in 1:nrit){			
			cat(paste("Random iteration: ",r,"/", nrit, "...\n", sep=""))
						
			correlation=0
			no_considered_reads=0
	
			##Create a random order for all given regions.
			random=sample(1:length(GRange.Reads), length(GRange.Reads))		
			
			####################################
			#The loops: Second one is  for calculating correlations for subsets of the current total set of given regions.
			####################################
			for(l in 1:nit){						
				cat(paste("Processing subset ", l, "/", nit, "...\n", sep=""))
				
				subset_start_one=(l-1)*subset_size+1			
				subset_start_two=((l-1)*subset_size+1)+distinctSets_size
				
				if(l<nit){	
					subset_stop_one=l*subset_size
					subset_stop_two=(l*subset_size)+distinctSets_size												
					no_considered_reads=c(no_considered_reads, l*subset_size)
				}
				if(l==nit){	
					subset_stop_one=distinctSets_size							
					subset_stop_two=length(GRange.Reads)					
					no_considered_reads=c(no_considered_reads, distinctSets_size)
				}		
				
				##Calculate coverage							
				if(l==1){
					genome_value_a = countOverlaps(Granges.genomeVec, GRange.Reads[c(random[subset_start_one:subset_stop_one])]) 
					genome_value_b = countOverlaps(Granges.genomeVec, GRange.Reads[c(random[subset_start_two:subset_stop_two])]) 
				}
				else{
					genome_value_a = genome_value_a + countOverlaps(Granges.genomeVec, GRange.Reads[c(random[subset_start_one:subset_stop_one])]) 
					genome_value_b = genome_value_b + countOverlaps(Granges.genomeVec, GRange.Reads[c(random[subset_start_two:subset_stop_two])]) 
				}				
								
				##Calculate correlation for each interval							
				if(empty_bins){
					if(!rank){correlation=c(correlation, cor(genome_value_a, genome_value_b))}
					if(rank){correlation=c(correlation, cor(rank(genome_value_a), rank(genome_value_b)))}
					
				}
				if(!empty_bins){
					if(!rank){correlation=c(correlation, cor(genome_value_a[genome_value_a!=0 & genome_value_b!=0], genome_value_b[genome_value_a!=0 & genome_value_b!=0]))}
					if(rank){correlation=c(correlation, cor(rank(genome_value_a[genome_value_a!=0 & genome_value_b!=0]), rank(genome_value_b[genome_value_a!=0 & genome_value_b!=0])))}
				}
			}			
	
			##Store the results of each iteration step.
			if(r==1){
				output_matrix=matrix(ncol=2, nrow=length(correlation))
				output_matrix[,1]=no_considered_reads
				output_matrix[,2]=correlation
			}
			if(r>1){
				output_matrix[,2]=output_matrix[,2]+correlation
			}					
		
		}

		##Correct for the number of Monte-Carlo iterations
		output_matrix[,2]=output_matrix[,2]/nrit
		
		if(a==1){
			distinct=output_matrix
			maxTruCor=c(output_matrix[length(correlation),1], output_matrix[length(correlation),2])
		}
		if(a==2){
			maxEstCor=c(output_matrix[length(correlation),1], output_matrix[length(correlation),2])
			saturationObj=list(distinctSets=distinct, estimation=output_matrix, numberReads=(length(GRange.Reads)/2), maxEstCor=maxEstCor, maxTruCor=maxTruCor)
		}
	}
	gc()				
	return(saturationObj)
}
