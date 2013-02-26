##########################################################################
##Merges two MEDIPS SETs created based on the same window_size and based on the same reference genome and based on the same chromosomes.
##########################################################################
##Input:	Two MEDIPS SETs
##Param:	MSet1, MSet2, name
##Output:	MEDIPS SET
##Modified:	11/15/2011
##Author:	Lukas Chavez

MEDIPS.mergeSets <- function(MSet1=NULL, MSet2=NULL, name="Merged Set"){
	
	##Proof of correctness
	#######################
	if(class(MSet1)!="MEDIPSset"){stop("You have to state a MEDIPSset object!")}
	if(class(MSet2)!="MEDIPSset"){stop("You have to state a second MEDIPSset object!")}
	
	if(genome_name(MSet1)!=genome_name(MSet2)){stop("MEDIPSset objects MSet1 and MSet2 have different reference genomes!")}
	if(length(chr_names(MSet1))!=length(chr_names(MSet2))){stop("MEDIPSset objects MSet1 and MSet2 have different chromosome sets!")}
	for(i in 1:length(chr_names(MSet1))){
		if(chr_names(MSet1)[i]!=chr_names(MSet2)[i]){stop("MEDIPSset objects MSet1 and MSet2 have different chromosome sets!")}
	}
	if(window_size(MSet1)!=window_size(MSet2)){stop("MEDIPSset objects MSet1 and MSet2 have different window sizes!")}
	
	MEDIPSsetObj = new('MEDIPSset', sample_name=name, genome_name=genome_name(MSet1), number_regions=number_regions(MSet1)+number_regions(MSet2), chr_names=chr_names(MSet1), chr_lengths=chr_lengths(MSet1), genome_count=genome_count(MSet1)+genome_count(MSet2), window_size=window_size(MSet1))
    	
	gc()
	return(MEDIPSsetObj) 	
}
