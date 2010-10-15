#####
#Function that reads a region file.
#Region file is a simple tab-delimited text files (chr | start | stop | strand).
#####

MEDIPS.readAlignedSequences <-
function(file=NULL, BSgenome=NULL, numrows=-1){
		
	## Species dependent definitions.
	if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
	if(is.null(file)){stop("Must specify a region file.")}
		
	## Read region file	
	regions=NULL
	
	fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
	path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1],collapse="/") 
	if(path==""){
		path=getwd()
	}
	
	cat(paste("Reading file ", fileName, " in ", path, "...\n", sep=""))		
	
	if(!fileName%in%dir(path)){
		stop(paste("File", fileName, " not found in", path, sep =" "))				
	}
	
	regions=read.table(file, sep='\t', header=FALSE, row.names=NULL, nrows=numrows, colClasses=c("character", "numeric", "numeric", "character"))			
		
	## Sort chromosomes
	if(length(unique(regions[,1]))>1){chromosomes=mixedsort(unique(regions[,1]))}
	if(length(unique(regions[,1]))==1){chromosomes=unique(regions[,1])}
	
	## Get chromosome lengths for all chromosomes within data set.
	cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))		
	dataset=get(ls(paste("package:", BSgenome, sep="")))
	chr_lengths=as.numeric(sapply(chromosomes, function(x){as.numeric(length(dataset[[x]]))}))
		
	## creating MEDIPSset object
    	MEDIPSsetObj = new('MEDIPSset', sample_name=fileName, genome_name=BSgenome, regions_chr=regions[,1], regions_start=regions[,2], regions_stop=regions[,3], regions_strand=regions[,4], number_regions=length(regions[,1]), chr_names=chromosomes, chr_lengths=chr_lengths)
    	return(MEDIPSsetObj)    
}
