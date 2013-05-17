##########################################################################
##Function merges neighbouring genomic coordinates
##########################################################################
##Input:	table containing row-wise genomic coordinates
##Param:	frames, distance
##Output:	merged table of genomic coordinates
##Modified:	11/23/2011
##Author:	Lukas Chavez

MEDIPS.mergeFrames = function(frames=NULL, distance=1){
	
	if(is.null(frames)){stop("Must specify a frame set.")}	
	
	chromosomes=as.character(unique(frames[,1]))
	
	output=NULL
	
	for(i in 1:length(chromosomes)){
		start=as.numeric(frames[,2][as.character(frames[,1])==chromosomes[i]])
		stop=as.numeric(frames[,3][as.character(frames[,1])==chromosomes[i]])
		stop_ext=stop+(distance)
			
		start_ID=vector(length=length(start), mode="numeric") 
		start_ID[]=1					
			
		stop_ID=vector(length=length(stop), mode="numeric") 
		stop_ID[]=-1
			
		positions=vector(length=(length(start)+length(stop)), mode="numeric") 
		positions[]=append(start, stop_ext)
		
		positions_orig=vector(length=(length(start)+length(stop)), mode="numeric")
		positions_orig[]=append(start, stop)
		
		IDs=vector(length=(length(start)+length(stop)), mode="numeric") 
		IDs[]=append(start_ID, stop_ID)
		
		IDs=IDs[order(positions)]
		positions_orig=positions_orig[order(positions)]
		positions=positions[order(positions)]
		
		merged_IDs=cumsum(IDs)		
				
		new_stop=positions_orig[merged_IDs==0]
		indices_stop=which(merged_IDs==0)

		new_start = positions_orig[1]
		new_start = c(new_start, positions_orig[c(indices_stop+1)])
		new_start = new_start[-length(new_start)]

		output=rbind(output, cbind(chromosomes[i], new_start, new_stop))
	}

	output = cbind(output, paste("ID", seq(1, nrow(output)), sep="_"))

	print(paste("Number of merged frames: ", length(output[,1]), sep=""), quote=F)
	output=list(chr=as.character(output[,1]), start=format(as.numeric(output[,2]), scientific=F, trim=T), stop=format(as.numeric(output[,3]), scientific=F, trim=T), ID=as.character(output[,4]))
	return(as.data.frame(output,stringsAsFactors=F))	
}

