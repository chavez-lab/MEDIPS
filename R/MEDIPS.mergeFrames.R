MEDIPS.mergeFrames = function(frames=NULL){
	
	if(is.null(frames)){stop("Must specify a frame set.")}	
	
	chromosomes=as.character(unique(frames[,1]))
	
	output=NULL
	
	for(i in 1:length(chromosomes)){
		start=as.numeric(frames[,2][as.character(frames[,1])==chromosomes[i]])
		stop=as.numeric(frames[,3][as.character(frames[,1])==chromosomes[i]])
			
		start_ID=vector(length=length(start), mode="numeric") 
		start_ID[]=1					
			
		stop_ID=vector(length=length(stop), mode="numeric") 
		stop_ID[]=-1
			
		positions=vector(length=(length(start)+length(stop)), mode="numeric") 
		positions[]=append(start, stop)
		
		IDs=vector(length=(length(start)+length(stop)), mode="numeric") 
		IDs[]=append(start_ID, stop_ID)
		
		IDs=IDs[order(positions)]
		positions=positions[order(positions)]		
		
		merged_IDs=cumsum(IDs)		
				
		new_stop=positions[merged_IDs==0]
		indices_stop=which(merged_IDs==0)

		if(length(indices_stop)!=1){
				
			indices_stop=indices_stop+1
				
			indices_stop=indices_stop[1:(length(indices_stop)-1)]
			new_start=positions[indices_stop]
			new_start=c(positions[1], new_start)
		
			tmp_out=cbind(chromosomes[i], new_start, new_stop)
			output=rbind(output, tmp_out)
		}
		else{
			output=rbind(output, cbind(chromosomes[i], positions[1], positions[length(positions)]))
		}
		
	}
	
	print(paste("Number of merged frames: ", length(output[,1]), sep=""), quote=F)
	output=list(chr=as.character(output[,1]), new_start=format(as.numeric(output[,2]), scientific=F), new_stop=format(as.numeric(output[,3]), scientific=F))
	return(as.data.frame(output,stringsAsFactors=F))	
}

