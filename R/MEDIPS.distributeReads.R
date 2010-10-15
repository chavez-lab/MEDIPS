###############
##Function takes a 
##	> reads starts vector
##	> reads stop vector
##	> genome position vector
##	> extend parameter
##extends the reads by the extend, calculates the number of reads overlapping the specified positions and
##returns a vector that contains the number of overlapping reads at each specified position.
###############

MEDIPS.distributeReads <- function(reads_start=NULL, reads_stop=NULL, reads_strand=NULL, positions=NULL, extend=0){
		
	if(extend!=0){
		reads_start[reads_strand=="-"] = reads_stop[reads_strand=="-"]-extend
		reads_stop[reads_strand=="+"] = reads_start[reads_strand=="+"]+extend
	}	
	
	reads_start_id=vector(length=length(reads_start), mode="numeric")
	reads_start_id[]=1
	
	reads_stop_id=vector(length=length(reads_stop), mode="numeric")
	reads_stop_id[]=-1
	
	positions_id=vector(length=length(positions), mode="numeric")
	positions_id[]=0	
	
	ctmatrix_pos=vector(length=length(reads_start)+length(reads_stop)+length(positions), mode="numeric")
	ctmatrix_id=vector(length=length(reads_start)+length(reads_stop)+length(positions), mode="numeric")
		
	ctmatrix_pos[]=append(append(reads_start, reads_stop), positions)
	ctmatrix_id[]=append(append(reads_start_id, reads_stop_id), positions_id)
		
	ctmatrix_id=ctmatrix_id[order(ctmatrix_pos)]
	count=cumsum(ctmatrix_id)	
		
	return(count[ctmatrix_id==0])
}
