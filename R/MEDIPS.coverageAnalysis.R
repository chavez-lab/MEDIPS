####################################
#Function to calculate the number of covered genomic sequence patterns (e.g. CpGs) by given genomic regions.
####################################

MEDIPS.coverageAnalysis<-                                  
function(data=NULL, coverages=c(1,2,3,4,5,10), no_iterations=10, no_random_iterations=1, extend=NULL){
	
	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")
	regions_chr=regions_chr(data)
	regions_start=regions_start(data)
	regions_stop=regions_stop(data)
	regions_strand=regions_strand(data)
	pattern_chr=pattern_chr(data)
	pattern_pos=pattern_pos(data)	
	chr_names=chr_names(data)
	chr_lengths=chr_lengths(data)
	number_regions=number_regions(data)
	if(is.null(extend)){extend=extend(data)}
	pattern=seq_pattern(data)
	number_pattern=number_pattern(data)
	
	###################################
	#Calculate the subset size based on the number of iteration steps and on the total number of given regions. 
	###################################
	subset_size=floor(number_regions/no_iterations)
			
	output_matrix=matrix(ncol=(length(coverages)+1), nrow=no_iterations+2)
	output_matrix[,]=0
	output_matrix[1,]=c(0, coverages)
				
	#####
	#Preprocessing regions and pattern.
	#####
	chr_lengths=c(0, chr_lengths)
	supersize=cumsum(chr_lengths)		
	for(i in 1:(length(chr_names))){		
		cat(paste("Preprocessing ", chr_names[i], "...\n", sep=""))			
		if(extend!=0){
			regions_start[regions_strand=="-"] = regions_stop[regions_strand=="-"]-extend
			regions_stop[regions_strand=="+"] = regions_start[regions_strand=="+"]+extend
		}		
		regions_start[regions_chr==chr_names[i]] = (regions_start[regions_chr==chr_names[i]]+supersize[i])
		regions_stop[regions_chr==chr_names[i]] = (regions_stop[regions_chr==chr_names[i]]+supersize[i])
		pattern_pos[pattern_chr==chr_names[i]] = pattern_pos[pattern_chr==chr_names[i]]+supersize[i]		
	}	

	####################################
	#The loop 
	#1. extracts in each iteration step an increasing (according to the subset size) set of random selected reads.
	#2. performes a preprocessing for the subsequent
	#3. findInterval and
	#4. stores at each iteration step the number and depth of covered pattern (e.g. CpGs).
	####################################
	for(r in 1:no_random_iterations){	
		cat(paste("Random iteration: ",r,"/", no_random_iterations, "...\n", sep=""))		
		#random field assigns each region to a random number from 1 to number of regions
		random=sample(1:number_regions, number_regions) 
		
		for(i in 1:no_iterations){
			cat(paste("Processing sub-set: ",i,"/", no_iterations, "...\n", sep=""))		
	
			#####Extract current considered regions
			if(i<no_iterations){	
				regions_subset_start=regions_start[random<=(i*subset_size)]
				regions_subset_stop=regions_stop[random<=(i*subset_size)]
			}
			if(i==no_iterations){
				regions_subset_start=regions_start
				regions_subset_stop=regions_stop
			}
			
			current_subset_size=length(regions_subset_start)
			
			##Set the subset size within each iteration step (only necessary one time for each iteration step). 
			if(r==1){output_matrix[i+2,1]=length(regions_subset_start)}
			
			#####Transform selected regions into supersized vectors.
			regions_transformed_pos=vector(length=(2*current_subset_size), mode="numeric")
			regions_transformed_id=vector(length=(2*current_subset_size), mode="numeric")
			regions_transformed_pos[1:current_subset_size]=regions_subset_start
			regions_transformed_id[1:current_subset_size]=1			
			regions_transformed_pos[(current_subset_size+1):(2*current_subset_size)]=regions_subset_stop
			regions_transformed_id[(current_subset_size+1):(2*current_subset_size)]=-1

			#####Sort the one-vector representations of the regions 
			regions_transformed_pos_ordered=regions_transformed_pos[order(regions_transformed_pos)]
			regions_transformed_id_ordered=regions_transformed_id[order(regions_transformed_pos)]

			#####delete temporary vectors
			rm(regions_subset_start, regions_subset_stop, regions_transformed_pos, regions_transformed_id)
	
			#####Calculate the open_intervals vector 
			open_intervals=cumsum(regions_transformed_id_ordered)	
	
			#####Apply findInterval()
			results=findInterval(pattern_pos, regions_transformed_pos_ordered)
	
			#####For each questioned depth of coverage
			for(j in 1:length(coverages)){
				if(r==1){output_matrix[i+2,j+1]=(length(open_intervals[results][open_intervals[results]>=coverages[j]]))}
				if(r>1){output_matrix[i+2,j+1]=(output_matrix[i+2,j+1]+(length(open_intervals[results][open_intervals[results]>=coverages[j]])))}				
			}			
		}
	}

	#####Correct for the number of Monte-Carlo iterations
	for(i in 1:no_iterations){
		for(j in 1:length(coverages)){
			output_matrix[i+2,j+1]=output_matrix[i+2,j+1]/no_random_iterations
		}
	}
	
	coveredPos=matrix(ncol=length(coverages), nrow=3)
	coveredPos[1,]=coverages
	coveredPos[2,]=output_matrix[no_iterations+2, 2:(length(coverages)+1)]
	coveredPos[3,]=round((coveredPos[2,]/number_pattern), digits=2)
	
	coverageObject=list(matrix=output_matrix, maxPos=number_pattern, pattern=pattern, coveredPos=coveredPos)
	gc()
	return(coverageObject)
}
