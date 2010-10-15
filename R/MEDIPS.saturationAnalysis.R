####################################
#Function to simulate the re-producibility of MeDIP-Seq experiments by sampling subsets of the regions data and comparing data accordance.
####################################

MEDIPS.saturationAnalysis<-
function(data=NULL, no_iterations=10, no_random_iterations=1, empty_bins=TRUE, rank=FALSE, extend=400, bin_size=NULL){
	
	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")
	regions_chr=regions_chr(data)
	regions_start=regions_start(data)
	regions_stop=regions_stop(data)
	regions_strand=regions_strand(data)
	chr_names=chr_names(data)
	chr_lengths=chr_lengths(data)
	number_regions=number_regions(data)
	if(is.null(extend)){extend=extend(data)}
	if(is.null(bin_size)){bin_size=bin_size(data)}	
	chr_lengths=c(0, chr_lengths)
	supersize=cumsum(chr_lengths)
	genome_length=supersize[length(supersize)]		
		
	#####
	#Preprocessing, regions.
	#####
	for(i in 1:(length(chr_names))){		
		cat(paste("Preprocessing ", chr_names[i], "...\n", sep=""))		
		if(extend!=0){
			regions_start[regions_strand=="-"] = regions_stop[regions_strand=="-"]-extend
			regions_stop[regions_strand=="+"] = regions_start[regions_strand=="+"]+extend
		}		
		regions_start[regions_chr==chr_names[i]] = (regions_start[regions_chr==chr_names[i]]+supersize[i])
		regions_stop[regions_chr==chr_names[i]] = (regions_stop[regions_chr==chr_names[i]]+supersize[i])
	}
	
	##Calculate the number of genomic windows.
	no_windows=ceiling(genome_length/bin_size)	
	
	##Initialization of the genome vectors.
	cat("Initialization of the temporary genome vector....\n")
	genome_pos=as.numeric(seq(1, (no_windows*bin_size), bin_size))
	genome_value_a=vector(length=no_windows, mode="numeric")
	genome_value_b=vector(length=no_windows, mode="numeric")
	
	##Define the genome vector update functions.
	update_vec_one=function(x){genome_value_a[results_start[random[x]]:results_stop[random[x]]]<<-genome_value_a[results_start[random[x]]:results_stop[random[x]]]+1}
	update_vec_two=function(x){genome_value_b[results_start[random[x]]:results_stop[random[x]]]<<-genome_value_b[results_start[random[x]]:results_stop[random[x]]]+1}
		
	####################################
	#The first process is to calculate true saturation based on two distinct subsets of given regions.
	#The second process is to estimate a saturation based on an artifically doubled set of given regions.
	####################################
	for(a in 1:2){		
		if(a==1){
			distinctSets_size=floor(number_regions/2)
			subset_size=floor(distinctSets_size/no_iterations)
			cat("Saturation analysis...\n")
		}		
		if(a==2){
			regions_chr=rbind(regions_chr, regions_chr)
			regions_start=rbind(regions_start, regions_start)
			regions_stop=rbind(regions_stop, regions_stop)
			distinctSets_size=number_regions
			no_iterations=no_iterations*2
			cat("Estimated saturation analysis...\n")
		}			
						
		####################################
		#The loops: First one is  for repeating the random process and results are averaged at the end.
		####################################
		for(r in 1:no_random_iterations){			
			cat(paste("Random iteration: ",r,"/", no_random_iterations, "...\n", sep=""))
					
			##Clear the genome vector.		
			genome_value_a[]=0
			genome_value_b[]=0
			correlation=0
			no_considered_reads=0
	
			##Distribute regions over genome.
			#print("Find intervals...")
			results_start=findInterval(regions_start, genome_pos)
			results_stop=findInterval(regions_stop, genome_pos)
				
			##Create a random order for all given regions.
			random=sample(1:length(regions_start), length(regions_start))		
		
			####################################
			#The loops: Second one is  for calculating correcaltions for subsets of the current total set of given regions.
			####################################
			for(l in 1:no_iterations){						
				cat(paste("Processing subset ", l, "/", no_iterations, "...\n", sep=""))
								
				subset_start_one=(l-1)*subset_size+1			
				subset_start_two=((l-1)*subset_size+1)+distinctSets_size
				
				if(l<no_iterations){	
					subset_stop_one=l*subset_size			
					subset_stop_two=(l*subset_size)+distinctSets_size												
					no_considered_reads=c(no_considered_reads, l*subset_size)
				}
				if(l==no_iterations){	
					subset_stop_one=distinctSets_size							
					subset_stop_two=length(regions_chr)					
					no_considered_reads=c(no_considered_reads, distinctSets_size)
				}		
			
				apply(data.matrix(subset_start_one:subset_stop_one), 1, update_vec_one)
				apply(data.matrix(subset_start_two:subset_stop_two), 1, update_vec_two)

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
	
			##Store the results of the iteration step.
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
		output_matrix[,2]=output_matrix[,2]/no_random_iterations
		
		if(a==1){
			distinct=output_matrix
			maxTruCor=c(output_matrix[length(correlation),1], output_matrix[length(correlation),2])
		}
		if(a==2){
			maxEstCor=c(output_matrix[length(correlation),1], output_matrix[length(correlation),2])
			saturationObj=list(distinctSets=distinct, estimation=output_matrix, numberReads=number_regions(data), maxEstCor=maxEstCor, maxTruCor=maxTruCor)
		}
	}
	gc()				
	return(saturationObj)
}
