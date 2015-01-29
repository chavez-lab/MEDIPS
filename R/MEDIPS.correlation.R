##########################################################################
##Takes a list of MEDIPS SETs, calculates correlation between all sets and
##reports a correlation matrix.
##########################################################################
##Input:	List of MEDIPS SETs
##Output:	correlation matrix
##Modified:	04/20/2012
##Author:	Lukas Chavez

MEDIPS.correlation <- function(MSets=NULL, plot=T, method="pearson"){
	
	n=length(MSets)

	cor.matrix=matrix(ncol=n, nrow=n)
	c.names=NULL

	for(i in 1:n){
		c.names=c(c.names, sample_name(MSets[[i]]))
	
		for(j in i:n){
			#print(paste("Comparing sample ", i, " to sample ", j, "...", sep=""))
			##Proof of correctness
			#######################
			if(class(MSets[[i]])!="MEDIPSset"){stop("You have to state a MEDIPSset object!")}
			if(class(MSets[[j]])!="MEDIPSset"){stop("You have to state a second MEDIPSset object!")}
			if(genome_name(MSets[[i]])!=genome_name(MSets[[j]])){stop("MEDIPSset objects MSet1 and MSet2 have different reference genomes!")}
			if(length(chr_names(MSets[[i]]))!=length(chr_names(MSets[[j]]))){stop("MEDIPSset objects MSet1 and MSet2 have different chromosome sets!")}
			for(k in 1:length(chr_names(MSets[[i]]))){
				if(chr_names(MSets[[i]])[k]!=chr_names(MSets[[j]])[k]){stop("MEDIPSset objects MSet1 and MSet2 have different chromosome sets!")}
			}
			if(window_size(MSets[[i]])!=window_size(MSets[[j]])){stop("MEDIPSset objects MSet1 and MSet2 have different window sizes!")}		
		
			c=cor(MSets[[i]]@genome_count, MSets[[j]]@genome_count, method=method)
			cor.matrix[i,j]=c
			if(plot & i!=j){
				png(paste("Scatter_", sample_name(MSets[[i]]), "_vs_", sample_name(MSets[[j]]), ".png", sep=""))
				plot(log2(MSets[[i]]@genome_count), log2(MSets[[j]]@genome_count), pch=".", main=paste("Scater_", sample_name(MSets[[i]]), "_vs_", sample_name(MSets[[j]]),  sep=""), sub=paste(method, " correlation: ", round(cor.matrix[i,j], digits=2), sep=""), xlab=paste(sample_name(MSets[[i]]), " log2(counts)", sep=""), ylab=paste(sample_name(MSets[[j]]), " log2(counts)", sep=""))
				dev.off()
			}
		}
	}

	colnames(cor.matrix)=c.names
	rownames(cor.matrix)=c.names

	gc()
	return(cor.matrix)		
}
