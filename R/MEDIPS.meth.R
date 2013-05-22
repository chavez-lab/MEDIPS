##########################################################################
##Function calculates genome wide methylation signals for each provided set and
##calculates differential coverage between groups of MEDIPS SETs (if provided) and
##calculates CNV scores between groups of INPUT SETs (if provided)
##########################################################################
##Input:	Groups of MEDIPS and INPUT SETs
##Param:	MSet1, MSet2, CSet, ISet1, ISet2, chr, p.adj, diff.method, prob.method, type
##Output:	Result table
##Requires:	DNAcopy
##Modified:	11/22/2011
##Author:	Lukas Chavez, Joern Dietrich

MEDIPS.meth = function(
		MSet1 = NULL, 
		MSet2 = NULL, 
		CSet = NULL, 
		ISet1 = NULL, 
		ISet2 = NULL, 
		chr = NULL, 
		p.adj="bonferroni", 
		diff.method="ttest", 
		prob.method="poisson",
		CNV=FALSE,
		MeDIP=TRUE,
		type="rpkm",
		minRowSum=1
		)
{
	nMSets1 = length(MSet1)	
	nMSets2 = length(MSet2)	
	nISets1 = length(ISet1)	
	nISets2 = length(ISet2)	
	if(is.list(CSet)) CSet=CSet[[1]]
	if(!is.list(MSet1)) MSet1=c(MSet1)
	if(!is.list(MSet2)) MSet2=c(MSet2)
	if(!is.list(ISet1)) ISet1=c(ISet1)
	if(!is.list(ISet2)) ISet2=c(ISet2)
	
	##Proof of correctness
	#######################
	if(MeDIP){
		if(class(CSet)!="COUPLINGset"){stop("You have to state a COUPLINGset object!")}
	}

	controlSet = MSet1[[1]]

	for(i in 1:nMSets1){
		if(class(MSet1[[i]])!="MEDIPSset" & class(MSet1[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
		if(length(genome_count(MSet1[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
	}
	if(!is.null(MSet2)){
    		for(i in 1:nMSets2){
			if(class(MSet2[[i]])!="MEDIPSset" & class(MSet2[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
			if(length(genome_count(MSet2[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
    		}
	}	
	if(!is.null(ISet1)){
		for(i in 1:nISets1){
			if(class(ISet1[[i]])!="MEDIPSset" & class(ISet1[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
			if(length(genome_count(ISet1[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
		}
	}
	if(!is.null(ISet2)){
		for(i in 1:nISets2){
			if(class(ISet2[[i]])!="MEDIPSset" & class(ISet2[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
			if(length(genome_count(ISet2[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
		}
	}


	##Data preparation
	###################
	##Calculate genomic coordinates
	##################################
	if(class(controlSet)=="MEDIPSset"){
		window_size = window_size(controlSet)
		no_chr_windows = ceiling(chr_lengths(controlSet)/window_size(controlSet))
		supersize_chr = cumsum(no_chr_windows)	
		GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr_names(controlSet), chr_lengths(controlSet), window_size(controlSet))
	}
	else if(class(controlSet)=="MEDIPSroiSet"){
		GRanges.genome = rois(controlSet)
		window_size=as.numeric(width(GRanges.genome))
	} 
	
	##Create data frame for all genomic windows and MEDIPS SETs
	if(MeDIP){
		base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), stop=end(GRanges.genome), CF=genome_CF(CSet), stringsAsFactors=F)
	}
	else{
		base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), stop=end(GRanges.genome), stringsAsFactors=F)
	}

	rm(controlSet)
	gc()

	counts.medip = NULL
	rpkm.medip = NULL
	rms = NULL	
	prob = NULL	
	counts.input = NULL
	rpkm.input = NULL
	
	##Add counts
	if(!is.null(MSet1)){
		for(i in 1:nMSets1){
			cat(paste("Preprocessing MEDIPS SET ", i, " in MSet1...\n", sep=""))
			counts.medip = cbind(counts.medip, MSet1=genome_count(MSet1[[i]]))
			rpkm.medip = cbind(rpkm.medip, ((genome_count(MSet1[[i]])*10^9)/(window_size*number_regions(MSet1[[i]]))))
			if(MeDIP){
				ccObj = MEDIPS.calibrationCurve(MSet=MSet1[[i]], CSet=CSet, input=F)
				rms = cbind(rms, MEDIPS.rms(MSet1[[i]], CSet, ccObj=ccObj))
				if(prob.method=="poisson"){prob = cbind(prob, MEDIPS.pois(MSet1[[i]], CSet, ccObj=ccObj))}
				else if(prob.method=="negBinomial"){prob = cbind(prob, MEDIPS.negBin(MSet1[[i]], CSet, ccObj=ccObj))}		
				else{stop(paste("Method ", prob.method, " not supported.",  sep=""))}
			}		
		}
	}	
	if(!is.null(MSet2)){
		 for(i in 1:nMSets2){
			cat(paste("Preprocessing MEDIPS SET ", i, " in MSet2...\n", sep=""))
			counts.medip = cbind(counts.medip, MSet2=genome_count(MSet2[[i]]))
			rpkm.medip = cbind(rpkm.medip, ((genome_count(MSet2[[i]])*10^9)/(window_size*number_regions(MSet2[[i]]))))
			if(MeDIP){
				ccObj = MEDIPS.calibrationCurve(MSet=MSet2[[i]], CSet=CSet, input=F)
				rms = cbind(rms, MEDIPS.rms(MSet2[[i]], CSet, ccObj=ccObj))
				if(prob.method=="poisson"){prob = cbind(prob, MEDIPS.pois(MSet2[[i]], CSet, ccObj=ccObj))}
				else if(prob.method=="negBinomial"){prob = cbind(prob, MEDIPS.negBin(MSet2[[i]], CSet, ccObj=ccObj))}
				else{stop(paste("Method ", prob.method, " not supported.",  sep=""))}			
			}
	  	 }
	}	
	if(!is.null(ISet1)){		
		for(i in 1:nISets1){
			cat(paste("Preprocessing INPUT SET ", i, " in ISet1...\n", sep=""))
			counts.input = cbind(counts.input, ISet1=genome_count(ISet1[[i]]))
			rpkm.input = cbind(rpkm.input, ((genome_count(ISet1[[i]])*10^9)/(window_size*number_regions(ISet1[[i]]))))		
		   }
	}	
	if(!is.null(ISet2)){
		for(i in 1:nISets2){
			cat(paste("Preprocessing INPUT SET ", i, " in ISet2...\n", sep=""))
			counts.input = cbind(counts.input, ISet2=genome_count(ISet2[[i]]))
			rpkm.input = cbind(rpkm.input, ((genome_count(ISet2[[i]])*10^9)/(window_size*number_regions(ISet2[[i]]))))
		 	}
	}		
	
	##Extract data for selected chromosome
	#######################################
	if(!is.null(chr)){
		fi=base[,1]%in%chr
		
		cat("Extracting data for", chr, "...\n", sep=" ")
		if(length(fi)==0){stop("Stated chromosome does not exist in the COUPLING SET.")}
		if(!is.null(counts.medip)){
			counts.medip = counts.medip[fi,]
			rpkm.medip = rpkm.medip[fi,]
			rms = rms[fi,]
			prob = prob[fi,]
		}
		if(!is.null(counts.input)){
			counts.input = counts.input[fi,]
			rpkm.input = rpkm.input[fi,]
		}
		base = base[fi,]
		cat(nrow(base), "windows on", chr, "\n",sep=" ")		
	}
	
	##Set colnames and transform to data.frames
	############################################
	col.names.count = NULL
	col.names.rpkm = NULL
	col.names.rms = NULL
	col.names.prob = NULL	
	if(nMSets1!=0){
		for(i in 1:nMSets1){
			col.names.count = c(col.names.count, paste(sample_name(MSet1[[i]]), ".counts", sep=""))
			col.names.rpkm = c(col.names.rpkm, paste(sample_name(MSet1[[i]]), ".rpkm", sep=""))
			if(MeDIP){
				col.names.rms = c(col.names.rms, paste(sample_name(MSet1[[i]]), ".rms", sep=""))
				col.names.prob = c(col.names.prob, paste(sample_name(MSet1[[i]]), ".prob", sep=""))
			}
		}		
	}
	if(nMSets2!=0){
		for(i in 1:nMSets2){
			col.names.count = c(col.names.count, paste(sample_name(MSet2[[i]]), ".counts", sep=""))
			col.names.rpkm = c(col.names.rpkm, paste(sample_name(MSet2[[i]]), ".rpkm", sep=""))
			if(MeDIP){
				col.names.rms = c(col.names.rms, paste(sample_name(MSet2[[i]]), ".rms", sep=""))
				col.names.prob = c(col.names.prob, paste(sample_name(MSet2[[i]]), ".prob", sep=""))
			}
		}		
	}
	if(nMSets1!=0 | nMSets2!=0){
		
		counts.medip = data.frame(counts.medip)
		colnames(counts.medip) = col.names.count
		rpkm.medip =  data.frame(rpkm.medip)
		colnames(rpkm.medip) = 	col.names.rpkm
		if(MeDIP){
			rms = data.frame(rms)
			colnames(rms) = col.names.rms
			prob =  data.frame(prob)
			colnames(prob) = col.names.prob
		}
	}		
	
	col.names.count.input = NULL
	col.names.rpkm.input = NULL
	if(nISets1!=0){
		for(i in 1:nISets1){
			col.names.count.input = c(col.names.count.input, paste(sample_name(ISet1[[i]]), ".counts", sep=""))
			col.names.rpkm.input = c(col.names.rpkm.input, paste(sample_name(ISet1[[i]]), ".rpkm", sep=""))
		}		
	}
	if(nISets2!=0){
		for(i in 1:nISets2){
			col.names.count.input = c(col.names.count.input, paste(sample_name(ISet2[[i]]), ".counts", sep=""))
			col.names.rpkm.input = c(col.names.rpkm.input, paste(sample_name(ISet2[[i]]), ".rpkm", sep=""))
		}		
	}
	
	if(nISets1!=0 | nISets2!=0){		
		
		counts.input = data.frame(counts.input)
		colnames(counts.input) = 	col.names.count.input
		rpkm.input = data.frame(rpkm.input)	
		colnames(rpkm.input) = 		col.names.rpkm.input
	}
		
	
	##If two groups of MEDIPS SETs are given 
	##calculate differential coverage
	##################################
	if(!is.null(MSet1) & !is.null(MSet2)){
		cat(paste("Differential coverage analysis...\n", sep=" "))
		
		##Correct for test selection if necessary
		if((nMSets1<3 | nMSets2<3) & diff.method=="ttest"){
			stop("Method 'ttest' is not valid for less than 3 replicates per group. Method 'edgeR' can be applied in this case.")
		}
	
		if(diff.method=="edgeR"){
			##Extract number of reads per sample
			if(!is.null(MSet1)){
				n.r.M1 = NULL
				for(i in 1:nMSets1){
					n.r.M1 = c(n.r.M1, number_regions(MSet1[[i]]))		
				}
			}	
			if(!is.null(MSet2)){
		 		n.r.M2 = NULL
		 		for(i in 1:nMSets2){
					n.r.M2 = c(n.r.M2, number_regions(MSet2[[i]]))	
	  	 		}
			}			
			diff.results.list = MEDIPS.diffMeth(base=base, values=counts.medip, diff.method="edgeR", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, n.r.M1=n.r.M1, n.r.M2=n.r.M2, MeDIP=MeDIP, minRowSum=minRowSum)
		}
		else if(diff.method=="ttest"){
			if(type=="rpkm"){
				diff.results.list = MEDIPS.diffMeth(base=base, values=rpkm.medip, diff.method="ttest", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, MeDIP=MeDIP, minRowSum=minRowSum)
			}
			else if(type=="rms"){
				if(MeDIP){
					diff.results.list = MEDIPS.diffMeth(base=base, values=rms, diff.method="ttest", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, MeDIP=MeDIP, minRowSum=minRowSum)
				}
				else{
					stop("Invalid specification for parameter type because parameter MeDIP is FALSE (no rms values have been calculated).")
				}
			}
			else{
				stop("Unknown specification for parameter type.")
			}
		}
		else{stop("Selected method for calculating differential coverage not supported")}
		
		cat("Please note, log2 ratios are reported as log2(MSet1/MSet2).\n")
		diff.results = diff.results.list$diff.results
		diff.index = diff.results.list$diff.index
				
		rm(diff.results.list)
		gc()
	
	}
	else{
		cat("No differential coverage will be calculated- only one group of MEDIPS SETs given.\n")
	}
	
	##If two groups of INPUT SETs are given 
	##calculate CNV on mean per set
	##################################
	if(CNV){
		if(!is.null(ISet1) & !is.null(ISet2)){
			cat(paste("CNV analysis...\n", sep=" "))		
			cnv.combined = MEDIPS.cnv(base=base, rpkm.input=rpkm.input, nISets1=nISets1, nISets2=nISets2)
		}
		else{
			cat("Cannot perform CNV analysis- please specify two groups of INPUT SETs!\n")
		}
	}
	
	##Create results table	
	##################################
	cat(paste("Creating results table...\n", sep=" "))
	if(!is.null(counts.medip)){
		if(MeDIP){
			results = data.frame(base, counts.medip, rpkm.medip, rms, prob, stringsAsFactors=F)
		}
		else{
			results = data.frame(base, counts.medip, rpkm.medip, stringsAsFactors=F)
		}
	}
	if(!is.null(counts.input)){
		if(!is.null(counts.medip))
		{
			results = data.frame(results, counts.input, rpkm.input, stringsAsFactors=F)
		}
		else{
			results = data.frame(base, counts.input, rpkm.input, stringsAsFactors=F)
		}
	}
	
	##Add mean counts, rpkm, probs columns
	
	#if(nMSets1>1){
	set1idx=1:(nMSets1)
	#counts.mean.C=apply(FUN=mean,X=counts.medip[,set1idx,drop=F],MARGIN=1)
	#rpkm.mean.C=apply(FUN=mean,X=rpkm.medip[,set1idx,drop=F],MARGIN=1)	
	counts.mean.C=numeric(dim(counts.medip)[1])
	rpkm.mean.C=numeric(dim(rpkm.medip)[1])
	for (i in set1idx){
		counts.mean.C=counts.mean.C+counts.medip[,i]
		rpkm.mean.C=rpkm.mean.C+rpkm.medip[,i]
	}
	counts.mean.C=counts.mean.C/nMSets1
	rpkm.mean.C=rpkm.mean.C/nMSets1
	if(MeDIP){
		#rms.mean.C = apply(FUN=mean,X=rms[,set1idx,drop=F],MARGIN=1)	
		#prob.mean.C = apply(FUN=mean,X=prob[,set1idx,drop=F],MARGIN=1)	
		rms.mean.C =numeric(dim(rms)[1])
		prob.mean.C =numeric(dim(prob)[1])
		for (i in set1idx){
			rms.mean.C=rms.mean.C+rms[,i]
			prob.mean.C=prob.mean.C+prob[,i]
		}
		rms.mean.C=rms.mean.C/nMSets1
		prob.mean.C=prob.mean.C/nMSets1
		results = data.frame(results, MSets1.counts.mean=counts.mean.C, MSets1.rpkm.mean=rpkm.mean.C, MSets1.rms.mean=rms.mean.C, MSets1.prob.mean=prob.mean.C, stringsAsFactors=F)
		rm(counts.mean.C,rpkm.mean.C,set1idx,rms.mean.C,prob.mean.C)
	}
	else{
		results = data.frame(results, MSets1.counts.mean=counts.mean.C, MSets1.rpkm.mean=rpkm.mean.C, stringsAsFactors=F)
		rm(counts.mean.C,rpkm.mean.C,set1idx)
	}
	#}

	if(nMSets2>0){
		set2idx=(nMSets1+1):(nMSets1+nMSets2)
		#counts.mean.T=apply(FUN=mean,X=counts.medip[,set2idx,drop=F],MARGIN=1)
		#rpkm.mean.T=apply(FUN=mean,X=rpkm.medip[,set2idx,drop=F],MARGIN=1)	
		counts.mean.T=numeric(dim(counts.medip)[1])
		rpkm.mean.T=numeric(dim(rpkm.medip)[1])
		for (i in set2idx){
			counts.mean.T=counts.mean.T+counts.medip[,i]
			rpkm.mean.T=rpkm.mean.T+rpkm.medip[,i]
		}	
		counts.mean.T=counts.mean.T/nMSets2
		rpkm.mean.T=rpkm.mean.T/nMSets2

		if(MeDIP){
			#rms.mean.T = apply(FUN=mean,X=rms[,set2idx,drop=F],MARGIN=1)	
			#prob.mean.T = apply(FUN=mean,X=prob[,set2idx,drop=F],MARGIN=1)	
			rms.mean.T =numeric(dim(rms)[1])
			prob.mean.T =numeric(dim(prob)[1])
			for (i in set2idx){
				rms.mean.T=rms.mean.T+rms[,i]
				prob.mean.T=prob.mean.T+prob[,i]
			}
			rms.mean.T=rms.mean.T/nMSets2
			prob.mean.T=prob.mean.T/nMSets2
			results = data.frame(results, MSets2.counts.mean=counts.mean.T, MSets2.rpkm.mean=rpkm.mean.T, MSets2.rms.mean=rms.mean.T, MSets2.prob.mean=prob.mean.T, stringsAsFactors=F)
			rm(counts.mean.T,rpkm.mean.T,set2idx,rms.mean.T,prob.mean.T)
		}
		else{
			results = data.frame(results, MSets2.counts.mean=counts.mean.T, MSets2.rpkm.mean=rpkm.mean.T, stringsAsFactors=F)
			rm(counts.mean.T,rpkm.mean.T,set2idx)
		}
	}

	if(nISets1>1){
		setI1idx=1:(nISets1)
		#counts.input.mean.C = apply(FUN=mean,X=counts.input[,setI1idx,drop=F],MARGIN=1)
		#rpkm.input.mean.C = apply(FUN=mean,X=rpkm.input[,setI1idx,drop=F],MARGIN=1)
		counts.input.mean.C = counts.input[,setI1idx[1]]
		rpkm.input.mean.C = rpkm.input[,setI1idx[1]]
		for (i in setI1idx[-1]){
		  counts.input.mean.C =counts.input.mean.C+counts.input[,i]
		  rpkm.input.mean.C =rpkm.input.mean.C+rpkm.input[,i]
		}
		counts.input.mean.C =counts.input.mean.C/nISets1
		rpkm.input.mean.C =rpkm.input.mean.C/nISets1
		results = data.frame(results, ISets1.counts.mean=counts.input.mean.C, ISets1.rpkm.mean=rpkm.input.mean.C, stringsAsFactors=F)
		rm(counts.input.mean.C,rpkm.input.mean.C,setI1idx)
	}
	if(nISets2>1){
		setI2idx=(nISets1+1):(nISets1+nISets2)
		#counts.input.mean.T = apply(FUN=mean,X=counts.input[,setI2idx,drop=F],MARGIN=1)
		#rpkm.input.mean.T = apply(FUN=mean,X=rpkm.input[,setI2idx,drop=F],MARGIN=1)
		counts.input.mean.T = counts.input[,setI2idx[1]]
		rpkm.input.mean.T = rpkm.input[,setI2idx[1]]
		for (i in setI2idx[-1]){
		  counts.input.mean.T =counts.input.mean.T+counts.input[,i]
		  rpkm.input.mean.T =rpkm.input.mean.T+rpkm.input[,i]
		}
		counts.input.mean.T =counts.input.mean.T/nISets2
		rpkm.input.mean.T =rpkm.input.mean.T/nISets2

		results = data.frame(results, ISets2.counts.mean=counts.input.mean.T, ISets2.rpkm.mean=rpkm.input.mean.T, stringsAsFactors=F)
		rm(counts.input.mean.T,rpkm.input.mean.T,setI2idx)
	}
	
	if(MeDIP){rm(base, counts.medip, rpkm.medip, rms, prob, counts.input, rpkm.input)}
	else{rm(base, counts.medip, rpkm.medip, counts.input, rpkm.input)}
	gc()
	
	##Add diff.meth results
	if(nMSets1!=0 & nMSets2!=0){
		cat(paste("Adding differential coverage results...\n", sep=" "))
		dummy.results = matrix(ncol=ncol(diff.results), nrow=nrow(results))
		
		if(diff.method=="edgeR"){
			c.names = colnames(diff.results)
			diff.results <- matrix(unlist(diff.results), ncol=ncol(diff.results), byrow=FALSE)
			colnames(diff.results)=c.names
			rm(c.names)
		}
		
		dummy.results[diff.index,] = diff.results
		colnames(dummy.results)=colnames(diff.results)		
		results = data.frame(results, dummy.results, stringsAsFactors=F)
			
		rm(diff.results, dummy.results, diff.index)
		gc()
	}
	
	##Add CNV results
	if(!is.null(ISet1) & !is.null(ISet2)){
		if(CNV){
		cat(paste("Adding CNV results...\n", sep=" "))
		dummy.results = matrix(ncol=1, nrow=(nrow(results)))
		for(i in 1:nrow(cnv.combined)){
			dummy.results[cnv.combined[i,1]:cnv.combined[i,2]] = cnv.combined[i,3]		
		}
		colnames(dummy.results)="CNV.log2.ratio"
		results = data.frame(results, dummy.results, stringsAsFactors=F)
		
		rm(dummy.results)	
		gc()
	}
	}
	
	rownames(results) = seq(1, nrow(results))
    	gc()
	return(results)
		
}

