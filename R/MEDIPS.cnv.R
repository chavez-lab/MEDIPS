##########################################################################
##Function calcuates CNV via DNAcopy for two groups of INPUT SETs
##########################################################################
##Input:	Two groups of INPUT SETs
##Param:	base, rpkm.input, nISets1, nISets2
##Output:	results of CNV analysis
##Requires:	DNAcopy
##Modified:	11/22/2011
##Author:	Lukas Chavez, Joern Dietrich

MEDIPS.cnv = function(base=base, rpkm.input=NULL, nISets1=NULL, nISets2=NULL)
{
	if(nISets1>1){
		control.input.mean = rowMeans(rpkm.input[,1:nISets1])
	}
	else{
		control.input.mean = rpkm.input[,1]
	}
	
	if(nISets2>1){
                treat.input.mean = rowMeans(rpkm.input[,(nISets1+1):ncol(rpkm.input)])
        }
        else{
		treat.input.mean = rpkm.input[,nISets1+1]
       	}


	log2ratios.inputs = log2(control.input.mean/treat.input.mean)
	cnv.ID = "ISets.mean"
	
	nperm=10000
	alpha=0.01 
	max.ones=floor(nperm*alpha)+1
	default.DNAcopy.bdry=DNAcopy::getbdry(nperm=10000, eta=0.05,max.ones=max.ones)
  	CNA.object <- DNAcopy::CNA(log2ratios.inputs, base[,1], floor((base[,2]+base[,3])/2), data.type="logratio", sampleid=cnv.ID, presorted=TRUE)
	
	smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object)
	
	segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, verbose=2,sbdry = default.DNAcopy.bdry)				
	
	cnv.combined = cbind(segment.smoothed.CNA.object$segRows$startRow,
						 segment.smoothed.CNA.object$segRows$endRow, 
						 segment.smoothed.CNA.object$output$seg.mean)
	rm(CNA.object, smoothed.CNA.object, segment.smoothed.CNA.object, control.input.mean, treat.input.mean, log2ratios.inputs)
	gc()
	
	return(cnv.combined)
}

##########################################################################
##Function to add results from CNV analysis 
##########################################################################
##Input:	Result matrix and input sets 
##Param:	ISet1,ISet2,results,uniq,chr,cnv.Frame
##Output:	List of annotation matrix 
##Requires:	BSgenome,DNAcopy
##Modified:	11/22/2011
##Author:	Joern Dietrich
MEDIPS.addCNV<-function(ISet1, ISet2, results, cnv.Frame=1000){

	nISets1 = length(ISet1)	
	nISets2 = length(ISet2)	
	
	if(!is.list(ISet1)) ISet1=c(ISet1)
	if(!is.list(ISet2)) ISet2=c(ISet2)
	
	cnv.rpkm.input=NULL
	
	for(i in 1:nISets1){
		file=paste(path_name(ISet1[[i]]),sample_name(ISet1[[i]]),sep="/");
		tmp.cnv = MEDIPS.createSet(file=file, 
				BSgenome=genome_name(ISet1[[i]]), 
				chr.select=chr_names(ISet1[[i]]),
				extend=extend(ISet1[[i]]),
				shift=shifted(ISet1[[i]]), 
				uniq=uniq(ISet1[[i]]), 
				window_size=cnv.Frame)
		
		cnv.rpkm.input = cbind(cnv.rpkm.input, ((genome_count(tmp.cnv)*10^9)/(cnv.Frame*number_regions(tmp.cnv))))	
	}
	
	for(i in nISets2){
		file=paste(path_name(ISet2[[i]]),sample_name(ISet2[[i]]),sep="/");
		tmp.cnv = MEDIPS.createSet(file=file, 
				BSgenome=genome_name(ISet2[[i]]),
				chr.select=chr_names(ISet2[[i]]), 
				extend=extend(ISet2[[i]]), 
				shift=shifted(ISet2[[i]]),
				uniq=uniq(ISet2[[i]]), 
				window_size=cnv.Frame)
		
		cnv.rpkm.input = cbind(cnv.rpkm.input, ((genome_count(tmp.cnv)*10^9)/(cnv.Frame*number_regions(tmp.cnv))))	
	}
							           							           
	no_chr_windows = ceiling(chr_lengths(tmp.cnv)/cnv.Frame)
	supersize_chr = cumsum(no_chr_windows)
	cvn.GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr_names(tmp.cnv), chr_lengths(tmp.cnv), cnv.Frame)
	base = data.frame(chr=as.vector(seqnames(cvn.GRanges.genome)), start=start(cvn.GRanges.genome), stop=end(cvn.GRanges.genome), stringsAsFactors=F)
	
	cat(paste("CNV analysis...\n", sep=" "))
	cnv.combined = MEDIPS.cnv(base=base, rpkm.input=cnv.rpkm.input, nISets1=nISets1, nISets2=nISets2)
	
	dummy.results = matrix(ncol=1, nrow=(nrow(base)))
	for(i in 1:nrow(cnv.combined)){
		dummy.results[cnv.combined[i,1]:cnv.combined[i,2]] = cnv.combined[i,3]		
	}
	colnames(dummy.results)="CNV"
    	dummy.results=cbind(dummy.results, base)
	
	results=MEDIPS.setAnnotation(results, dummy.results, cnv=T)
	
	rm(dummy.results)	
	gc()
	
	return(results)
}



