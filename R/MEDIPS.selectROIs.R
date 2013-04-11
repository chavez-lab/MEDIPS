##########################################################################
##Function selects subsets of a results table returned by MEDIPS.meth w.r.s
##given ROIs and (optional) summarized methylation values of selected windows
##########################################################################
##Input:	results table from MEDIPS.meth, genomic coordinates of ROIs as list or as GRange object
##		In the latter case the Id must be in the first column of values
##Param:	results, rois, columns, summarize
##Output:	subset of results, or mean methylation values for each ROI and for selected columns 
##Modified:	11/24/2011
##Author:	Lukas Chavez

MEDIPS.selectROIs=function(results=NULL, rois=NULL, columns=NULL, summarize=F){
	
	##Get data column ids to summarize for ROIs
	#column.ids=NULL
	#for(i in 1:length(columns)){
	#	column.ids = c(column.ids, grep(columns[i], colnames(results)))
	#}
	if(!is.null(columns))
		column.ids = which(colnames(results) %in% columns)
	else
		column.ids = 4:dim(results)[2]

	

	##Convert results to GRange object
	results.GRange = GRanges(seqnames=results[,1], ranges=IRanges(start=as.numeric(results[,2]), end=as.numeric(results[,3])))		
	elementMetadata(results.GRange) = results[, c(column.ids)]
	colnames(elementMetadata(results.GRange)) = colnames(results)[column.ids]
	
	if(class(rois)!="GRanges"){
		##Convert rois to GRange object
		rois.Grange = GRanges(seqnames=rois[,1], ranges=IRanges(start=as.numeric(rois[,2]), end=as.numeric(rois[,3])), ids=rois[,4])		
	}	
	
	if(summarize){
		##For each roi, select according result frames and 
		##calculate mean for all specified columns
		###################################################
		m=IRanges::as.matrix(findOverlaps(rois.Grange, results.GRange))
		ind=m[,2]
		g=m[,1]
		rm(m)		
		gc()
   		mean.rois.data=apply(IRanges::as.data.frame(values(results.GRange)),2,function(x){
       			l=split(x[ind],g)
       			return(unlist(lapply(l, mean, na.rm=T)))
   		})		
		##Convert Grange rois result object to base
		g = unique(g)
		base = data.frame(chr=as.character(as.vector(seqnames(rois.Grange[g]))), start=start(rois.Grange[g]), end=end(rois.Grange[g]), stringsAsFactors=F)
		ids = as.character(values(rois.Grange[g])[,1])		
		out=cbind(base, mean.rois.data)
		rownames(out)=ids
		
		return(out)		
	}
	else{
		##Preselect from results only ranges that overlap with rois 
		results.rois = subsetByOverlaps(results.GRange, rois.Grange)
	
		##Convert Grange result object to base and data matrices	
		base = data.frame(chr=as.character(as.vector(seqnames(results.rois))), start=start(results.rois), end=end(results.rois))	
		results.rois.data = IRanges::as.data.frame(results.rois)[,6:ncol(IRanges::as.data.frame(results.rois))]		
		colnames(results.rois.data)=colnames(results)[column.ids]
		return(cbind(base, results.rois.data))		
	}	
}
    	   	
	

