##########################################################################
##Function selects subsets of a results table returned by MEDIPS.meth w.r.s
##given ROIs and (optional) summarized methylation values of selected windows
##########################################################################
##Input:	results table from MEDIPS.meth, genomic coordinates of ROIs as list or as GRange object
##		In the latter case the Id must be in the first column of values
##Param:	results, rois, columns, summarize
##Output:	subset of results, or mean methylation values for each ROI and for selected columns 
##Modified:	05/15/2013
##Author:	Lukas Chavez, Matthias Lienhard

MEDIPS.selectROIs=function(results=NULL, rois=NULL, columns=NULL, summarize=NULL){
	
	if (! (is.null(summarize) || summarize=="avg" || summarize=="minP")){
		stop("summarize must be \"avg\", \"minP\" or NULL")
	} 
		
	if(!is.null(columns)){
		column.ids = which(colnames(results) %in% columns)
		}else{
		column.ids = 4:dim(results)[2]
		}

	##Convert results to GRange object
	results.GRange = GRanges(seqnames=results[,1], ranges=IRanges(start=as.numeric(results[,2]), end=as.numeric(results[,3])))		
	elementMetadata(results.GRange) = results[, c(column.ids)]
	colnames(elementMetadata(results.GRange)) = colnames(results)[column.ids]
	
	if(class(rois)!="GRanges"){
		##Convert rois to GRange object
		rois.Grange = GRanges(seqnames=rois[,1], ranges=IRanges(start=as.numeric(rois[,2]), end=as.numeric(rois[,3])), ids=rois[,4])		
	}else{
		rois.Grange=rois
		if(is.null(rois.Grange$ids)) stop("Granges object for ROIs needs to contain \"ids\"")
	}	
	if(is.null(summarize)){
		##Preselect from results only ranges that overlap with rois 
		m=IRanges::as.matrix(findOverlaps(rois.Grange, results.GRange))

		##Convert Grange result object to base and data matrices	
		results.rois.data = results[m[,2],c(1,2,3,column.ids)]
		results.rois.data$ROI=rois.Grange$ids[m[,1]]
		return(results.rois.data)		
	}else if(summarize=="avg"){
		##For each roi, select according result frames and 
		##calculate mean for all specified columns
		###################################################
		m=IRanges::as.matrix(findOverlaps(rois.Grange, results.GRange))
   		mean.rois.data=apply(IRanges::as.data.frame(values(results.GRange)),2,function(x){
       			l=split(x[m[,2]],m[,1])
       			return(unlist(lapply(l,mean)))
   		})		
		##Convert Grange rois result object to base
		g = unique(m[,1])
		base = data.frame(chr=as.character(as.vector(seqnames(rois.Grange[g]))), start=start(rois.Grange[g]), end=end(rois.Grange[g]), stringsAsFactors=F)
		ids = rois.Grange$ids[g]
		if(dim(base)[1]==1){
			avgROI=data.frame(c(base, mean.rois.data))
		}else{avgROI=cbind(base, mean.rois.data)}
		avgROI$ROI=ids		
		return(avgROI)
	}else if(summarize =="minP"){
		#find ROIs in data
		m=IRanges::as.matrix(findOverlaps(rois.Grange, results.GRange))
		#find col containing the pvalue
		pval_idx=grep("p.value", names(results)[c(1,2,3,column.ids)])[1]
	        if(is.na(pval_idx)){
			stop("Error: no p.value column selected while summarize is \"minP\"")
        	}
		#group windows by roi
                roi_list=split(results[m[,2],c(1,2,3,column.ids)],m[,1])
		#set the name of the rois
		names(roi_list)=rois.Grange$ids[as.numeric(names(roi_list))]
		#find rois containing windows with pvalue
		nonEmpty=lapply(X=roi_list, FUN=nrow)>0 &
		unlist(lapply(X=roi_list,FUN=function(x){return(any(!is.na(x[,pval_idx])))}))
		#for each ROI: select most significant window 
       		minPvalROI=unsplit(lapply(X=roi_list[nonEmpty], FUN=function(x){return(as.vector(x[which.min(x[,pval_idx]),]))} ),1:sum(nonEmpty))
		#set name of rois in result table
		minPvalROI$ROI=names(roi_list)[nonEmpty]
		return(minPvalROI)
	}
}
    	   	
	

