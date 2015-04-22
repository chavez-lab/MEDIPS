##########################################################################
##Selects putative DMRs from result table returned by MEDIPS.meth()
## w.r.t. to specified parameters
##########################################################################
##Input:	result table returned by MEDIPS.meth (or by MEDIPS.selectROIs)
##Param:	
##Output:	row-wise subset of input
##Modified:	11/28/2011
##Author:	Lukas Chavez

MEDIPS.selectSig = function(results=NULL, p.value=0.01, adj=T, ratio=NULL, bg.counts=NULL, CNV=F){
	
	if(is.null(results)){stop("Must specify a result table.")}		
	cat(paste("Total number of windows: ", nrow(results), "\n", sep=""))		

	##Background counts- preparation
        ################################
        if(!is.null(bg.counts)){
                if(is.numeric(bg.counts)){
                        bg.t = bg.counts
                }
                else{
                     	column.bg = grep(bg.counts, colnames(results))
                        if(length(column.bg)==0){
                                stop("Specified background column not available.")
                        }
                        bg.t = quantile(results[, column.bg], probs=0.95)
                }
        }

	##Filter for windows tested for differential  methylation
	##########################################################	
	results=results[!is.na(results[,grep("p.value", colnames(results))[1]]),]
	cat(paste("Number of windows tested for differential methylation: ", nrow(results), "\n", sep=""))
	
	
	#Filter for p.values
	#######################
	if(adj){value_pvalue = grep("adj.p.value", colnames(results))}
	else{value_pvalue = grep("p.value", colnames(results))[1]} 
		
	results=results[results[,value_pvalue]<=p.value,]
	
	if(adj){cat(paste("Remaining number of windows with adjusted p.value<=", p.value, ": ", nrow(results), "\n", sep=""))}
	else{cat(paste("Remaining number of windows with p.value<=", p.value, ": ", nrow(results), "\n", sep=""))}

	##Ratio and CNV specific filter
	##############################
	##Filter for ratio
	if(!is.null(ratio)& CNV==F){
		
		column_ratio = grep("score.log2.ratio", colnames(results))
		if(length(column_ratio)==0){
			column_ratio = grep("edgeR.logFC", colnames(results))
		}

		results=results[results[,column_ratio]>=log2(ratio) | results[,column_ratio]<=(log2(1/ratio)),]
		cat(paste("Remaining number of windows with ratio >=",ratio," (or <=", round(1/ratio, digits=2), ", respectively): ", nrow(results), "\n", sep=""))		
	}	
	if(!is.null(ratio) & CNV==T){
		
		column_ratio = grep("score.log2.ratio", colnames(results))
                if(length(column_ratio)==0){
                        column_ratio = grep("edgeR.logFC", colnames(results))   
               	}
		
		column_CNVratio = grep("CNV.log2.ratio", colnames(results))
		if(length(column_CNVratio)==0){
			stop("No CNV results available.")
                }
		else if(length(column_CNVratio)>1){
			column_CNVratio = column_CNVratio[1]
		}
		results[is.na(results[,column_CNVratio]),column_CNVratio]=0
		results=results[abs(results[,column_ratio]-results[,column_CNVratio])>=abs(log2(ratio)),]
		cat(paste("Remaining number of windows with CNV corrected ratio >=",ratio," (or <=", round(1/ratio, digits=2), ", respectively): ", nrow(results), "\n", sep=""))
	}
	if(CNV==T & is.null(ratio)){
		stop("Please specify a ratio when setting CNV=T.")
	}

	##Background counts 
        ########################################
        if(!is.null(bg.counts)){
                results=results[results$MSets1.counts.mean>=bg.t | results$MSets2.counts.mean>=bg.t,]
                cat(paste("Remaining number of windows where the mean count of at least one group is >=", round(bg.t, digits=2), ": ", nrow(results),"\n", sep=""))
        }

	gc()
	return(as.data.frame(results,  stringsAsFactors=F))
}


