##########################################################################
##Function provides two moudiles for calculating differential methylation
##########################################################################
##Input:	genomic coordinates and data for groups of MEDIPS SETs
##Param:	base, values, diff.method, nMSets1, nMSets2, p.adj
##Output:	results of differential methylation analysis plus index in genome wide table
##Requires:	edgeR
##Modified:	11/22/2011
##Author:	Lukas Chavez, Joern Dietrich

MEDIPS.diffMeth = function(base = NULL, values=NULL, diff.method="ttest", nMSets1=NULL, nMSets2=NULL, n.r.M1=n.r.M1, n.r.M2=n.r.M2, p.adj="bonferroni", MeDIP, minRowSum=12)
{
	##edgeR##
	#########
	if(diff.method=="edgeR"){		
			
		##Extract non-zero MeDIP count windows
		cat(paste("Extracting count windows with at least",minRowSum," reads...\n", sep=" "))
		filter= rowSums(values)>=minRowSum
		##Extract non-zero coupling factor rows

		if(MeDIP){
			cat(paste("Extracting non-zero coupling factor windows...\n", sep=" "))
			filter=filter & base[,4]!=0
		}

		cat(paste("Execute edgeR for count data of", sum(filter), "windows...\n", sep=" "))
		cat("(Neglecting parameter 'type')\n")
							
		edRObj.group=c(rep(1, nMSets1), rep(2, nMSets2))
		edRObj.length=c(n.r.M1, n.r.M2)
		d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group, lib.size=edRObj.length)

		#rm(values)
		#gc()

		d=edgeR::calcNormFactors(d)

		if(nMSets1!=1 | nMSets2!=1){
			d=edgeR::estimateCommonDisp(d)
			de.com=edgeR::exactTest(d,pair=c("2","1"))
		}
		else{
			cat("There is no replication, setting dispersion to 0.01.\n")
			cat("Please consider http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf section 2.9.\n")
			
			de.com = suppressWarnings(edgeR::exactTest(d, dispersion=0.01, pair=c("2","1")))
		}
		
		##Adjusting p.values for multiple testing
		cat(paste("Adjusting p.values for multiple testing...\n", sep=" "))
		colnames(de.com$table) = c("edgeR.logFC", "edgeR.logCPM", "edgeR.p.value")		
		diff.results = cbind(de.com$table, edgeR.adj.p.value=p.adjust(de.com$table$edgeR.p.value, p.adj)) 
		#diff.index = match(rownames(base.sub), rownames(base))
		
		rm(de.com, edRObj.group, edRObj.length, d)
		gc()
	}
		
	##ttest/Score##
	#########
	else{			
		##Extract non-zero MeDIP rows
		cat(paste("Extracting count windows with at least",minRowSum," reads...\n", sep=" "))

		filter= rowSums(values)>=minRowSum

		##Extract non-zero coupling factor windows
		if(MeDIP){
			cat(paste("Extracting non-zero coupling factor windows...\n", sep=" "))
			filter=filter & base[,4]!=0
		}
	
		#filter NA
		#if there is a na in a row, this row is sorted out?!
		filter=filter & !is.na(rowSums(values))

			
		cat(paste("Calculating score for", sum(filter), "windows...\n", sep=" "))
		ms1=1:nMSets1
		ms2=(nMSets1+1):(nMSets1+nMSets2)
		##Calculate ratios##
		####################
		ratio=rowSums(values[filter,ms1]+0.1)/rowSums(values[filter,ms2,drop=F]+0.1)*(nMSets2/nMSets1)
	
		##Calculate p.values##
		######################
		##Check for constant entries
		## Filter out constant entries, as they have no sd

		const = apply(X=values[filter,ms1,drop=F],MARGIN=1,FUN=min) - apply(X=values[filter,ms1,drop=F],MARGIN=1,FUN=max) 	+
			apply(X=values[filter,ms2,drop=F],MARGIN=1,FUN=min) - apply(X=values[filter,ms2,drop=F],MARGIN=1,FUN=max) 	== 0

		t.test.p.value = rep(NA, sum(filter))				

		t.test.p.value[!const]=matTtest(values[filter,][!const,],groups=c(rep(1,nMSets1),rep(2,nMSets2)))$p.value

		##Calculate the final score##
		#############################
		score = (-log10(t.test.p.value)*10)*log(ratio)
		
		##Adjusting p.values for multiple testing
		cat(paste("Adjusting p.values for multiple testing...\n", sep=" "))	
			
		diff.results = cbind(score.log2.ratio=log2(ratio), score.p.value=t.test.p.value, score.adj.p.value=p.adjust(t.test.p.value, p.adj), score=score)
			
		rm(const, ratio, t.test.p.value, score)
		
	}		
	return(list(diff.results=diff.results, diff.index=which(filter)))
}

