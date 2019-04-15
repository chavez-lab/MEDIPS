##########################################################################
##Function provides two moudiles for calculating differential methylation
##########################################################################
##Input:	genomic coordinates and data for groups of MEDIPS SETs
##Param:	base, values, diff.method, nMSets1, nMSets2, p.adj, n.r.M1, n.r.M2, p.adj, MeDIP, minRowSum
##Output:	results of differential methylation analysis plus index in genome wide table
##Requires:	edgeR
##Modified:	1/9/2015
##Author:	Lukas Chavez, Joern Dietrich

MEDIPS.diffMeth = function(base = NULL, values=NULL, diff.method="ttest", nMSets1=NULL, nMSets2=NULL, n.r.M1=n.r.M1, n.r.M2=n.r.M2, p.adj="bonferroni", MeDIP, minRowSum=10, diffnorm="tmm")
{
	##edgeR##
	#########
	if(diff.method=="edgeR"){

		##Extract non-zero MeDIP count windows
		message("Extracting count windows with at least ", minRowSum, " reads...", appendLF=T)
		filter= rowSums(values)>=minRowSum

		##Extract non-zero coupling factor rows
		if(MeDIP){
			message("Extracting non-zero coupling factor windows...", appendLF=T)
			filter=filter & base[,4]!=0
		}

		message("Execute edgeR for count data of ", sum(filter), " windows...", appendLF=T)
		message("(Neglecting parameter 'type')", appendLF=T)

		message("Creating a DGEList object...", appendLF=T)
		edRObj.group=c(rep(1, nMSets1), rep(2, nMSets2))
		edRObj.length=c(n.r.M1, n.r.M2)
		#d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group, lib.size=edRObj.length)
		d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group)

		if(diffnorm=="tmm"){
			message("Apply trimmed mean of M-values (TMM) for library sizes normalization...", appendLF=T)
			d=edgeR::calcNormFactors(d)
		}
		if(diffnorm=="quantile" | diffnorm=="none"){
			message("Skipping trimmed mean of M-values (TMM) library size normalization...", appendLF=T)
			d=edgeR::calcNormFactors(d, method="none")
		}
	    if(diffnorm!="tmm" & diffnorm!="quantile" & diffnorm!="none"){
	    	stop("diffnorm method unknown.")
	    }

		if(nMSets1!=1 | nMSets2!=1){
			message("Estimating common dispersion...", appendLF=T)
			d=edgeR::estimateCommonDisp(d)
			message("Estimating tagwise dispersion...", appendLF=T)
			d=edgeR::estimateTagwiseDisp(d)
			message("Calculating differential coverage...", appendLF=T)
			de.com=edgeR::exactTest(d,pair=c("2","1"))
		}
		else{
			warning("There is no replication, setting dispersion to bcv^2 where bcv=0.01.\n")
			warning("Please consider http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf section 2.9.\n")
			bcv = 0.01
			message("Calculating differential coverage...", appendLF=T)
			de.com = suppressWarnings(edgeR::exactTest(d, dispersion=bcv^2, pair=c("2","1")))
		}

		##Adjusting p.values for multiple testing
		message("Adjusting p.values for multiple testing...", appendLF=T)
		colnames(de.com$table) = c("edgeR.logFC", "edgeR.logCPM", "edgeR.p.value")
		diff.results = cbind(de.com$table, edgeR.adj.p.value=p.adjust(de.com$table$edgeR.p.value, p.adj))

		rm(de.com, edRObj.group, edRObj.length, d)
		gc()
	}

	##ttest/Score##
	#########
	else{
		##Extract non-zero MeDIP rows
		message("Extracting count windows with at least ",minRowSum, " reads...", appendLF=T)

		filter= rowSums(values)>=minRowSum

		##Extract non-zero coupling factor windows
		if(MeDIP){
			message("Extracting non-zero coupling factor windows...", appendLF=T)
			filter=filter & base[,4]!=0
		}

		#filter NA
		#if there is a na in a row, this row is sorted out?!
		filter=filter & !is.na(rowSums(values))


		message("Calculating score for ", sum(filter), " windows...", appendLF=T)
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
		message("Adjusting p.values for multiple testing...", appendLF=T)

		diff.results = cbind(score.log2.ratio=log2(ratio), score.p.value=t.test.p.value, score.adj.p.value=p.adjust(t.test.p.value, p.adj), score=score)

		rm(const, ratio, t.test.p.value, score)

	}
	return(list(diff.results=diff.results, diff.index=which(filter)))
}
