##########
#Function to plot the results of the sequence coverage analysis.
##########
##Input:	A seqCoverage result object.
##Param:	seqCoverageObj, main (characrter), type (pie | hist), cov.level
##Output:	Plot
##Modified:	11/10/2011
##Author:	Lukas Chavez

MEDIPS.plotSeqCoverage <-
function(seqCoverageObj=NULL, main=NULL, type="pie", cov.level = c(0,1,2,3,4,5), t="Inf"){
	
	if(is.null(seqCoverageObj)){stop("Must specify a coverage object.")}

	cov.res=	seqCoverageObj$cov.res
	numberPos=	length(cov.res)	
	pattern=	seqCoverageObj$pattern
	numberReads=	seqCoverageObj$numberReads
	numberReadsWO=	seqCoverageObj$numberReadsWO
	
	if(is.null(main)){main = paste("Total number of ", pattern, "'s: ", numberPos, sep="") }
	
	sub = paste(numberReadsWO, " of ", numberReads, " reads (", round((numberReadsWO/numberReads*100), digits=2) , "%) do not cover a pattern", sep="")
 	
	if(type=="pie"){
	
		cat("Creating summary...\n")
		results = NULL
		for(i in 1:length(cov.level)){
			if(i==1){
				results = c(results, length(cov.res[cov.res<=cov.level[i]]))
			}
			else{
				results = c(results, length(cov.res[cov.res>cov.level[i-1] & cov.res<=cov.level[i]]))
			}
		}
		results = c(results, length(cov.res[cov.res>cov.level[length(cov.level)]]))	
		
		#Set the labels for the output plot
		labels = NULL
		for(i in 1:length(cov.level)){
			if(i==1){
				labels = c(labels, paste("<=", cov.level[i], "x (", round((results[i]/numberPos)*100, digits=2), "%)", sep=""))
			}
			else{
				labels = c(labels, paste(cov.level[i-1]+1, "-", cov.level[i], "x (", round((results[i]/numberPos)*100, digits=2), "%)", sep=""))
		 	}
		}
		labels = c(labels, paste(">", cov.level[length(cov.level)], "x (", round((results[length(results)]/numberPos)*100, digits=2), "%)", sep=""))
			
	
		#Set the colors for the output plot
		col = NULL
		for(i in 1:(length(cov.level))){	
			if(i==1){col="red"}
			if(i==2){col=c(col,"green")}
			if(i==3){col=c(col,"blue")}
			if(i==4){col=c(col,"orange")}
			if(i==5){col=c(col,"brown")}
			if(i>5){col=c(col,"darkgrey")}
		}
		##Plot	
		pie(
			results, 
			labels=labels, 
			col=col, 
			main=main,
			sub=sub
		)
	}
	else if(type=="hist"){
		hist(cov.res[cov.res<=t], 
			main=main,
			sub=sub
		)
	}
	else{stop("Plot type not supported.")}
	
}
