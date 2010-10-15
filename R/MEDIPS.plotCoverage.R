##########
#Function to plot the results of the coverage analysis.
##########
MEDIPS.plotCoverage <-
function(coverageObj=NULL){
	
	if(is.null(coverageObj)){stop("Must specify a coverage object.")}
	coverage_matrix=coverageObj$matrix
	maxPos=coverageObj$maxPos
	pattern=coverageObj$pattern
	
	#Set the colors for the output plot
	for(i in 1:((length(coverage_matrix[1,]))-1)){	
		if(i==1){colors="red"}
		if(i==2){colors=c(colors,"green")}
		if(i==3){colors=c(colors,"blue")}
		if(i==4){colors=c(colors,"orange")}
		if(i==5){colors=c(colors,"brown")}
		if(i>5){colors=c(colors,"black")}
	}
	
	##Plot	
	legend=paste("min. ", coverage_matrix[1,2],"-fold",sep="")
	plot(coverage_matrix[2:length(coverage_matrix[,1]),1], coverage_matrix[2:length(coverage_matrix[,1]),2], type="l", col=colors[1], xlab="Number of reads", ylab=paste("Number of covered ", pattern, "s", sep=""), main="Coverage analysis", sub=paste("Total number of ", pattern, "s: ", maxPos, sep=""))
	for(i in 1:((length(coverage_matrix[1,]))-2)){
		lines(coverage_matrix[2:length(coverage_matrix[,1]),1], coverage_matrix[2:length(coverage_matrix[,1]),i+2], col=colors[i+1])
		legend=c(legend, paste("min. ", coverage_matrix[1,i+2],"-fold",sep=""))		
	}
	legend(1, coverage_matrix[length(coverage_matrix[,1]), 2], legend, fill=colors)
	abline(maxPos, 0)
	
}
