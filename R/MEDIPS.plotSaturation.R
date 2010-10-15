##########
#Function to plot the results of the saturation analysis
##########
MEDIPS.plotSaturation <-
function(saturationObj=NULL){
	
	if(is.null(saturationObj)){stop("Must specify a saturation object.")}
	distinctSets=saturationObj$distinct
	estimation=saturationObj$estimation
	maxEstCor=saturationObj$maxEstCor
	maxTruCor=saturationObj$maxTruCor
		
	#Set the colors for the output plot	
	colors=c("red", "blue")
	#legend=c(paste("Saturation for ", maxTruCor[1], " reads", sep=""), paste("Estimated saturation for ", maxEstCor[1], " reads", sep=""))
	legend=c("Saturation", "Estimated saturation")
		
	##Plot	
	plot(estimation[,1], estimation[,2], type="l", col=colors[1], xlab="Number of reads", ylab="Pearson correlation coefficient", main="Saturation analysis", sub=paste("Saturation cor: ", round(maxTruCor[2], digits=2), "; Estimated cor: ", round(maxEstCor[2], digits=2), sep=""))
	lines(distinctSets[,1], distinctSets[,2], col=colors[2])		
	legend(1, maxEstCor[2], legend, fill=c("blue", "red"))
	abline(1, 0)	
}
