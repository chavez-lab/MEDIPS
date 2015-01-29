##########
#Function to plot the results of the saturation analysis.
##########
##Input:	A saturation result object.
##Param:	saturationObj, main (characrter)
##Output:	Plot
##Modified:	11/10/2011
##Author:	Lukas Chavez

MEDIPS.plotSaturation <-
function(saturationObj=NULL, main="Saturation analysis"){
	
	if(is.null(saturationObj)){stop("Must specify a saturation object.")}
	distinctSets=saturationObj$distinct
	estimation=saturationObj$estimation
	maxEstCor=saturationObj$maxEstCor
	maxTruCor=saturationObj$maxTruCor
		
	#Set the colors for the output plot	
	colors=c("red", "blue")
	legend=c("Saturation", "Estimated saturation")
		
	##Plot	
	plot(ylim=c(0,1), estimation[,1], estimation[,2], type="l", col=colors[1], xlab="Number of reads", ylab="Pearson correlation coefficient", main=main, sub=paste("Saturation cor: ", round(maxTruCor[2], digits=2), "; Estimated cor: ", round(maxEstCor[2], digits=2), sep=""))
	lines(distinctSets[,1], distinctSets[,2], col=colors[2])		
	legend("bottomright", maxEstCor[2], legend, fill=c("blue", "red"))
	abline(1, 0)	
}
