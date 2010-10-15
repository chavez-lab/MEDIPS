##########
##Function to plot the calibration plot
##########
MEDIPS.plotCalibrationPlot <- function(data=NULL, xrange=NULL, linearFit=FALSE, plot_chr="all", rpm=F){
	
	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")	
	
	signal=genome_raw(data)
	coupling=genome_CF(data)
	seq_pattern=seq_pattern(data)	
	
	##Get calibration curve informations
	#####################################
	calcurve_mean_signals=calcurve_mean_signals(data)
	calcurve_mean_coupling=calcurve_mean_coupling(data)
	
	##Get Linear regression informations (necessary, if linearFit is selected)
	#####################################
	intercept=intercept(data)
	slope=slope(data)
		
	##Check, if a subset of chromosomes has been selected
	######################################################
	if(plot_chr!="all"){
		cat(paste("Extract data for",plot_chr, "...\n", sep=" "))
		genome_chr=genome_chr(data)
		if(length(genome_chr[genome_chr==plot_chr])==0){stop("Stated calibration chromosome does not exist within the genome vector.")}
		signal=signal[genome_chr==plot_chr]
		coupling=coupling[genome_chr==plot_chr]
		cat(paste("Plotting calibration plot for", plot_chr, "...\n", sep=" "))
	}
	if(plot_chr=="all"){cat("Plotting calibration plot for all chromosomes. It is recommended to call a png() function before!\n")}
	
	
	##Check, if linearFit has been selected
	#######################################
	if(linearFit){
		x_values=seq(0,max(signal),1)
		x_values_weighted=intercept+x_values*slope
	}
	
	##Preparations for raw data plots
	#################################
	if(!rpm){		
		descSignal="#reads/bin"
		if(!is.null(xrange)){			
			coupling=coupling[signal<=xrange]
			signal=signal[signal<=xrange]			
		}			
	}
	##Preparations for rpm data plots
	#################################
	else{
		descSignal="rpm"
		if(linearFit){
			x_values=x_values/(number_regions(data)/1000000)			
		}
		signal=signal/(number_regions(data)/1000000)
		calcurve_mean_signals=calcurve_mean_signals(data)/(number_regions(data)/1000000)						
		if(!is.null(xrange)){			
			coupling=coupling[signal<=xrange]						
			signal=signal[signal<=xrange]			
		}		
	}	
			
	
	##Plot
	#######
	plot(signal, coupling, pch=".", main="Calibration plot", xlab=paste(descSignal, "", sep=""), ylab=paste(seq_pattern, " coupling factor", sep=""), col="lightblue")		
	lines(calcurve_mean_signals,  calcurve_mean_coupling, col="red")	
	if(linearFit){
		lines(x_values, x_values_weighted, col="green")		
		legend(1, max(coupling), c("Genomic window", "Mean for coupling level", "Estimated linear fit"), fill=c("lightblue", "red", "green"))	
	}
	else{
		legend(1, max(coupling), c("Genomic window", "Mean for coupling level"), fill=c("lightblue", "red"))	
	}		
}
