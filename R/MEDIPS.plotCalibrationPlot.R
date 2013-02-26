##########################################################################
##Plots the calibration plot
##########################################################################
##Input:	MEDIPS SET and COUPLING SET
##Param:	MSet, CSet, plot_chr, rpkm, main, xrange, input
##Output:	Calibration Plot
##Modified:	11/15/2011
##Author:	Lukas Chavez

MEDIPS.plotCalibrationPlot <- function(MSet=NULL, CSet=NULL, plot_chr="chr1", rpkm=F, main="Calibration plot", xrange=T, input=FALSE){
	
	##Proof data accordance....
	##########################
	if(class(MSet)!="MEDIPSset") stop("Must specify a MEDIPSset object!")	
	if(class(CSet)!="COUPLINGset") stop("Must specify a COUPLINGset object!")
	if(window_size(MSet)!=window_size(CSet)) stop("MEDIPSset and COUPLINGset have different window sizes!")
	if(length(chr_names(MSet))!=length(chr_names(CSet))) stop("MEDIPSset and COUPLINGset contain a different number of chromosomes!")
	for(i in 1:length(chr_names(MSet))){
		if(chr_names(MSet)[i]!=chr_names(CSet)[i]){stop("MEDIPSset and COUPLINGset contain different chomosomes!")}
	}		
	
	signal=genome_count(MSet)
	coupling=genome_CF(CSet)
	seq_pattern=seq_pattern(CSet)
	chr_lengths = chr_lengths(MSet)
	window_size = window_size(MSet)
	number_regions = number_regions(MSet)
	chromosomes=chr_names(MSet)
	
	##Calculate calibration curve
	#####################################
	ccObj = MEDIPS.calibrationCurve(MSet=MSet, CSet=CSet, input)
	mean_signal = ccObj$mean_signal
	coupling_level = ccObj$coupling_level
	intercept = ccObj$intercept
	slope = ccObj$slope
	rm(ccObj)
	gc()	
		
	##Check, if a subset of chromosomes has been selected
	######################################################
	if(plot_chr!="all"){
		cat(paste("Extracting data for",plot_chr, "...\n", sep=" "))
		
		##Calculate genomic coordinates
		##################################
		no_chr_windows = ceiling(chr_lengths/window_size)
		supersize_chr = cumsum(no_chr_windows)	
		genome_chr = as.vector(seqnames(MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)))
				
		if(length(genome_chr[genome_chr==plot_chr])==0){stop("Stated calibration chromosome does not exist within the MEDIPS SET.")}
		signal = signal[genome_chr==plot_chr]
		coupling = coupling[genome_chr==plot_chr]
		cat(paste("Plotting calibration plot for", plot_chr, "...\n", sep=" "))
	}
	if(plot_chr=="all"){cat("Plotting calibration plot for all chromosomes. [It is recommended to redirect the output to a graphic device.]\n")}
	
	
	##Calculate linear model
	#########################
	if(!input){
		y_values=seq(0,max(coupling),1)
		x_values_weighted=intercept+slope*y_values
	}
		
	##Preparations for raw data plots
	#################################
	if(!rpkm){		
		descSignal="#reads/window"
		if(xrange){						
			coupling=coupling[signal<=(max(mean_signal)*2)]
			signal=signal[signal<=(max(mean_signal)*2)]			
		}			
	}
	##Preparations for rpkm data plots
	#################################
	else{
		descSignal="rpkm"
		
		if(!input){x_values_weighted = (x_values_weighted*10^9)/(window_size*number_regions(MSet))}
		
		signal = (signal*10^9)/(window_size*number_regions(MSet))
		mean_signal = (mean_signal*10^9)/(window_size*number_regions(MSet))
		
		if(xrange){			
			coupling = coupling[signal<=(max(mean_signal)*2)]
			signal = signal[signal<=(max(mean_signal)*2)]					
		}		
	}			
		
	##Plot
	#######
	plot(signal, coupling, pch=".", main=main, xlab=paste(descSignal, "", sep=""), ylab=paste(seq_pattern, " coupling factor", sep=""), col="lightblue")		
	lines(mean_signal,  coupling_level, col="red")	
	if(!input){lines(x_values_weighted, y_values, col="green")}
	if(!input){legend(1, max(coupling), c("Genomic window", "Mean for coupling level", "Estimated linear fit"), fill=c("lightblue", "red", "green"))}
	else{legend(1, max(coupling), c("Genomic window", "Mean for coupling level"), fill=c("lightblue", "red"))}
	
}
