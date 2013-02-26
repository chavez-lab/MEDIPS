##########################################################################
##Function processes a MEDIPSset and weights the counts by the estimated coupling factor dependent signal of the calibration curve.
##Normalized counts are transformed into a rpkm scale.
##########################################################################
##Input:	MEDIPS SET, CS
##Output:	rms vector
##Modified:	11/16/2011
##Author:	Lukas Chavez

MEDIPS.rms <- function(MSet=NULL, CSet=NULL, ccObj=NULL){
	
	signal = genome_count(MSet)
	coupling = genome_CF(CSet)
		
	if(is.null(ccObj)){
		##Calculate calibration curve
		#####################################
		ccObj = MEDIPS.calibrationCurve(MSet=MSet, CSet=CSet, input=F)
	}

	##Weight signals by linear regression obtained parameters
	####################
	cat("Calculating relative methylation score...\n")
	
        estim=numeric(length(ccObj$mean_signal))
	##For the low range of the calibration curve (< max_index) divide counts by observed mean count
	low_range=1:ccObj$max_index
        estim[low_range]=ccObj$mean_signal[low_range]
	##For the higher range of the calibration curve (>=max_index) divide counts by estimated count
	high_range=(ccObj$max_index+1):length(estim)
#################################################was buggy: slope*signal
	estim[high_range]=ccObj$intercept + (ccObj$slope * ccObj$coupling_level[high_range]) 
#################################################
	#rms normalization:
	signal=signal/estim[coupling+1]
	##For cf==0 set signal to zero
	signal[coupling==0]=0
	##Transform weighted signals into reads per million
	if(class(MSet)=="MEDIPSset"){
		signal = (signal*10^9)/(window_size(MSet)*number_regions(MSet))
	}
	else{
		signal = (signal*10^9)/(as.numeric(width(rois(MSet)))*number_regions(MSet))
	}
	#ist reads per million nicht eher das?
#################################################
	#signal = (signal*10^9)/sum(genome_count(MSet))
#################################################
	return(signal)
}
