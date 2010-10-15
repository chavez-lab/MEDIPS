###############
##Function takes a MEDIPSset object, calculates the calibration curve and the parameters intercept and slope.
###############
MEDIPS.calibrationCurve <- function(data=NULL){
	
	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")
	
	cal_chr="all"
	signal=genome_raw(data)
	coupling=genome_CF(data)	
			
	cat("Calculate calibration curve...\n")	
	maxCoup=floor(max(coupling))
			
	mean_signal=vector(mode="numeric",length=maxCoup+1)
	mean_coupling=vector(mode="numeric",length=maxCoup+1)

	out=.C("calibration",
         as.integer(maxCoup),
	 as.integer(coupling),
	 length(coupling),
	 as.integer(signal),
	 length(signal),
	 a=as.numeric(mean_signal),
	 b=as.integer(mean_coupling),
	 PACKAGE="MEDIPS"
	 )
	
 	mean_signal=out$a
 	mean_coupling=out$b
	calCurve_signal=mean_signal[mean_signal!=maxCoup+5 | mean_coupling!=maxCoup+5]
	calCurve_CC=mean_coupling[mean_signal!=maxCoup+5 | mean_coupling!=maxCoup+5] 	
 	
	##Calculate max.CC index
	highest_signal=0
	count_decrease_steps=0
	max_CCurve_signal_index=NULL
	former=9999999	
	
	for(i in 1:length(calCurve_signal)){
		##Test, if calCurve_signal decreases
		if(i!=1 & calCurve_signal[i]<former){
			count_decrease_steps=count_decrease_steps+1
			if(count_decrease_steps==3){
				max_CCurve_signal_index=i-3
				break
			}	
		}		
		if(calCurve_signal[i]>=highest_signal){
			highest_signal=calCurve_signal[i]
			count_decrease_steps=0
		}		
		former=calCurve_signal[i]
	}
	
	##Linear curve by linear least squares regression
	##################
	cat("Performing linear regression...\n")
	curve_signal_lowCpG=NULL
	curve_CF_lowCpG=NULL
	for(i in 1:length(calCurve_signal)){	
		if(i<=max_CCurve_signal_index){		
			curve_signal_lowCpG=c(curve_signal_lowCpG, calCurve_signal[i])
			curve_CF_lowCpG=c(curve_CF_lowCpG, calCurve_CC[i])
		}			
	}

	fit=lm(curve_CF_lowCpG ~ curve_signal_lowCpG)
	intercept=fit$coefficients[[1]]
	slope=fit$coefficients[[2]]
		
	if(slope==0){stop("Estimated slope is 0!?")}
	
	MEDIPSsetObj = new('MEDIPSset', cali_chr=as.character(cal_chr), calcurve_mean_signals=mean_signal, calcurve_mean_coupling=mean_coupling, intercept=intercept, slope=slope, genome_CF=genome_CF(data), fragmentLength=fragmentLength(data), distFunction=distFunction(data), distFile=distFile(data), seq_pattern=seq_pattern(data), pattern_chr=pattern_chr(data), pattern_pos=pattern_pos(data), number_pattern=number_pattern(data), genome_chr=genome_chr(data), genome_pos=genome_pos(data), genome_raw=genome_raw(data), extend=extend(data), bin_size=bin_size(data), sample_name=sample_name(data), genome_name=genome_name(data), regions_chr=regions_chr(data), regions_start=regions_start(data), regions_stop=regions_stop(data), regions_strand=regions_strand(data), number_regions=number_regions(data), chr_names=chr_names(data), chr_lengths=chr_lengths(data))
	gc()
	return(MEDIPSsetObj)
}
