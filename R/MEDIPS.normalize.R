###############
##Function takes a MEDIPSset and weights all signals by the estimated mean signal of its coupling factor,
##transforms the normalized signals into a reads/million scale and finally transforms the data range into the interval of [0:1000].
###############

##Perform full normalization on a given MEDIPS.SET.
MEDIPS.normalize <- function(data=NULL){
	
	if(class(data)!="MEDIPSset") stop("Must specify a MEDIPSset object.")
	
	signal=genome_raw(data)
	coupling=genome_CF(data)		
	intercept=intercept(data)
	slope=slope(data)
			
	##Weight signals by linear regression obtained parameters
	####################
	cat("Weight raw signals by the estimated parameters...\n")
	estimated_meanSignal=(coupling-intercept)/slope	
	if(length(estimated_meanSignal[estimated_meanSignal<=0])!=0){
		signal[estimated_meanSignal>0]=signal[estimated_meanSignal>0]/estimated_meanSignal[estimated_meanSignal>0]
		minWeight=min(estimated_meanSignal[estimated_meanSignal>0])
		signal[estimated_meanSignal<=0]=signal[estimated_meanSignal<=0]/minWeight
	}
	else{signal=signal/estimated_meanSignal}		
	
	##Transform weighted signals into reads per million
	####################
	cat("Transform weighted signals into reads per million (rpm)...\n")
	signal=signal/(number_regions(data)/1000000)
	
	MEDIPSsetObj = new('MEDIPSset', genome_norm=signal, cali_chr=cali_chr(data), calcurve_mean_signals=calcurve_mean_signals(data), calcurve_mean_coupling=calcurve_mean_coupling(data), calcurve_var=calcurve_var(data), intercept=intercept(data), slope=slope(data), genome_CF=genome_CF(data), fragmentLength=fragmentLength(data), distFunction=distFunction(data), distFile=distFile(data), seq_pattern=seq_pattern(data), pattern_chr=pattern_chr(data), pattern_pos=pattern_pos(data), number_pattern=number_pattern(data), genome_chr=genome_chr(data), genome_pos=genome_pos(data), genome_raw=genome_raw(data), extend=extend(data), bin_size=bin_size(data), sample_name=sample_name(data), genome_name=genome_name(data), regions_chr=regions_chr(data), regions_start=regions_start(data), regions_stop=regions_stop(data), regions_strand=regions_strand(data), number_regions=number_regions(data), chr_names=chr_names(data), chr_lengths=chr_lengths(data))
}

##Only transform given data into log scale and shift into interval
##################################################################
MEDIPS.transform <- function(data=NULL){
	##Log2 of signals except signal=0
	####################
	data[!is.na(data)] = log2(data[!is.na(data)])
				
	##Shift signals into positive value range
	####################
	minsignal=min(data[data!=-Inf & !is.na(data)])
	data[data!=-Inf & !is.na(data)]=data[data!=-Inf & !is.na(data)]+abs(minsignal)
	
	##Transform values into the interval [0:1000]
	####################
	maxsignal=max(data[data!=Inf & !is.na(data)])
	data[!is.na(data)]=(data[!is.na(data)]/maxsignal)*1000

	##Eliminate -Inf --> 0
	#######################
	data[data==-Inf]=0

	gc()
	return(data)
}

