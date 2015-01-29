##########################################################################
##Plots the calibration plot
##########################################################################
##Input:	MEDIPS SET and COUPLING SET
##Param:	MSet, CSet, plot_chr, rpkm, main, xrange, input
##Output:	Calibration Plot
##Modified:	03/18/2013
##Author:	Lukas Chavez, Matthias Lienhard

MEDIPS.plotCalibrationPlot <- function(MSet=NULL, ISet=NULL, CSet=NULL, plot_chr="all", rpkm=T, main="Calibration plot", xrange=T){
	
	##Proof data accordance....
	##########################
	input=F

	if(!(is.null(MSet)) & (class(MSet)!="MEDIPSset" & class(MSet)!="MEDIPSroiSet")) stop("MSet must be a MEDIPSset or MEDIPSroiSet object!")	
	if(!(is.null(ISet)) & (class(ISet)!="MEDIPSset" & class(ISet)!="MEDIPSroiSet")) stop("ISet must be a MEDIPSset or MEDIPSroiSet object!")	
	if(is.null(MSet) & is.null(ISet) ) stop("ISet and/or MSet must be specified")	

	#print("Checking CSet")
	if(class(CSet)!="COUPLINGset") stop("Must specify a COUPLINGset object!")
	
	#print("Checking MSet")
	if(!is.null(MSet) & class(MSet)=="MEDIPSset"){
	  if(window_size(MSet)!=window_size(CSet)) stop("MSet and COUPLINGset have different window sizes!")
	  if(length(chr_names(MSet))!=length(chr_names(CSet))) stop("MSet and COUPLINGset contain a different number of chromosomes!")
	  for(i in 1:length(chr_names(MSet))){
		if(chr_names(MSet)[i]!=chr_names(CSet)[i]){stop("MSset and COUPLINGset contain different chomosomes!")}
	  }
	}
	
	if(!is.null(MSet) & class(MSet)=="MEDIPSroiSet"){
		if( length(genome_count(MSet)) != length(genome_CF(CSet)) ){stop("MSet and COUPLINGset have different number of regions of interest!")}
		if( length(chr_names(MSet)) != length(chr_names(CSet)) ){stop("MSet and COUPLINGset contain a different number of chromosomes!")}
          	for(i in 1:length(chr_names(MSet))){
                	if(chr_names(MSet)[i]!=chr_names(CSet)[i]){stop("MSset and COUPLINGset contain different chomosomes!")}
          	}
	}

	#print("Checking ISet")
	if(!is.null(ISet) & class(ISet)=="MEDIPSset"){
	  if(window_size(ISet)!=window_size(CSet)) stop("ISet and COUPLINGset have different window sizes!")
	  if(length(chr_names(ISet))!=length(chr_names(CSet))) stop("ISet and COUPLINGset contain a different number of chromosomes!")
	  for(i in 1:length(chr_names(ISet))){
		if(chr_names(ISet)[i]!=chr_names(CSet)[i]){stop("ISet and COUPLINGset contain different chomosomes!")}
	  }
	}

	if(!is.null(ISet) & class(ISet)=="MEDIPSroiSet"){
		if( length(genome_count(ISet)) != length(genome_CF(CSet)) ){stop("MSet and COUPLINGset have different number of regions of interest!")}
                if( length(chr_names(ISet)) != length(chr_names(CSet)) ){stop("MSet and COUPLINGset contain a different number of chromosomes!")}
                for(i in 1:length(chr_names(ISet))){
                        if(chr_names(ISet)[i]!=chr_names(CSet)[i]){stop("MSset and COUPLINGset contain different chomosomes!")}
               	}
	}

	if(!is.null(MSet)){
		signal=genome_count(MSet)
		chr_lengths = chr_lengths(MSet)
		if(class(MSet)=="MEDIPSset"){window_size = window_size(MSet)}
		if(class(MSet)=="MEDIPSroiSet"){window_size = width(rois(MSet))}
		number_regions = number_regions(MSet)
		chromosomes=chr_names(MSet)
	}else{
		signal=genome_count(ISet)
		chr_lengths = chr_lengths(ISet)
		if(class(ISet)=="MEDIPSset"){window_size = window_size(ISet)}
                if(class(ISet)=="MEDIPSroiSet"){window_size = width(rois(ISet))}
		number_regions = number_regions(ISet)
		chromosomes=chr_names(ISet)
	}
	coupling=genome_CF(CSet)
	seq_pattern=seq_pattern(CSet)
		
	##Calculate calibration curve
	#####################################
	if (!is.null(MSet))
	 	ccObj_MSet = MEDIPS.calibrationCurve(MSet=MSet, CSet=CSet, input=F)
	if (!is.null(ISet))
		ccObj_ISet = MEDIPS.calibrationCurve(MSet=ISet, CSet=CSet, input=T)

		
	##Check, if a subset of chromosomes has been selected
	######################################################
	if(plot_chr!="all" & (class(MSet)=="MEDIPSset" | class(ISet)=="MEDIPSset")){
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

	if(plot_chr!="all" & (class(MSet)=="MEDIPSroiSet" | class(ISet)=="MEDIPSroiSet")){stop("Selecting chromosomes nor supported for regions of interest.")}

	if(plot_chr=="all"){cat("Plotting calibration plot for all chromosomes. It is recommended to redirect the output to a graphic device.\n")}

	if (!rpkm) {
        	descSignal = "#reads/window"
	}else{
		descSignal = "rpkm"
		if(!is.null(MSet)){
		  if(class(MSet)=="MEDIPSroiSet"){stop("Transformation to rpkm not supported for regions of interest.")}
		  f=10^9/(window_size * number_regions(MSet))
		  signal=signal*f
		  ccObj_MSet$slope=ccObj_MSet$slope*f
		  ccObj_MSet$intercept=ccObj_MSet$intercept*f
		  ccObj_MSet$mean_signal=ccObj_MSet$mean_signal*f
		}
		if(!is.null(ISet)){
		  if(class(ISet)=="MEDIPSroiSet"){stop("Transformation to rpkm not supported for regions of interest.")}
		  signal=signal*10^9/(window_size * number_regions(ISet))
         	  fI=10^9/(window_size * number_regions(ISet))
		  ccObj_ISet$mean_signal=ccObj_ISet$mean_signal/fI
		}
	}

	##Preparations for raw data plots
	#################################
	if(xrange){
	  if(!is.null(MSet))
	    range=c(0,max(ccObj_MSet$mean_signal)*5)
	  else
	    range=c(0,max(ccObj_ISet$mean_signal)*5)
	}else
	  range=c(0,max(signal))
	  
	##Plot
	#######
	plot(coupling,signal, pch=".", main=main, ylab=descSignal,ylim=range, xlab=paste(seq_pattern, " coupling factor", sep=""), col="lightblue")		
	for(i in 0:max(coupling)){
	  t=table(signal[coupling==i])
	  points(x=rep(i,length(t)), y=as.numeric(names(t)),lwd=log(t, 10), pch=4, col="lightblue")
	}
	if(!is.null(MSet))
	  llab="MeDIP reads in genomic window"
	else
	  llab="Input reads in genomic window"
	lcol="lightblue"
	if(! is.null(MSet)){
		lines(ccObj_MSet$coupling_level, ccObj_MSet$mean_signal, col="red")
		abline(a=ccObj_MSet$intercept, b=ccObj_MSet$slope, col="green")
		llab=c(llab,"MeDIP read density for coupling level","Estimated linear fit for MeDIP set")
		lcol=c(lcol,"red", "green")
	}
	if(! is.null(ISet)){
		lines(ccObj_ISet$coupling_level, ccObj_ISet$mean_signal, col="blue")
		llab=c(llab,"Input read density for coupling level" )
		lcol=c(lcol,"blue")
	}

	legend("topright", legend=llab, fill=lcol)
	
}
