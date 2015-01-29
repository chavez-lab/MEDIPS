##########################################################################
##Function to export a genome vector as wiggle file for a suitable genome browser (e.g. UCSC).
##########################################################################
##Input:	MEDIPS SET and/or COUPLING SET
##Param:	Set, CSet, (output)file, format, descr
##Output:	Wiggle file
##Modified:	11/01/2015
##Author:	Lukas Chavez

MEDIPS.exportWIG <-
function(Set=NULL, CSet=NULL, file=NULL, format="rpkm", descr=""){
	
	##Proof for correctness
	########################
	if(is.null(file)){stop("Must specify an output file!")}
	if(format != "count" & format!= "rpkm" & format != "pdensity" & format != "rms"){stop("Format pareeter must be either count, rpkm or pdensity!")}
	
	##Set sample names and output data
	###################################
	if(format=="pdensity"){		
		if(class(CSet)!="COUPLINGset"){stop("For exporting pattern densities, you have to specify a COUPLINGset object at CSet!")}
		window_size=window_size(CSet)
		chr_names=chr_names(CSet)
		chr_lengths=chr_lengths(CSet)
		chromosomes=chr_names(CSet)
		output_data = genome_CF(CSet)
		header_out=paste("track type=wiggle_0 name=\"", descr, "\" description=\"", descr, "\" visibility=full autoScale=on color=0,200,100 maxHeightPixels=100:50:20 graphType=bar priority=20", sep="")			
	}
	else if(format=="count"){	
		if(class(Set)!="MEDIPSset"){stop("For exporting count data, you have to specify a MEDIPSset object at Set!")}
		window_size=window_size(Set)
		chr_names=chr_names(Set)
		chr_lengths=chr_lengths(Set)
		chromosomes=chr_names(Set)		
		output_data = genome_count(Set)
		header_out=paste("track type=wiggle_0 name=\"", descr, "\" description=\"", descr, "\" visibility=full autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20", sep="")	
	}
	else if(format=="rpkm"){
		if(class(Set)!="MEDIPSset"){stop("For exporting rpkm data, you have to specify a MEDIPSset object at Set!")}
		window_size=window_size(Set)
		chr_names=chr_names(Set)
		chr_lengths=chr_lengths(Set)
		chromosomes=chr_names(Set)				
		output_data = round((genome_count(Set)*10^9)/(window_size*number_regions(Set)), digits=2)
		header_out=paste("track type=wiggle_0 name=\"", descr, "\" description=\"", descr, "\" visibility=full autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20", sep="")	
	}
	else if(format=="rms"){
		if(class(Set)!="MEDIPSset"){stop("For exporting rms data, you have to specify a MEDIPSset object at Set!")}
                if(class(CSet)!="COUPLINGset"){stop("For exporting rms data, you have to specify a COUPLINGset at CSet!")}
		window_size=window_size(Set)
                chr_names=chr_names(Set)
                chr_lengths=chr_lengths(Set)
                chromosomes=chr_names(Set)
		ccObj = MEDIPS.calibrationCurve(MSet=Set, CSet=CSet, input=F)
		output_data = MEDIPS.rms(Set, CSet, ccObj=ccObj)
                header_out=paste("track type=wiggle_0 name=\"", descr, "\" description=\"", descr, "\" visibility=full autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20", sep="")
	}
	
	##Calculate genomic coordinates
	##################################
	no_chr_windows = ceiling(chr_lengths/window_size)
	supersize_chr = cumsum(no_chr_windows)	
	genome_chr = as.vector(seqnames(MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)))
		
	write.table(header_out, file=file, sep="", quote=F, row.names=F, col.names=F)
	
	for(i in 1:length(chr_names)){
		cat(paste("Writing data for ", chr_names[i], "...\n", sep=""))
		chr_header=paste("fixedStep chrom=", chr_names[i]," start=1 step=", window_size, " span=",  window_size, sep="")
		write.table(chr_header, file=file, sep="", quote=F, row.names=F, col.names=F, append=T)
		temp_data=output_data[genome_chr==chr_names[i]]
		temp_data=temp_data[-length(temp_data)]
		write.table(temp_data, file=file, sep="", quote=F, row.names=F, col.names=F, append=T)
	}
	
	
}
