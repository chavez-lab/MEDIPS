#### MEDIPS class definition ###
setClass(Class = 'MEDIPSset', 
	representation = representation(
					sample_name='character',
					path_name='character',
					genome_name='character', 
					number_regions='numeric', 
					chr_names='character',					
					chr_lengths='numeric', 
					genome_count='numeric',					
					window_size='numeric', 
					extend='numeric',
					shifted='numeric',
					uniq='numeric'				
					),
	prototype = prototype(),
	validity = function(object){
		if(FALSE){stop("")}
		return(TRUE)
	}
)

### show ###
setMethod('show', signature='MEDIPSset', definition=function(object) {
	cat("S4 Object of class MEDIPSset")
	cat("\n=======================================\n")
	cat(paste("Regions file: ", sample_name(object), "\n"))
	cat(paste("File path: ", path_name(object), "\n"))
	cat(paste("Genome: ", genome_name(object), '\n'))
	cat(paste("Number of regions: ", number_regions(object), '\n'))
	cat("Chromosomes: ")
	cat(chr_names(object))
	cat("\nChromosome lengths: ")
	cat(chr_lengths(object))	
	cat(paste("\nGenome wide window size: ", window_size(object), '\n'))
	cat(paste("Reads extended to: ", extend(object), '\n'))
	cat(paste("Reads shifted by: ", shifted(object), '\n'))
	cat(paste("Parameter uniq: ", uniq(object), '\n'))	
	}
)

#### MEDIPSroi class definition ###
setClass(Class = 'MEDIPSroiSet', 
	representation = representation(
					sample_name='character',
					path_name='character',
					genome_name='character', 
					number_regions='numeric', 
					chr_names='character',					
					chr_lengths='numeric', 
					genome_count='numeric',					
					bin_number='numeric', 
					extend='numeric',
					shifted='numeric',
					uniq='numeric',
					ROI='GRanges'	
					),
	prototype = prototype(),
	validity = function(object){
		if(FALSE){stop("")}
		return(TRUE)
	}
)

### show ###
setMethod('show', signature='MEDIPSroiSet', definition=function(object) {
	cat("S4 Object of class MEDIPSset")
	cat("\n=======================================\n")
	cat(paste("Regions file: ", sample_name(object), "\n"))
	cat(paste("File path: ", path_name(object), "\n"))
	cat(paste("Genome: ", genome_name(object), '\n'))
	cat(paste("Number of regions: ", number_regions(object), '\n'))
	cat("Chromosomes: ")
	cat(chr_names(object))
	cat("\nChromosome lengths: ")
	cat(chr_lengths(object))	
	cat(paste("\nNumber of bins per ROI: ", bin_number(object), '\n'))
	cat(paste("Reads extended to: ", extend(object), '\n'))
	cat(paste("Reads shifted by: ", shifted(object), '\n'))
	cat(paste("Parameter uniq: ", uniq(object), '\n'))
	cat("ROIs: \n")
	print(rois(object))
	}
)

#### COUPLING class definition ###
setClass(Class = 'COUPLINGset', 
	representation = representation(
					seq_pattern='character',
					genome_name='character', 
					genome_CF='numeric',
					number_pattern='numeric',
					window_size='numeric', 
					chr_names='character',
					chr_lengths='numeric'
					),
	prototype = prototype(),
	validity = function(object){
		if(FALSE){stop("")}
		return(TRUE)
	}
)

### show ###
setMethod('show', signature='COUPLINGset', definition=function(object) {
	cat("S4 Object of class COUPLINGset")
	cat("\n=======================================\n")
	cat(paste("Sequence pattern: ", seq_pattern(object), "\n"))
	cat(paste("Genome: ", genome_name(object), "\n"))
	cat("Chromosomes: ")
	cat(chr_names(object))
	cat("\nChromosome lengths: ")
	cat(chr_lengths(object))	
	cat(paste("\nNumber of sequence pattern in genome: ", number_pattern(object), '\n'))
	cat(paste("Window size: ", window_size(object), '\n'))	
	}
)

#### Methods for extracting slots ###
# Genome
setGeneric('genome_name', function(object) standardGeneric('genome_name'))
setMethod('genome_name', 'MEDIPSset', function(object) object@genome_name)
setMethod('genome_name', 'MEDIPSroiSet', function(object) object@genome_name)
setMethod('genome_name', 'COUPLINGset', function(object) object@genome_name)

# window size
setGeneric('window_size', function(object) standardGeneric('window_size'))
setMethod('window_size','MEDIPSset', function(object) object@window_size)
setMethod('window_size','COUPLINGset', function(object) object@window_size)

# bin number
setGeneric('bin_number', function(object) standardGeneric('bin_number'))
setMethod('bin_number','MEDIPSroiSet', function(object) object@bin_number)

# chr names
setGeneric('chr_names', function(object) standardGeneric('chr_names'))
setMethod('chr_names','COUPLINGset', function(object) object@chr_names)
setMethod('chr_names','MEDIPSset', function(object) object@chr_names)
setMethod('chr_names','MEDIPSroiSet', function(object) object@chr_names)

# chr lengths
setGeneric('chr_lengths', function(object) standardGeneric('chr_lengths'))
setMethod('chr_lengths','MEDIPSset', function(object) object@chr_lengths)
setMethod('chr_lengths','MEDIPSroiSet', function(object) object@chr_lengths)
setMethod('chr_lengths','COUPLINGset', function(object) object@chr_lengths)

# pattern
setGeneric('seq_pattern', function(object) standardGeneric('seq_pattern'))
setMethod('seq_pattern','COUPLINGset', function(object) object@seq_pattern)

# coupling factor per window
setGeneric('genome_CF', function(object) standardGeneric('genome_CF'))
setMethod('genome_CF','COUPLINGset', function(object) object@genome_CF)

# Coupling factors
setGeneric('genome_CF', function(object) standardGeneric('genome_CF'))
setMethod('genome_CF','COUPLINGset', function(object) object@genome_CF)

# number pattern 
setGeneric('number_pattern', function(object) standardGeneric('number_pattern'))
setMethod('number_pattern','COUPLINGset', function(object) object@number_pattern)

# sample name
setGeneric('sample_name', function(object) standardGeneric('sample_name'))
setMethod('sample_name','MEDIPSset', function(object) object@sample_name)
setMethod('sample_name','MEDIPSroiSet', function(object) object@sample_name)

# path name
setGeneric('path_name', function(object) standardGeneric('path_name'))
setMethod('path_name','MEDIPSset', function(object) object@path_name)
setMethod('path_name','MEDIPSroiSet', function(object) object@path_name)

# number of regions
setGeneric('number_regions', function(object) standardGeneric('number_regions'))
setMethod('number_regions','MEDIPSset', function(object) object@number_regions)
setMethod('number_regions','MEDIPSroiSet', function(object) object@number_regions)

# genome count data per window
setGeneric('genome_count', function(object) standardGeneric('genome_count'))
setMethod('genome_count','MEDIPSset', function(object) object@genome_count)
setMethod('genome_count','MEDIPSroiSet', function(object) object@genome_count)

# extend
setGeneric('extend', function(object) standardGeneric('extend'))
setMethod('extend','MEDIPSset', function(object) object@extend)
setMethod('extend','MEDIPSroiSet', function(object) object@extend)

# shifted
setGeneric('shifted', function(object) standardGeneric('shifted'))
setMethod('shifted','MEDIPSset', function(object) object@shifted)
setMethod('shifted','MEDIPSroiSet', function(object) object@shifted)

# uniq
setGeneric('uniq', function(object) standardGeneric('uniq'))
setMethod('uniq','MEDIPSset', function(object) object@uniq)
setMethod('uniq','MEDIPSroiSet', function(object) object@uniq)

# rois
setGeneric('rois', function(object) standardGeneric('rois'))
setMethod('rois','MEDIPSroiSet', function(object) object@ROI)

