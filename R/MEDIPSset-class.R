#### MEDIPS class definition ###
setClass(Class = 'MEDIPSset', 
	representation = representation(
					regions_chr='character', 
					regions_start='numeric', 
					regions_stop='numeric',
					regions_strand='character', 
					number_regions='numeric', 
					pattern_chr='character', 
					pattern_pos='numeric', 
					number_pattern='numeric', 
					genome_chr='character', 
					genome_pos='numeric', 
					genome_CF='numeric', 
					genome_raw='numeric', 
					genome_norm='numeric', 
					genome_name='character', 
					bin_size='numeric', 
					extend='numeric', 
					fragmentLength='numeric', 
					sample_name='character', 
					chr_lengths='numeric', 
					chr_names='character', 
					seq_pattern='character', 
					distFunction='character', 
					distFile='character', 
					calcurve_mean_signals='numeric', 
					calcurve_mean_coupling='numeric', 
					calcurve_var='numeric', 
					intercept='numeric', 
					slope='numeric', 
					cali_chr='character'
					),
	prototype = prototype(),
	validity = function(object){
		if(FALSE){stop("")}
		return(TRUE)
	}
)

#### Methods for extracting slots ###
# regions chr
setGeneric('regions_chr', function(object) standardGeneric('regions_chr'))
setMethod('regions_chr','MEDIPSset', function(object) object@regions_chr)
# regions start
setGeneric('regions_start', function(object) standardGeneric('regions_start'))
setMethod('regions_start','MEDIPSset', function(object) object@regions_start)
# regions stop
setGeneric('regions_stop', function(object) standardGeneric('regions_stop'))
setMethod('regions_stop','MEDIPSset', function(object) object@regions_stop)
# regions strand
setGeneric('regions_strand', function(object) standardGeneric('regions_strand'))
setMethod('regions_strand','MEDIPSset', function(object) object@regions_strand)
# number of regions
setGeneric('number_regions', function(object) standardGeneric('number_regions'))
setMethod('number_regions','MEDIPSset', function(object) object@number_regions)
# genome chr
setGeneric('genome_chr', function(object) standardGeneric('genome_chr'))
setMethod('genome_chr','MEDIPSset', function(object) object@genome_chr)
# pattern
setGeneric('seq_pattern', function(object) standardGeneric('seq_pattern'))
setMethod('seq_pattern','MEDIPSset', function(object) object@seq_pattern)
# pattern chr
setGeneric('pattern_chr', function(object) standardGeneric('pattern_chr'))
setMethod('pattern_chr','MEDIPSset', function(object) object@pattern_chr)
# pattern pos
setGeneric('pattern_pos', function(object) standardGeneric('pattern_pos'))
setMethod('pattern_pos','MEDIPSset', function(object) object@pattern_pos)
# number pattern 
setGeneric('number_pattern', function(object) standardGeneric('number_pattern'))
setMethod('number_pattern','MEDIPSset', function(object) object@number_pattern)
# genome position
setGeneric('genome_pos', function(object) standardGeneric('genome_pos'))
setMethod('genome_pos','MEDIPSset', function(object) object@genome_pos)
# genome coupling factor
setGeneric('genome_CF', function(object) standardGeneric('genome_CF'))
setMethod('genome_CF','MEDIPSset', function(object) object@genome_CF)
# genome raw data per bin
setGeneric('genome_raw', function(object) standardGeneric('genome_raw'))
setMethod('genome_raw','MEDIPSset', function(object) object@genome_raw)
# genome normalized data per bin
setGeneric('genome_norm', function(object) standardGeneric('genome_norm'))
setMethod('genome_norm','MEDIPSset', function(object) object@genome_norm)
# organism
setGeneric('genome_name', function(object) standardGeneric('genome_name'))
setMethod('genome_name','MEDIPSset', function(object) object@genome_name)
# bin size
setGeneric('bin_size', function(object) standardGeneric('bin_size'))
setMethod('bin_size','MEDIPSset', function(object) object@bin_size)
# extend
setGeneric('extend', function(object) standardGeneric('extend'))
setMethod('extend','MEDIPSset', function(object) object@extend)
# fragment length
setGeneric('fragmentLength', function(object) standardGeneric('fragmentLength'))
setMethod('fragmentLength','MEDIPSset', function(object) object@fragmentLength)
# sample name
setGeneric('sample_name', function(object) standardGeneric('sample_name'))
setMethod('sample_name','MEDIPSset', function(object) object@sample_name)
# chr lengths
setGeneric('chr_lengths', function(object) standardGeneric('chr_lengths'))
setMethod('chr_lengths','MEDIPSset', function(object) object@chr_lengths)
# chr names
setGeneric('chr_names', function(object) standardGeneric('chr_names'))
setMethod('chr_names','MEDIPSset', function(object) object@chr_names)
# distance function
setGeneric('distFunction', function(object) standardGeneric('distFunction'))
setMethod('distFunction','MEDIPSset', function(object) object@distFunction)
# distance file
setGeneric('distFile', function(object) standardGeneric('distFile'))
setMethod('distFile','MEDIPSset', function(object) object@distFile)
# calibration curve mean signals
setGeneric('calcurve_mean_signals', function(object) standardGeneric('calcurve_mean_signals'))
setMethod('calcurve_mean_signals','MEDIPSset', function(object) object@calcurve_mean_signals)
# calibration curve mean coupling
setGeneric('calcurve_mean_coupling', function(object) standardGeneric('calcurve_mean_coupling'))
setMethod('calcurve_mean_coupling','MEDIPSset', function(object) object@calcurve_mean_coupling)
# calibration curve variance
setGeneric('calcurve_var', function(object) standardGeneric('calcurve_var'))
setMethod('calcurve_var','MEDIPSset', function(object) object@calcurve_var)
# intercept
setGeneric('intercept', function(object) standardGeneric('intercept'))
setMethod('intercept','MEDIPSset', function(object) object@intercept)
# slope
setGeneric('slope', function(object) standardGeneric('slope'))
setMethod('slope','MEDIPSset', function(object) object@slope)
# cali_chr
setGeneric('cali_chr', function(object) standardGeneric('cali_chr'))
setMethod('cali_chr','MEDIPSset', function(object) object@cali_chr)

### show ###
setMethod('show', signature='MEDIPSset', definition=function(object) {
	cat("S4 Object of class MEDIPSset")
	cat("\n=======================================\n")
	cat("Regions information\n")
	cat("=======================================\n")
	cat(paste("Regions file: ", object@sample_name, "\n"))
	cat(paste("Organism: ", object@genome_name, '\n'))
	cat("Chromosomes: ")
	cat(object@chr_names)
	cat("\nChromosome lengths: ")
	cat(object@chr_lengths)
	cat(paste("\nNumber of regions: ", object@number_regions, '\n'))
	max_out_regions = min(3, length(object@regions_chr))
	cat("Regions chromosomes: ")
	cat(object@regions_chr[1:max_out_regions])
	cat("...\nRegions start positions: ")
	cat(object@regions_start[1:max_out_regions])
	cat("...\nRegions stop positions: ")
	cat(object@regions_stop[1:max_out_regions])
	cat("...\nRegions strand: ")
	cat(object@regions_strand[1:max_out_regions])
	cat("...\n=======================================\n")
	cat("Genome vector signals information\n")
	cat("=======================================\n")
	cat(paste("Genome wide bin size: ", object@bin_size, '\n'))
	cat(paste("Reads extended by: ", object@extend, '\n'))
	max_out_genome = min(3, length(object@regions_chr))
	cat("Genome vector chromosomes: ")
	cat(object@genome_chr[1:max_out_genome])
	cat("...\nGenome vector positions: ")
	cat(object@genome_pos[1:max_out_genome])
	cat("...\nGenome vector signals: ")
	cat(object@genome_raw[1:max_out_genome])
	cat("...\n=======================================\n")
	cat("Pattern information\n")
	cat("=======================================\n")
	cat(paste("Pattern: ", object@seq_pattern, '\n'))
	cat(paste("Number of patterns: ", object@number_pattern, '\n'))
	max_out_pattern = min(3, object@number_pattern)
	cat("Pattern chromosomes: ")
	cat(object@pattern_chr[1:max_out_pattern])
	cat("...\nPattern positions: ")
	cat(object@pattern_pos[1:max_out_pattern])
	cat("...\n=======================================\n")
	cat("Genome vector coupling factor information\n")
	cat("=======================================\n")	
	cat(paste("Distance function: ", object@distFunction, '\n'))
	cat(paste("Distance file: ", object@distFile, '\n'))
	cat(paste("Fragment length: ", object@fragmentLength, '\n'))
	cat("Genome vector coupling factors: ")
	max_out_CF = min(3, length(object@genome_CF))
	cat(object@genome_CF[1:max_out_CF])
	cat("...\n=======================================\n")
	cat("Calibration information\n")
	cat("=======================================\n")
	max_out_ccurve = min(3, length(object@calcurve_mean_signals))
	cat("Calibration curve mean signals: ")
	cat(object@calcurve_mean_signals[1:max_out_ccurve])
	cat("...\nCalibration curve mean coupling factors: ")
	cat(object@calcurve_mean_coupling[1:max_out_ccurve])
	cat("...\nCalibration curve variance: ")
	cat(object@calcurve_var[1:max_out_ccurve])
	cat(paste("...\nIntercept: ", object@intercept))
	cat(paste("\nSlope: ", object@slope))
	cat(paste("\nCalibration chromosome: ", object@cali_chr))
	cat("\n=======================================\n")
	cat("Genome vector normalized signal information\n")
	cat("=======================================\n")
	cat("Normalization output interval: [0:1000]\n")
	cat("Genome vector normalized signals: ")
	max_out_norm = min(3, length(object@genome_norm))
	cat(object@genome_norm[1:max_out_norm])
	cat("...\n")	
	}
)

