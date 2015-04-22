#######################################
##Read bed file
#######################################
##Input:	tab (|) separated bed file "chr | start | stop | name | score | strand | ..."
##Param:	allignment.file, path, extend, shift, uniq, dataset
##Output:	Granges object
##Requires:	GenomicRanges
##Modified:	11/10/2011
##Author:	Lukas Chavez, Joern Dietrich, Isaac Lopez Moyado 

getGRange <-
function (fileName, path = NULL, extend, shift, chr.select = NULL, 
    dataset = NULL, uniq = 1e-3, ROI = NULL) 
{
    ext = substr(fileName, nchar(fileName) - 3, nchar(fileName))
    bam = (ext == ".bam" | ext == ".BAM")
    bamindex = bam & file.exists(paste(path, "/", fileName, ".bai", 
        sep = ""))
    if (bam) {
        scanFlag = scanBamFlag(isUnmappedQuery = F)
        if (bamindex & (!is.null(chr.select) | !is.null(ROI))) {
            if (!is.null(ROI)) {
                cat("Reading bam alignment", fileName, "\n considering ROIs using bam index\n")
                if (!is.null(extend)) {
                  ROI[, 2] = ROI[, 2] - extend
                  ROI[, 3] = ROI[, 3] + extend
                }
                if (!is.null(shift)) {
                  ROI[, 2] = ROI[, 2] - shift
                  ROI[, 3] = ROI[, 3] - shift
                }
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
            else {
                cat("Reading bam alignment", fileName, "\n considering ", 
                  chr.select, " using bam index\n")
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
            scanParam = ScanBamParam(flag = scanFlag, what = c("rname", 
                "pos", "strand", "qwidth"), which = sel)
        }
        else {
            cat("Reading bam alignment", fileName, "\n")
            scanParam = ScanBamParam(flag = scanFlag, what = c("rname", 
                "pos", "strand", "qwidth"))
        }
        regions = scanBam(file = paste(path, fileName, sep = "/"), 
            param = scanParam)
        regions = do.call(rbind, lapply(regions, as.data.frame, 
            stringsAsFactors = F))
        regions = data.frame(chr = as.character(as.vector(regions$rname)), 
            start = as.numeric(as.vector(regions$pos)), stop = as.numeric(as.vector(regions$pos) + 
                as.vector(regions$qwidth) - 1), strand = as.character(as.vector(regions$strand)), 
            stringsAsFactors = F)
    }
    else {
        cat("Reading bed alignment", fileName, "\n")
        regions = read.table(paste(path, fileName, sep = "/"), 
            sep = "\t", header = FALSE, row.names = NULL, comment.char = "", 
            colClasses = c("character", "numeric", "numeric", 
                "NULL", "NULL", "character"))
        names(regions) = c("chr", "start", "stop", "strand")
    }
    if (!is.null(chr.select) & !bamindex) {
        cat("Selecting ", chr.select, "\n")
        regions = regions[regions[, 1] %in% as.vector(chr.select), 
            ]
    }
    cat("Total number of imported short reads: ", nrow(regions), 
        "\n", sep = "")
    regions = adjustReads(regions, extend, shift)
    cat("Creating GRange Object...\n")
    regions_GRange = GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, 
        end = regions$stop), strand = regions$strand)
   
    if(is.logical(uniq)){stop("Parameter 'uniq' is not logical anymore, please specify a p-value and see the MEDIPS vignette.")}
    if (uniq == 1) {
		cat("Keep only one representative of stacked reads mapping to the same genomic location.\n", sep = "")
		regions_GRange = unique(regions_GRange)
		cat("Number of remaining reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq < 1 & uniq > 0) {
		max_dup_number = qpois(1 - as.numeric(uniq), length(regions_GRange) / 
			sum(as.numeric(seqlengths(dataset)[chr.select])))
		max_dup_number = max(1, max_dup_number)
		cat("Keep at most ", max_dup_number, 
			" reads mapping to the same genomic location\n", sep = "")
		uniq_regions = unique(regions_GRange)
		dup_number = countMatches(uniq_regions, regions_GRange)
		dup_number[dup_number > max_dup_number] = max_dup_number
		regions_GRange = rep(uniq_regions, times = dup_number)
		cat("Number of remaining reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq == 0) {
		cat("Do not correct for potential PCR artefacts (keep all reads).\n", sep = "")
	} else {
		stop("Must specify a valid value for parameter uniq. Please check MEDIPS vignette.")
	}
	strand(regions_GRange) = "*"
	return(regions_GRange)
}

getPairedGRange <-
function (fileName, path = NULL, extend, shift, chr.select = NULL, 
    dataset = NULL, uniq = 1e-3, ROI = NULL, bwa = FALSE) 
{
    if (bwa) {
        cat("Warning: processing of bwa alignment files for paired-end sequencing data is still under development.\n")
    }
    ext = substr(fileName, nchar(fileName) - 3, nchar(fileName))
    bam = (ext == ".bam" | ext == ".BAM")
    bamindex = bam & file.exists(paste(path, "/", fileName, ".bai", 
        sep = ""))
    if (bam) {
        scanFlag = scanBamFlag(isPaired = T, isProperPair = TRUE, 
            hasUnmappedMate = FALSE, isUnmappedQuery = F, isFirstMateRead = T, 
            isSecondMateRead = F)
        if (bamindex & (!is.null(chr.select) | !is.null(ROI))) {
            if (!is.null(ROI)) {
                cat("Reading bam alignment", fileName, "\n considering ROIs using bam index\n")
                if (!is.null(extend)) {
                  ROI[, 2] = ROI[, 2] - extend
                  ROI[, 3] = ROI[, 3] + extend
                }
                if (!is.null(shift)) {
                  ROI[, 2] = ROI[, 2] - shift
                  ROI[, 3] = ROI[, 3] - shift
                }
                sel = GRanges(ROI[, 1], IRanges(start = ROI[, 
                  2], end = ROI[, 3]))
            }
            else {
                cat("Reading bam alignment", fileName, "\n considering ", 
                  chr.select, " using bam index\n")
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
            if (bwa) {
                scanParam = ScanBamParam(flag = scanFlag, what = c("rname", 
                  "pos", "strand", "qwidth", "isize", "mpos"), 
                  which = sel)
            }
            else {
                scanParam = ScanBamParam(flag = scanFlag, what = c("rname", 
                  "pos", "strand", "qwidth", "isize"), which = sel)
            }
        }
        else {
            cat("Reading bam alignment", fileName, "\n")
            if (bwa) {
                scanParam = ScanBamParam(flag = scanFlag, what = c("rname", 
                  "pos", "strand", "qwidth", "isize", "mpos"))
            }
            else {
                scanParam = ScanBamParam(flag = scanFlag, what = c("rname", 
                  "pos", "strand", "qwidth", "isize"))
            }
        }
        regions = scanBam(file = paste(path, fileName, sep = "/"), 
            param = scanParam)
        regions = do.call(rbind, lapply(regions, as.data.frame, 
            stringsAsFactors = F))
    }
    else {
        stop("BED files in paired end mode not supported.\n")
    }
    if (!is.null(chr.select) & !bamindex) {
        cat("Selecting", chr.select, "\n")
        regions = regions[regions[, 1] %in% as.vector(chr.select), 
            ]
    }
    cat("Total number of imported first mate reads in properly mapped pairs: ", 
        nrow(regions), "\n", sep = "")
    cat("scanBamFlag: isPaired = T, isProperPair=TRUE , hasUnmappedMate=FALSE, ", 
        "isUnmappedQuery = F, isFirstMateRead = T, isSecondMateRead = F\n", 
        sep = "")
    cat("Mean insertion size: ", mean(abs(regions$isize)), " nt\n", 
        sep = "")
    cat("SD of the insertion size: ", sd(abs(regions$isize)), 
        " nt\n", sep = "")
    cat("Max insertion size: ", max(abs(regions$isize)), " nt\n", 
        sep = "")
    cat("Min insertion size: ", min(abs(regions$isize)), " nt\n", 
        sep = "")
    if (bwa) {
        cat("Paired-end alignment (BAM) files generated by bwa are processed in a different way compared to BAM files generated by bowtie, because in the bwa output the first mate can be either the 'left' or the 'right' mate regardless of their alignment to the plus or the minus strand.\n")
        qwidth = regions[, "qwidth"]
        regions = data.frame(chr = as.character(as.vector(regions$rname)), 
            start = as.numeric(as.vector(regions$pos)), stop = as.numeric(as.vector(regions$mpos)), 
            strand = as.character(as.vector(regions$strand)), 
            isize = as.numeric(as.vector(regions$isize)), stringsAsFactors = F)
        regionsToRev = regions$start > regions$stop
        tmp = regions[regionsToRev, ]$start
        regions[regionsToRev, ]$start = regions[regionsToRev, 
            ]$stop
        regions[regionsToRev, ]$stop = tmp
        regions[, "stop"] = regions[, "stop"] + qwidth - 1 + 
            extend
    }
    else {
        regions = data.frame(chr = as.character(as.vector(regions$rname)), 
            start = as.numeric(as.vector(regions$pos)), stop = as.numeric(as.vector(regions$pos) + 
                as.vector(regions$qwidth) - 1), strand = as.character(as.vector(regions$strand)), 
            isize = as.numeric(as.vector(regions$isize)), stringsAsFactors = F)
        plus = regions$strand == "+"
        regions[plus, "stop"] = regions[plus, "start"] + regions[plus, 
            "isize"] + extend
        regions[!plus, "start"] = regions[!plus, "stop"] + regions[!plus, 
            "isize"] - extend
    }
    regions[, "stop"] = regions[, "stop"] + shift
    regions[, "start"] = regions[, "start"] + shift
    cat("Creating GRange Object...\n")
    regions_GRange = GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, 
        end = regions$stop), strand = regions$strand)
	
	if(is.logical(uniq)){stop("Parameter 'uniq' is not logical anymore, please specify a p-value and see the MEDIPS vignette.")}
	
	if (uniq == 1) {
		cat("Keep only one representative of stacked reads mapping to the same genomic location.\n", sep = "")
		regions_GRange = unique(regions_GRange)
		cat("Number of remaining short reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq < 1 & uniq > 0) {
		max_dup_number = qpois(1 - as.numeric(uniq), length(regions_GRange) / 
			sum(as.numeric(seqlengths(dataset)[chr.select])))
		max_dup_number = max(1, max_dup_number)
		cat("Keep at most ", max_dup_number, 
			" first mate reads mapping to the same genomic location\n", sep = "")		
		uniq_regions = unique(regions_GRange)
		dup_number = countMatches(uniq_regions, regions_GRange)
		dup_number[dup_number > max_dup_number] = max_dup_number
		regions_GRange = rep(uniq_regions, times = dup_number)
		cat("Number of remaining short reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq == 0) {
		cat("Do not correct for potential PCR artefacts (keep all reads).\n", sep = "")
	} else {
		stop("Must specify a valid value for parameter uniq. Please check MEDIPS vignette.")
	}
	strand(regions_GRange) = "*"
	return(regions_GRange)
}

scanBamToGRanges <- function(...) {
	dat <- scanBam(...)[[1]]
	keep <- !is.na(dat$pos)
	GRanges(seqnames=dat$rname[keep],
	ranges=IRanges(start=dat$pos[keep], width=nchar(dat$seq[keep])),
	strand=dat$strand[keep], isize=dat$isize[keep],
	mrnm=dat$mrnm[keep], flag=dat$flag[keep])
}

#######################################
##Read wig file
#######################################
##Input:	wiggle file
##Param:	wiggle.file, path, chromosomes, genome
##Output:	MEDIPSsetObj
##Requires:	rtracklayer
##Modified:	29/10/2012
##Author:	Matthias Lienhard
getMObjectFromWIG <- function(fileName, path, chr.select=NULL,BSgenome){
        cat("Reading wiggle file",fileName,"\n")	
	if(!is.null(chr.select)){
          cat("Select chromosomes",chr.select,"\n")
          sel=GRanges(chr.select,IRanges(1, 536870912))
	  #this function will warn, if type is not bigwig
	  wiggle=rtracklayer::import(paste(path,fileName,sep="/"), asRangedData=FALSE, which=sel)	
	}else{
	  wiggle=rtracklayer::import(paste(path,fileName,sep="/"), asRangedData=FALSE)
	  chr.select=names(seqlengths(wiggle))
	}
	
	dataset=get(ls(paste("package:", BSgenome, sep="")))
	chr_lengths=as.numeric(seqlengths(dataset)[chr.select])
	genome_count=values(wiggle)[,1]
	window_size=width(wiggle)[1]
	#check that all chromosomes are completely covered
        wiggle_chrL=runLength(seqnames(wiggle))
	wiggle_chrN=as.character(runValue(seqnames(wiggle)))
	m=match(chr.select,wiggle_chrN)
	if(any(is.na(m))){
	    stop("ERROR: wiggle file must cover all selected chromosomes\nNot covered chr: ",paste(chr.select[is.na(m)],sep=", "),"\n")
	}
	if(any(wiggle_chrL[m] != ceiling(chr_lengths/window_size))){
	    stop("ERROR: wiggle file must completly cover all selected chromosomes\nNot covered chr: ",paste(wiggle_chrL[m] != ceiling(chr_lengths/window_size),sep=", "),"\n")
	}
	
	#check that values are integer
	tol = .Machine$double.eps^0.5
	if(any(abs(genome_count - round(genome_count))>tol)){
	  stop("ERROR: wiggle file must contain counts of genomic windows, but found floting numbers\n")
	}
	return(new('MEDIPSset', sample_name=fileName,
				                        path_name=path,
				                        genome_name=BSgenome, 
							number_regions=sum(genome_count),
							chr_names=chr.select, 
							chr_lengths=chr_lengths,
							genome_count=genome_count, 
							extend=0,
							shifted=0, 
							window_size=window_size,
							uniq=NA))
}


#####################################
##get types function
#####################################
##Input:	file
##Param:	file, sep
##Output:	vector object
##Modified:	12/10/2011
##Author:	Joern Dietrich
getTypes<-function(file,sep="\t"){
	data=read.table(file,nrows=1,sep=sep,comment.char='')
	types=apply(data,2,typeof)
	return(types)
}

#####################################
##set types function
#####################################
##Input:	vector
##Param:	types
##Output:	vector object
##Modified:	12/10/2011
##Author:	Joern Dietrich
setTypes<-function(types){
	types[2:3]="numeric"
	types[4:5]="NULL"
	return(types)
}

#####################################
##adjust GRange function
#####################################
##Input:	dataframe
##Param:	regions, extend, shift
##Output:	dataframe object
##Modified:	12/10/2011
##Author:	Lukas Chavez, Joern Dietrich
adjustReads<-function(regions, extend, shift){
	if(extend!=0){		
		cat("Extending reads...\n")
		extend.c = pmax(0,extend-regions$stop+regions$start)
		regions$stop[regions$strand=="+"]=regions$stop[regions$strand=="+"]+extend.c[regions$strand=="+"]
		regions$start[regions$strand=="-"]=regions$start[regions$strand=="-"]-extend.c[regions$strand=="-"]	
	}
	
	if(shift!=0){
		cat("Shifting reads...\n")
		regions$start[regions$strand=="+"] = regions$start[regions$strand=="+"]+shift
		regions$stop[regions$strand=="+"] = regions$stop[regions$strand=="+"]+shift
		regions$start[regions$strand=="-"] = regions$start[regions$strand=="-"]-shift
		regions$stop[regions$strand=="-"] = regions$stop[regions$strand=="-"]-shift	
	}
	return(regions)
}
