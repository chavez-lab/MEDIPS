#######################################
##Read bed file
#######################################
##Input:	tab (|) separated bed file "chr | start | stop | name | score | strand | ..."
##Param:	allignment.file, path, extend, shift, uniq, dataset
##Output:	Granges object
##Requires:	GenomicRanges
##Modified:	06/24/2016
##Author:	Lukas Chavez, Joern Dietrich, Isaac Lopez Moyado

getGRange <-
function (fileName, path = NULL, extend, shift, chr.select = NULL,
    dataset = NULL, uniq = 1e-3, ROI = NULL, isSecondaryAlignment = FALSE, simpleCigar=TRUE)
{
    ext = substr(fileName, nchar(fileName) - 3, nchar(fileName))
    bam = (ext == ".bam" | ext == ".BAM")
    bamindex = bam & file.exists(paste(path, "/", fileName, ".bai",
        sep = ""))
    if (bam) {
        scanFlag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = isSecondaryAlignment)
        if (bamindex & (!is.null(chr.select) | !is.null(ROI))) {
            if (!is.null(ROI)) {
                message("Reading bam alignment ", fileName, "\n considering ROIs using bam index", appendLF=T)
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
                message("Reading bam alignment ", fileName, "\n considering ",
                  paste(chr.select, collapse="\t"), " using bam index", appendLF=T)
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
            scanParam = ScanBamParam(flag = scanFlag, simpleCigar= simpleCigar, what = c("rname",
                "pos", "strand", "qwidth", "isize", "mpos"), which = sel)
        }
        else {
            message("Reading bam alignment ", fileName, appendLF=T)
            scanParam = ScanBamParam(flag = scanFlag, simpleCigar= simpleCigar, what = c("rname",
                "pos", "strand", "qwidth", "isize", "mpos"))
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
        message("Reading bed alignment ", fileName, appendLF=T)
        regions = read.table(paste(path, fileName, sep = "/"),
            sep = "\t", header = FALSE, row.names = NULL, comment.char = "",
            colClasses = c("character", "numeric", "numeric",
                "NULL", "NULL", "character"))
        names(regions) = c("chr", "start", "stop", "strand")
    }
    if (!is.null(chr.select) & !bamindex) {
        message("Selecting ", paste(chr.select, collapse="\t"), appendLF=T)
        regions = regions[regions[, 1] %in% as.vector(chr.select),
            ]
    }
    message("Total number of imported short reads: ", nrow(regions), appendLF=T)
    regions = adjustReads(regions, extend, shift)
    message("Creating GRange Object...", appendLF=T)
    regions_GRange = GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start,
        end = regions$stop), strand = regions$strand)

    if(is.logical(uniq)){stop("Parameter 'uniq' is not logical anymore, please specify a p-value and see the MEDIPS vignette.")}
    if (uniq == 1) {
		message("Keep at most one 1 read mapping to the same genomic location.", appendLF=T)
		regions_GRange = unique(regions_GRange)
		message("Number of remaining reads: ", length(regions_GRange), appendLF=T)
	} else if (uniq < 1 & uniq > 0) {
		max_dup_number = qpois(1 - as.numeric(uniq), length(regions_GRange) /
			sum(as.numeric(seqlengths(dataset)[chr.select])))
		max_dup_number = max(1, max_dup_number)
		message("Keep at most ", max_dup_number,
			" read(s) mapping to the same genomic location", appendLF=T)
		uniq_regions = unique(regions_GRange)
		dup_number = countMatches(uniq_regions, regions_GRange)
		dup_number[dup_number > max_dup_number] = max_dup_number
		regions_GRange = rep(uniq_regions, times = dup_number)
		message("Number of remaining reads: ", length(regions_GRange), appendLF=T)
	} else if (uniq == 0) {
		message("Do not correct for potential PCR artefacts (keep all reads).", appendLF=T)
	} else {
		stop("Must specify a valid value for parameter uniq. Please check MEDIPS vignette.")
	}
	strand(regions_GRange) = "*"
	return(regions_GRange)
}

getPairedGRange <-
function (fileName, path = NULL, extend, shift, chr.select = NULL,
    dataset = NULL, uniq = 1e-3, ROI = NULL, isSecondaryAlignment = FALSE, simpleCigar=TRUE)
{
    ext = substr(fileName, nchar(fileName) - 3, nchar(fileName))
    bam = (ext == ".bam" | ext == ".BAM")
    bamindex = bam & file.exists(paste(path, "/", fileName, ".bai",
        sep = ""))
    if (bam) {
        scanFlag = scanBamFlag(isPaired = T, isProperPair = TRUE,
            hasUnmappedMate = FALSE, isUnmappedQuery = F, isFirstMateRead = T,
            isSecondMateRead = F, isSecondaryAlignment = isSecondaryAlignment)
        if (bamindex & (!is.null(chr.select) | !is.null(ROI))) {
            if (!is.null(ROI)) {
                message("Reading bam alignment ", fileName, "\n considering ROIs using bam index", appendLF=T)
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
                message("Reading bam alignment", fileName, "\n considering ",
                  paste(chr.select, collapse="\t"), " using bam index", appendLF=T)
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
        scanParam = ScanBamParam(flag = scanFlag, simpleCigar= simpleCigar, what = c("rname",
		    "pos", "strand", "qwidth", "isize", "mpos"), which = sel)
        }
        else {
            message("Reading bam alignment ", fileName, appendLF="TRUE")
            scanParam = ScanBamParam(flag = scanFlag, simpleCigar = simpleCigar, what = c("rname",
                  "pos", "strand", "qwidth", "isize", "mpos"))
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
        message("Selecting", paste(chr.select, collapse="\t"), appendLF=T)
        regions = regions[regions[, 1] %in% as.vector(chr.select),
            ]
    }
    message("Total number of imported first mate reads in properly mapped pairs: ",
        nrow(regions), appendLF=T)
    message("scanBamFlag: isPaired = T, isProperPair=TRUE , hasUnmappedMate=FALSE, ",
        "isUnmappedQuery = F, isFirstMateRead = T, isSecondMateRead = F",
        appendLF=T)
    message("Mean insertion size: ", mean(abs(regions$isize)), " nt",
        appendLF=T)
    message("SD of the insertion size: ", sd(abs(regions$isize)),
        " nt", appendLF=T)
    message("Max insertion size: ", max(abs(regions$isize)), " nt",
        appendLF=T)
    message("Min insertion size: ", min(abs(regions$isize)), " nt",
        appendLF=T)

   qwidth = regions[, "qwidth"]
   regions = data.frame(chr = as.character(as.vector(regions$rname)),
   start = as.numeric(as.vector(regions$pos)), stop = as.numeric(as.vector(regions$mpos)),
   strand = as.character(as.vector(regions$strand)),
   isize = as.numeric(as.vector(regions$isize)), stringsAsFactors = F)

   regionsToRev = regions$start > regions$stop
   regions[regionsToRev, ]$start = regions[regionsToRev,]$stop
   regions[, "stop"] = regions[, "start"] + abs(regions[, "isize"]) - 1

   if(extend!=0){warning("The extend parameter will be neglected, because the actual DNA fragment length is known in paired-end data.\n")}
   if(shift!=0){warning("The shift parameter will be neglected, because the actual DNA fragment position is known in paired-end data.\n")}

    message("Creating GRange Object...", appendLF=T)
    regions_GRange = GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start,
        end = regions$stop), strand = regions$strand)

	if(is.logical(uniq)){stop("Parameter 'uniq' is not logical anymore, please specify a p-value and see the MEDIPS vignette.")}

	if (uniq == 1) {
		message("Keep at most 1 read mapping to the same genomic location.", appendLF=T)
		regions_GRange = unique(regions_GRange)
		message("Number of remaining short reads: ", length(regions_GRange),
			appendLF=T)
	} else if (uniq < 1 & uniq > 0) {
		max_dup_number = qpois(1 - as.numeric(uniq), length(regions_GRange) /
			sum(as.numeric(seqlengths(dataset)[chr.select])))
		max_dup_number = max(1, max_dup_number)
		message("Keep at most ", max_dup_number,
			" first mate read(s) mapping to the same genomic location", appendLF=T)
		uniq_regions = unique(regions_GRange)
		dup_number = countMatches(uniq_regions, regions_GRange)
		dup_number[dup_number > max_dup_number] = max_dup_number
		regions_GRange = rep(uniq_regions, times = dup_number)
		message("Number of remaining short reads: ", length(regions_GRange),
			appendLF=T)
	} else if (uniq == 0) {
		message("Do not correct for potential PCR artefacts (keep all reads).", appendLF=T)
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
        message("Reading wiggle file ", fileName, appendLF=T)
	if(!is.null(chr.select)){
          message("Select chromosomes ", paste(chr.select, collapse="\t"), appendLF=T)
          sel=GRanges(chr.select,IRanges(1, 536870912))
	  #this function will warn, if type is not bigwig
	  wiggle=rtracklayer::import(paste(path,fileName,sep="/"), which=sel)
	}else{
	  wiggle=rtracklayer::import(paste(path,fileName,sep="/"))
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
		message("Extending reads...", appendLF=T)
		extend.c = pmax(0,extend-regions$stop+regions$start)
		regions$stop[regions$strand=="+"]=regions$stop[regions$strand=="+"]+extend.c[regions$strand=="+"]
		regions$start[regions$strand=="-"]=regions$start[regions$strand=="-"]-extend.c[regions$strand=="-"]
	}

	if(shift!=0){
		message("Shifting reads...", appendLF=T)
		regions$start[regions$strand=="+"] = regions$start[regions$strand=="+"]+shift
		regions$stop[regions$strand=="+"] = regions$stop[regions$strand=="+"]+shift
		regions$start[regions$strand=="-"] = regions$start[regions$strand=="-"]-shift
		regions$stop[regions$strand=="-"] = regions$stop[regions$strand=="-"]-shift
	}
	return(regions)
}
