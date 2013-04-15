#######################################
##Read bed file
#######################################
##Input:	tab (|) separated bed file "chr | start | stop | name | score | strand | ..."
##Param:	allignment.file, path, extend, shift, uniq
##Output:	Granges object
##Requires:	GenomicRanges
##Modified:	11/10/2011
##Author:	Lukas Chavez, Joern Dietrich
getGRange <- function(fileName, path=NULL,extend, shift, chr.select=NULL, uniq=FALSE, ROI=NULL){
	if(!is.null(ROI) & !is.null(chr.select)){
		cat("Selecting ROIs from ",chr.select,"\n")		
		ROI=ROI[ROI[,1] %in% chr.select,]
	}
	ext=substr(fileName,nchar(fileName)-3,nchar(fileName))
	bam=( ext==".bam" | ext ==".BAM" )
	bamindex=bam & file.exists(paste(path,"/", fileName,".bai", sep=""))
	if (bam){

		scanFlag = scanBamFlag(isUnmappedQuery = F)
		
		if(bamindex & (!is.null(chr.select) | !is.null(ROI)) ){#read bam with index			
			if(!is.null(ROI)){
			    cat("Reading bam alignment",fileName,"\n considering ROIs using bam index\n")	
			    if(!is.null(extend)) {
				ROI[,2]=ROI[,2]-extend
				ROI[,3]=ROI[,3]+extend
			    }
			    if(!is.null(shift)){				
				ROI[,2]=ROI[,2]-shift
				ROI[,3]=ROI[,3]-shift
			    }
                            #sel=GRanges(ROI[,1],IRanges(start=ROI[,2], end=ROI[,3]))
			    sel=GRanges(chr.select,IRanges(1, 536870912))
			}else{
			    cat("Reading bam alignment",fileName,"\n considering ",chr.select," using bam index\n")		
                            sel=GRanges(chr.select,IRanges(1, 536870912))
			}
			scanParam=ScanBamParam(flag=scanFlag, what = c("rname", "pos", "strand", "qwidth"), which=sel)			
		}else {#read bam without index
			cat("Reading bam alignment",fileName,"\n")		
			scanParam=ScanBamParam(flag=scanFlag, what = c("rname", "pos", "strand", "qwidth"))
		}
		regions = scanBam(file=paste(path, fileName, sep="/"), param=scanParam)
		regions = do.call(rbind,lapply(regions, as.data.frame, stringsAsFactors=F))
		regions = data.frame(chr=as.character(as.vector(regions$rname)), start=as.numeric(as.vector(regions$pos)), stop=as.numeric(as.vector(regions$pos)+as.vector(regions$qwidth)-1), strand=as.character(as.vector(regions$strand)), stringsAsFactors=F)
		
	}else{#read from bed file
		cat("Reading bed alignment",fileName,"\n")		
		regions=read.table(paste(path, fileName, sep="/"), sep='\t', header=FALSE, row.names=NULL, comment.char='', colClasses=c("character", "numeric", "numeric", "NULL", "NULL", "character"))
		names(regions)	=c("chr", "start", "stop", "strand")
	}

	if(!is.null(chr.select)& !bamindex){
		cat("Selecting",chr.select,"\n")
		regions = regions[regions[,1] %in% as.vector(chr.select),]
	}		
	cat("Total number of imported short reads: ", nrow(regions), "\n", sep="")		
	
	regions = adjustReads(regions, extend, shift)
	
	cat("Creating GRange Object...\n")
	regions_GRange = GRanges(seqnames=regions$chr, ranges=IRanges(start=regions$start, end=regions$stop), strand=regions$strand)
	if(uniq){
		cat(paste("Extract unique regions...\n", sep=""))
		regions_GRange=unique(regions_GRange)
		cat("Number of unique short reads: ", length(regions_GRange), "\n", sep="")
	}
	strand(regions_GRange)="*"
	return(regions_GRange)
}

getPairedGRange <- function(fileName, path=NULL,extend, shift, chr.select=NULL, uniq=FALSE, ROI=NULL){
	if(!is.null(ROI) & !is.null(chr.select)){
		cat("Selecting ROIs from ",chr.select,"\n")		
		ROI=ROI[ROI[,1] %in% chr.select,]
	}
        ext=substr(fileName,nchar(fileName)-3,nchar(fileName))
        bam=( ext==".bam" | ext ==".BAM" )
        bamindex=bam & file.exists(paste(path,"/", fileName,".bai", sep=""))
        if (bam){
		scanFlag = scanBamFlag(isPaired=T, isProperPair=TRUE ,
			hasUnmappedMate=FALSE, isUnmappedQuery = F, isFirstMateRead = T,
			isSecondMateRead = F)
		#scanFlag = scanBamFlag(isPaired=T, isUnmappedQuery = F, isFirstMateRead = T, isSecondMateRead = F)
                if(bamindex & (!is.null(chr.select) |! is.null(ROI))){#read bam with index
			if(!is.null(ROI)){
			    cat("Reading bam alignment",fileName,"\n considering ROIs using bam index\n")	
			    if(!is.null(extend)) {
				ROI[,2]=ROI[,2]-extend
				ROI[,3]=ROI[,3]+extend
			    }
			    if(!is.null(shift)){				
				ROI[,2]=ROI[,2]-shift
				ROI[,3]=ROI[,3]-shift
			    }
                            sel=GRanges(ROI[,1],IRanges(start=ROI[,2], end=ROI[,3]))
			}else{
			    cat("Reading bam alignment",fileName,"\n considering ",chr.select," using bam index\n")		
                            sel=GRanges(chr.select,IRanges(1, 536870912))
			}
			scanParam=ScanBamParam(flag=scanFlag, what = c("rname", "pos", "strand", "qwidth", "isize"), which=sel)
                }else {#read bam without index
                        cat("Reading bam alignment",fileName,"\n")
			scanParam=ScanBamParam(flag=scanFlag, what = c("rname", "pos", "strand", "qwidth", "isize"))
                }
                regions = scanBam(file=paste(path, fileName, sep="/"), param=scanParam)
		regions = do.call(rbind,lapply(regions, as.data.frame, stringsAsFactors=F))		

        }else{#read from bed file
                cat("BED files in paired end mode not allowed.\n")
		return(NULL)
        }

	if(!is.null(chr.select)& !bamindex){
                cat("Selecting",chr.select,"\n")
                regions = regions[regions[,1] %in% as.vector(chr.select),]
        }

	cat("Total number of imported first mate reads: ", nrow(regions), "\n", sep="")
	#cat("scanBamFlag: isPaired = T, isUnmappedQuery = F, isFirstMateRead = T, isSecondMateRead = F\n", sep="")
	cat("scanBamFlag: isPaired = T, isProperPair=TRUE , hasUnmappedMate=FALSE, ",
		"isUnmappedQuery = F, isFirstMateRead = T, isSecondMateRead = F\n", sep="")

        if(extend!=0 | shift!=0){
		cat("In paired end mode no region adjusting allowed\n")
	}
	
	##Get mean distance and sd
	cat("Mean insertion size: ", mean(abs(regions$isize)), " nt\n", sep="")
	cat("SD of the insertion size: ", sd(abs(regions$isize)), " nt\n", sep="")
	cat("Max insertion size: ", max(abs(regions$isize)), " nt\n", sep="")
	cat("Min insertion size: ", min(abs(regions$isize)), " nt\n", sep="")
        
	prev=substr(fileName,1,nchar(fileName)-4)
	png(filename=paste(paste(path, prev, sep="/"), "IS.distribution.png", sep="."))
	hist(abs(regions$isize), breaks=100)
	dev.off()

	##Filter for F3 reads having a mate in reasonable distance
	t.l=quantile(abs(regions$isize), probs=c(0.99))[[1]]
	cat("Removing reads with insertion size>: ", t.l, "\n", sep="")	
	regions = regions[regions$isize<=t.l,]		
	cat("Number of remaining regions: ", nrow(regions), "\n", sep="")
		
	#Extend the regions according to their insertion size and create data frame
	cat("Extend first mate according to insert size + 35bp into sequencing direction.\n", sep="")
	mateRL = 35  #read length of the 2nd mate
	regions = data.frame(chr=as.character(as.vector(regions$rname)), start=as.numeric(as.vector(regions$pos)), stop=as.numeric(as.vector(regions$pos)+as.vector(regions$qwidth)-1), strand=as.character(as.vector(regions$strand)), isize=as.numeric(as.vector(regions$isize)), stringsAsFactors=F)
	regions[which(regions$strand=="+"), 3] = regions[which(regions$strand=="+"), 3] + regions[which(regions$strand=="+"), 5] + mateRL
	regions[which(regions$strand=="-"), 2] = regions[which(regions$strand=="-"), 2] + regions[which(regions$strand=="-"), 5] - mateRL

        cat("Creating GRange Object...\n")
        regions_GRange = GRanges(seqnames=regions$chr, ranges=IRanges(start=regions$start, end=regions$stop), strand=regions$strand)
        if(uniq){
                cat(paste("Extract unique regions...\n", sep=""))
                regions_GRange=unique(regions_GRange)
                cat("Number of unique short reads: ", length(regions_GRange), "\n", sep="")
        }
	strand(regions_GRange)="*"
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
	if(length(genome_count) == sum(ceiling(chr_lengths/window_size))){
	  pos=0
	  for(chr_idx in 1:length(chr.select)){
	    chr_in_wig=unique(as.vector(seqnames(wiggle[(pos+1):(pos+ceiling(chr_lengths[chr_idx]/window_size))] )))
	    if(length(chr_in_wig)!=1 | chr_in_wig != chr.select[chr_idx]){
		cat("ERROR: wiggle file must completely cover all selected chromosomes\nNon-conformance found in ",chr.select[chr_idx],"\n")
	  	return(NULL)
	    }
	    pos=pos+ceiling(chr_lengths[chr_idx]/window_size)
	  }
	}else{
	  stop("ERROR: wiggle file must completely cover all selected chromosomes\n")
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
##get GRange function
#####################################
##Input:	file
##Param:	filename, path, extend
##Output:	GRange object
##Modified:	12/10/2011
##Author:	Joern Dietrich
#getGRange<-function(fileName, path, extend, shift, chr.select){
#	 if(extend!=0 & shift!=0){
#                stop("One of the parameters extend or shift has to be 0!")
#	}
#	cat(paste("Reading file ", fileName, " in ", path, "...\n", sep=""))
#	GRange.Reads = suppressWarnings(try(MEDIPS.Bam2GRanges(fileName, path, extend, shift, chr.select),silent=T))
#	if(typeof(GRange.Reads)!="S4"){
#		GRange.Reads = suppressWarnings(try(MEDIPS.Bed2Granges(fileName, path, extend, shift, chr.select),silent=T))	
#	}
#	return(GRange.Reads)
#}

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
		####???????extend.c = abs(regions.l - extend)
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
