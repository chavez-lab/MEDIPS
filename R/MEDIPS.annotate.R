##########################################################################
##Function to get annotations by using biomaRt
##########################################################################
##Input:	Names of annotation, dataset and mart
##Param:	mart, dataset. annotation
##Output:	List of annotation matrix 
##Requires:	biomaRt
##Modified:	01/19/2016
##Author:	Joern Dietrich, Matthias Lienhard, Lukas Chavez
MEDIPS.getAnnotation<-function(host="www.ensembl.org",dataset=c("hsapiens_gene_ensembl","mmusculus_gene_ensembl")[1],annotation=c("TSS","EXON","GENE")[1],tssSz=c(-1000,500),chr=NULL){

	marts=biomaRt::listMarts(host=host)
	mart="ENSEMBL_MART_ENSEMBL"
	if(mart%in%marts[,1])
		cat("Getting annotation from database",as.character(marts[marts[,1]==mart,2]),"\n" )
	else
		stop(paste("cannot find mart \"",mart,"\" at ",host,sep=""))		
	biomart <- biomaRt::useMart(mart, dataset=dataset, host=host)
	dSets=biomaRt::listDatasets(biomart)
	if(dataset%in%dSets[,1])
		cat("Selecting dataset",as.character(dSets[dSets[,1]==dataset,2]),"\n" )
	else
		stop(paste("cannot find dataset \"",dataset,"\" at ",host,sep=""))		
	
	biomart <- biomaRt::useMart(mart, dataset,host=host)
	Annotation=NULL
	if(! is.null(chr)) {
		f="chromosome_name"
		for(i in 1:length(chr))#remove the chr prefix
			if(substr(chr[i],1,3)=="chr")
				chr[i]=substr(chr[i],4,nchar(chr[i]))
	}else{
		f=""
		chr=""
	} 

	if("GENE"%in%annotation){
		cat("...getting gene annotation\n")
		data <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),mart=biomart,filters=f, values=chr)
		data$chromosome_name=paste("chr",data$chromosome_name,sep="")
		names(data)=c("id", "chr", "start", "end")
		Annotation=list("Gene"=data)
	}
	
	if("TSS"%in%annotation){
		cat("...getting TSS annotation\n")
		data <- getBM(attributes=c("ensembl_transcript_id","chromosome_name","transcript_start","transcript_end","strand"),mart=biomart,filters=f, values=chr)
		fs=data$strand == 1 #forward strand
		data[fs ,"transcript_end"  ]=data[fs,"transcript_start"]+tssSz[2]
		data[fs ,"transcript_start"]=data[fs,"transcript_start"]+tssSz[1]
		data[!fs ,"transcript_start"]=data[!fs,"transcript_end"]-tssSz[2]
		data[!fs ,"transcript_end"  ]=data[!fs,"transcript_end"]-tssSz[1]
		data$ensembl_transcript_id=paste("TSS",data$ensembl_transcript_id,sep="_")
		data$chromosome_name=paste("chr",data$chromosome_name,sep="")
		names(data)=c("id", "chr", "start", "end","strand")
		Annotation=c(Annotation,list("TSS"=data[,1:4]))
	}
	
	if("EXON"%in%annotation){
		cat("...getting exon annotation\n")
		data <- getBM(attributes=c("ensembl_exon_id","chromosome_name","exon_chrom_start","exon_chrom_end"),mart=biomart,filters=f, values=chr)
		data$chromosome_name=paste("chr",data$chromosome_name,sep="")
		names(data)=c("id", "chr", "start", "end")
		Annotation=c(Annotation,list("EXON"=data))
	}
	
	
	return(Annotation)
}


##########################################################################
##Function to add annotation to result data regarding their position 
##########################################################################
##Input:	Result matrix or data.frame and list of annotation
##Param:	regions. annotation
##Output:	List of annotation matrix 
##Requires:	IRanges
##Modified:	11/01/2015
##Author:	Joern Dietrich, Matthias Lienhard
MEDIPS.setAnnotation<-function(regions, annotation, cnv=F){

	tmp.regions = GRanges(seqnames=regions[,1], ranges=IRanges(start=as.numeric(regions[,2]), end=as.numeric(regions[,3])))
	ans=NULL
	if(is.data.frame(annotation))annotation=list(annotation=annotation)
	for(n in names(annotation)){ #um an den namen für die message zu kommen
	anno=annotation[[n]]
	anno.data=GRanges(anno[,2], ranges=IRanges(start=anno[,3], end=anno[,4]),ID=anno[,1])
	overlapsM=IRanges::as.matrix(findOverlaps(tmp.regions,anno.data))
	splitL=split(overlapsM[,2],overlapsM[,1])
	maxEle=max(c(0,unlist(lapply(splitL,length)))) #vermeidet die Warnung
	if(maxEle==0){ #mache nichts für diese anntoation und mach mit dem Nächsten weiter
		message("no \"",n,"\"annotation overlap the provided regions")
	}else{
		tmp.ans=matrix(ncol=maxEle,nrow=length(tmp.regions))
		colnames(tmp.ans)=paste(c(1:maxEle),"_",names(anno)[1],sep="")
		j=rep(names(splitL),sapply(splitL,length))
		k=unlist(lapply(splitL,function(x){1:length(x)}))
		tmp.ans[matrix(c(as.integer(j),as.integer(as.vector(k))),ncol=2)]=as.character(values(anno.data)[overlapsM[,2],1])
		ans=cbind(ans,tmp.ans)
		}
	}
	if(cnv){
		ans=data.frame(regions, CNV.log2.ratio=as.numeric(ans), stringsAsFactors=F)
	}else{
		ans=cbind(regions, ans, stringsAsFactors=F)
	}
	return(ans)
}


