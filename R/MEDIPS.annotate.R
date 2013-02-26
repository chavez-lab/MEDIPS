MEDIPS.annotate=function(region=NULL, anno=NULL){
#MEDIPS.annotate=function(){
	
    .Deprecated("MEDIPS.setAnnotation")
    
#    cat(paste("Reading ROI file ", anno,"... \n"), sep="")
#    anno=read.table(anno, header=F)
#    anno=as.data.frame(anno, stringsAsFactors=F)
#    region=as.data.frame(region, stringsAsFactors=F)
    
#    regionchr=unique(region[,1])
#    result=NULL
#    for(chromosome in regionchr){
#        if(chromosome%in%unique(anno[,1])){
#	   cat(chromosome,"\n")
#	   sub_reg=region[region[,1]==chromosome,]
#	   sub_anno=anno[anno[,1]==chromosome,]
#	   REGION=cbind(as.numeric(factor(sub_reg[,1])),as.numeric(sub_reg[,2]),as.numeric(sub_reg[,3]))
#	   ANNO=cbind(as.numeric(factor(sub_anno[,1])),as.numeric(sub_anno[,2]),as.numeric(sub_anno[,3]),c(1:length(sub_anno[,4])))
# 	   sub_result=.Call("annotate",as.matrix(ANNO),as.matrix(REGION), PACKAGE=MEDIPS)
#	   if(dim(sub_result)[1]!=0){
#	        sub_result=cbind(rep(as.character(chromosome),times=length(sub_result[,1])),sub_result[,2],sub_result[,3],as.character(sub_anno[sub_result[,4],4]))
#		result=rbind(result,sub_result)
#          }
#	}
#    }
#   colnames(result)=c("chr", "start", "stop", "annotation")
#    return(as.data.frame(result, stringsAsFactors=F))
}



##########################################################################
##Function to get annotations by using biomaRt
##########################################################################
##Input:	Names of annotation, dataset and mart
##Param:	mart, dataset. annotation
##Output:	List of annotation matrix 
##Requires:	biomaRt
##Modified:	11/22/2011
##Author:	Joern Dietrich
MEDIPS.getAnnotation<-function(host="www.biomart.org",dataset=c("hsapiens_gene_ensembl","mmusculus_gene_ensembl")[1],annotation=c("TSS","EXON","GENE")[1],tssSz=c(-1000,500),chr=NULL){
	marts=biomaRt::listMarts(host=host)
	mart="ensembl"
	if(host!="www.biomart.org"){
		mart="ENSEMBL_MART_ENSEMBL"
	}
	if(mart%in%marts[,1])
		cat("Getting annotation from database",as.character(marts[marts[,1]==mart,2]),"\n" )
	else
		stop(paste("cannot find mart \"",mart,"\" at ",host,sep=""))		
	biomart <- biomaRt::useMart(mart,host=host)
	dSets=biomaRt::listDatasets(biomart)
	if(dataset%in%dSets[,1])
		cat("Selecting dataset",as.character(dSets[dSets[,1]==dataset,2]),"\n" )
	else
		stop(paste("cannot find dataset \"",dataset,"\" at ",host,sep=""))		
	
	biomart <- biomaRt::useMart(mart, dataset,host=host)
	Annotation=NULL
	f=""
	if(! is.null(chr)) {
		f="chromosome_name"
		for(i in 1:length(chr))#remove the chr prefix
			if(substr(chr[i],1,3)=="chr")
				chr[i]=substr(chr[i],4,nchar(chr[i]))
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
		Annotation=c(Annotation,list("TSS"=data[1:4,]))
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
##Modified:	11/22/2011
##Author:	Joern Dietrich
MEDIPS.setAnnotation<-function(regions, annotation, cnv=F){
	tmp.regions = GRanges(seqnames=regions[,1], ranges=IRanges(start=as.numeric(regions[,2]), end=as.numeric(regions[,3])))	
	ans=NULL
	if(is.data.frame(annotation))annotation=list(annotation=annotation)
	for(anno in annotation){
		anno.data=GRanges(anno[,2], ranges=IRanges(start=anno[,3], end=anno[,4]),ID=anno[,1])
		overlapsM=IRanges::as.matrix(findOverlaps(tmp.regions,anno.data))
		splitL=split(overlapsM[,2],overlapsM[,1])
		maxEle=max(unlist(lapply(splitL,length)))
		tmp.ans=matrix(ncol=maxEle,nrow=length(tmp.regions))
		colnames(tmp.ans)=paste(c(1:maxEle),"_",names(anno)[1],sep="")
		j=rep(names(splitL),sapply(splitL,length))
		k=unlist(lapply(splitL,function(x){1:length(x)}))
		tmp.ans[matrix(c(as.integer(j),as.integer(as.vector(k))),ncol=2)]=as.character(values(anno.data)[overlapsM[,2],1])
		ans=cbind(ans,tmp.ans)		
	}
	if(cnv)
		ans=data.frame(regions, CNV.log2.ratio=as.numeric(ans), stringsAsFactors=F)
	else
		ans=cbind(regions, ans, stringsAsFactors=F)
	return(ans)
}


