MEDIPS.annotate=function(region=NULL, anno=NULL){
    
    cat(paste("Reading ROI file ", anno,"... \n"), sep="")
    anno=read.table(anno, header=F)
    anno=as.data.frame(anno, stringsAsFactors=F)
    region=as.data.frame(region, stringsAsFactors=F)
    
    regionchr=unique(region[,1])
    result=NULL
    for(chromosome in regionchr){
        if(chromosome%in%unique(anno[,1])){
	   cat(chromosome,"\n")
	   sub_reg=region[region[,1]==chromosome,]
	   sub_anno=anno[anno[,1]==chromosome,]
	   REGION=cbind(as.numeric(factor(sub_reg[,1])),as.numeric(sub_reg[,2]),as.numeric(sub_reg[,3]))
	   ANNO=cbind(as.numeric(factor(sub_anno[,1])),as.numeric(sub_anno[,2]),as.numeric(sub_anno[,3]),c(1:length(sub_anno[,4])))
  	   sub_result=.Call("annotate",as.matrix(ANNO),as.matrix(REGION))
 	   if(dim(sub_result)[1]!=0){
  	        sub_result=cbind(rep(as.character(chromosome),times=length(sub_result[,1])),sub_result[,2],sub_result[,3],as.character(sub_anno[sub_result[,4],4]))
 		result=rbind(result,sub_result)
            }
	}
    }
    colnames(result)=c("chr", "start", "stop", "annotation")
    return(as.data.frame(result, stringsAsFactors=F))
}
