MEDIPS.CpGenrich <-function(data=data,extend=NULL){
  cat("Preprocessing...\n")
  CpG=RangedData(ranges=IRanges(start=regions_start(data),end=regions_stop(data)),strand=as.character(regions_strand(data)),space=as.character(regions_chr(data)))
  genome_name=genome_name(data)
  genome_chr=chr_names(data)
  chr_lengths=chr_lengths(data)
  dataset=get(ls(paste("package:", genome_name, sep="")))	
  if(!is.null(extend)){
      pos1=CpG[values(CpG)[,"strand"]=="+",] 
      pos2=CpG[values(CpG)[,"strand"]=="-",] 
      ranges(pos1)<- resize(ranges(pos1),extend+1)
      ranges(pos2)<- resize(ranges(pos2),extend+1, start = FALSE)
      CpG=IRanges::rbind(pos1,pos2)  
  }
  ranges(CpG) <-	 restrict(ranges(CpG),+1)
  total=length(genome_chr)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  cat("Calculating CpG density for given regions...\n")
  
  seq=matrix(unlist(IRanges::lapply(CpG,function(x){
  	Sys.sleep(0.1)
        i=which(sort(genome_chr)%in%names(x) )
	setTxtProgressBar(pb, i)
	ranges(x)<-restrict(ranges(x),end=chr_lengths[which(genome_chr %in% names(x))])
   	y=DNAStringSet(getSeq(dataset, names=space(x), start=start(x), end=end(x)))
	c(sum(vcountPattern("CG",y)),sum(vcountPattern("C",y)),sum(vcountPattern("G",y)),sum(width(y)),length(y))
       
  }),use.names=F),ncol=5,nrow=total,byrow=T)
  
  Value=colSums(seq)
  close(pb)
  unused=length(regions_start(data))-Value[5]
  if ( unused!=0 )cat(unused,"unused sequences, limits out of range\n")
  regions.CG=Value[1]
  regions.C=Value[2]
  regions.G=Value[3]
  all.genomic=Value[4]
  regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100
  regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))  
  cat(paste("Calculating CpG density for the reference genome...", genome_name, "\n", sep = " "))	
  CG <- DNAStringSet("CG")
  pdict0 <- PDict(CG)
  params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand"))
  genome.CG=sum(bsapply(params, pdict = pdict0))			
  params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand"), simplify=TRUE)
  alphabet=bsapply(params)
  genome.l=sum(as.numeric(alphabet))
 	
  genome.C=as.numeric(sum(alphabet[2,]))
  genome.G=as.numeric(sum(alphabet[3,]))
  genome.relH=genome.CG/genome.l*100
  genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);
	
  enrichment.score.relH=regions.relH/genome.relH	
  enrichment.score.GoGe=regions.GoGe/genome.GoGe	
 
  return(list(regions.CG=regions.CG, regions.C=regions.C, regions.G=regions.G, regions.relH=regions.relH, regions.GoGe=regions.GoGe, genome.C=genome.C, genome.G=genome.G, genome.CG=genome.CG, genome.relH=genome.relH, genome.GoGe=genome.GoGe, enrichment.score.relH=enrichment.score.relH, enrichment.score.GoGe=enrichment.score.GoGe))  
  
}
