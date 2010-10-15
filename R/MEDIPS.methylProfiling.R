MEDIPS.methylProfiling=function (data1 = NULL, data2 = NULL, input = NULL, ROI_file = NULL, frame_size = NULL, math = mean, step = NULL, select = 2, chr = NULL, transf=T)
{
    if(!is.null(ROI_file) & !is.null(chr)){stop("You cannot specify a ROI file and a selected chr.")}
    if (class(data1) != "MEDIPSset") stop("Must specify a MEDIPS.SET object at data1.")
    
     cat("Preprocessing...\n")
    data1@genome_raw = data1@genome_raw/(number_regions(data1)/1e+06)
      
    if(!is.null(input)){
        if (class(input) != "MEDIPSset") stop("Must specify a MEDIPS.SET object at input.")
        input = input@genome_raw/(number_regions(input)/1e+06)
    }
    if(!is.null(data2)){
        if (class(data2) != "MEDIPSset") stop("Must specify a MEDIPS.SET object at data2.")  
	if (!identical(chr_names(data1), chr_names(data2))) stop("MEDIPS.SETs contain different number of chromosomes!")
        if (bin_size(data1) != bin_size(data2)) stop("MEDIPS.SETs have different bin sizes!")
        data2@genome_raw = data2@genome_raw/(number_regions(data2)/1e+06)
    }
        
    if(!is.null(ROI_file) | !is.null(chr)) {
        if(!is.null(chr)){
            indices = which(genome_chr(data1) %in% chr)
            data1@chr_lengths = chr_lengths(data1)[chr_names(data1) %in% chr]
            if (length(indices) == 0) {stop("Stated chromosome not available in MEDIPS.SET")}
	}
        else{
            cat("Reading ROIs...\n")
	    if(is.matrix(ROI_file) | is.data.frame(ROI_file)){ROI=ROI_file}
	    else{ROI = read.table(ROI_file)}
	    if(dim(ROI)[2]<4){stop("Invalid ROI format!")}
            if (length(unique(ROI[, 1])) == 1) {ROI_chromosomes = unique(ROI[, 1])}
            else {
	        ROI_chromosomes = mixedsort(unique(ROI[, 1]))
	    }
	    if(!identical(chr_names(data1), as.character(ROI_chromosomes))) {
	        indices = which(genome_chr(data1) %in% unique(ROI[,1]))
                if (length(indices) == 0) {stop("MEDIPS.SET and ROIs have no chromosomes in common!")}
                data1@chr_lengths = chr_lengths(data1)[chr_names(data1) %in%ROI_chromosomes]
            }
            else {indices = NULL}
        }
        if (!is.null(indices)) {
            cat("Extract data according to given ROI...\n")
            data1@genome_chr = genome_chr(data1)[indices]
            data1@genome_pos = genome_pos(data1)[indices]
            data1@genome_raw = genome_raw(data1)[indices]
            data1@genome_norm = genome_norm(data1)[indices]
            data1@genome_CF = genome_CF(data1)[indices]
            data1@chr_names = unique(genome_chr(data1))
            if (!is.null(input)) {
                input=input[indices]
            }
	    if(!is.null(data2)){
	        data2@genome_raw = genome_raw(data2)[indices]
		data2@genome_norm = genome_norm(data2)[indices]		
	    }
	}
	if(!is.null(ROI_file)){
	    uni=unique(ROI[,1])
            exclude=uni[which(!uni%in%chr_names(data1))]
            ROI2=cbind(as.numeric(factor(as.character(ROI[,1]),exclude = exclude)),as.numeric(ROI[,2]),as.numeric(ROI[,3]))
	    bin_counts = sapply(chr_lengths(data1), function(x) {ceiling(x/bin_size(data1))})
            chr_binposition = sapply(1:length(bin_counts), function(x) {sum(bin_counts[1:x])})
            if(!is.null(data2)){cat("Differential methylation will be calculated on the ROI data set\n")}
	    else{cat("Methylation profile will be calculated on the ROI data set\n")}
	    
	    matrix=suppressWarnings(.Call("roiprofile",input,as.numeric(select),as.matrix(ROI2),as.integer(chr_binposition), data1,data2 ,environment(wilcox.test),wilcox.test,environment(var),var,environment(math),math,t.test,environment(t.test),as.numeric(factor(chr_names(data1)))))
 	    rownames(matrix)=ROI[,4]
	}
    }
    if(is.null(ROI_file)){
        if(!is.null(data2) & is.null(chr)){
	    if (is.null(step)) {
                cat("Differential methylation will be calculated on the full data set for all chromosomes with distinct frames of size ",frame_size,":\n")
            }
	    else{
                cat("Differential methylation will be calculated on the full data set for all chromosomes with frame size ",
                frame_size, " and step size",step,":\n")
            }
	}
	else if (is.null(chr)){
	    if (is.null(step)) {
                cat("Methylation profile will be calculated on the full data set for all chromosomes with distinct frames of size ",frame_size,":\n")
            }
	    else{
                cat("Methylation profile will be calculated on the full data set for all chromosomes with frame size ",
                frame_size, " and step size",step,":\n")
            }
	}
	else {
	    if(!is.null(data2)){
	        if (is.null(step)) {cat("Differential methylation will be calculated on the selected chromosome",chr,"with distinct frames of size ", frame_size,":\n")}
                else {cat("Differential methylation will be calculated on the selected chromosome ",chr, " with frame size ", frame_size, " and  step size",step,":\n")}
	    }
	    else{
	        if (is.null(step)) {cat("Methylation profile will be calculated on the selected chromosome",chr,"with distinct frames of size ", frame_size,":\n")}
            else {cat("Methylation profile will be calculated on the selected chromosome ",chr, " with frame size ", frame_size, " and  step size",step,":\n")}
	    }
	}
      	matrix=suppressWarnings(.Call("profile",as.numeric(select),input,as.double(step),as.double(frame_size),data1,data2,environment(wilcox.test),wilcox.test,environment(var),var,environment(math),math,t.test,environment(t.test)))
    }
    cat("Postprocessing...\n")
    colnames(matrix)=c("chr","start","stop","length","coupling","input","rpm_A","rpm_B","rms_A","rms_B","ams_A","ams_B","var_A","var_B","var_co_A","var_co_B","ratio","p.value.wilcox","p.value.ttest")
    matrix=as.data.frame(matrix)
    myreplace=function(x){
        new=sapply(x,function(y){
            if(!is.na(y)){
	        if( y == -1){y=as.numeric("Inf")}
	    }
            return(y);	     
        })
        return(new);
    }
    matrix2=apply(matrix,1,myreplace)
    matrix=as.data.frame(t(matrix2))
    if(!is.null(ROI_file)){
           matrix[,1]=ROI[,1]
    }
    else{
        matrix[,1]=chr_names(data1)[matrix[,1]]
    }
    if(transf){
	matrix[,9] = MEDIPS.transform(as.numeric(matrix[,9 ]))
	matrix[,11] = MEDIPS.transform(as.numeric(matrix[,11 ]))    
	if(!is.null(data2)){  
		matrix[,10] = MEDIPS.transform(as.numeric(matrix[,10 ]))    
		matrix[,12] = MEDIPS.transform(as.numeric(matrix[,12 ]))
    	}
    }

    gc()
    return(matrix)
}
