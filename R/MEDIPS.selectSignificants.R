##Function identifies specified columns (input, rpm, rms values...) and selects genomic regions that fulfill the defined constraints.

MEDIPS.selectSignificants = function(frames=NULL, input=T, control=T, up=1.333333, down=0.75, p.value=0.01, quant=0.9){
	
	if(is.null(frames)){stop("Must specify a frame set.")}	
				
	##Subtract right columns by header names
	########################################
	numberCols=length(colnames(frames))
	if(numberCols==0){stop("No colnames given. Can't extract right columns.")}
	value_input=NULL
	value_rpm_A=NULL
	value_rpm_B=NULL
	value_rms_A=NULL
	value_rms_B=NULL	
	for(i in 1:numberCols){
		if(colnames(frames)[i]=="input"){value_input=i}
		else if(colnames(frames)[i]=="ratio"){value_ratio=i}
		else if(colnames(frames)[i]=="rpm_A"){value_rpm_A=i}
		else if(colnames(frames)[i]=="rpm_B"){value_rpm_B=i}
		else if(colnames(frames)[i]=="rms_A"){value_rms_A=i}
		else if(colnames(frames)[i]=="rms_B"){value_rms_B=i}
		else if(colnames(frames)[i]=="p.value.wilcox"){value_pvalue.wilcox=i}
		else if(colnames(frames)[i]=="p.value.ttest"){value_pvalue.ttest=i}
	}	
	
	##Input processing
	if(input==T){
		if(!is.null(value_input)){
			print(paste("Processing input distribution..."), quote=F)
			t=quantile(as.numeric(frames[,value_input]), probs=c(quant), na.rm=T)[[1]]
			print(paste("Done."), quote=F)
		}
		else{stop("No input column identified.")}
	}	
	else{
		##Set the minimum rpm value to a fix number.
		t=input
	}	
		
	##Start filtering
	######################################	
	print(paste("Total number of frames: ", length(frames[,1]), sep=""), quote=F)
						
	##Filter for frames where A!=0 | B!=0	
	frames=frames[(as.numeric(frames[,value_rpm_A])!=0 | as.numeric(frames[,value_rpm_B])!=0),]
	print(paste("Number of frames where control or treatment !=0: ", length(frames[,1])), quote=F)
			
	#Filter for p.values
	frames=frames[(as.numeric(frames[,value_pvalue.wilcox])<=p.value | as.numeric(frames[,value_pvalue.ttest])<=p.value),]
	print(paste("Remaining number of frames with p.value<=", p.value, ": ", length(frames[,1]), sep=""), quote=F)
				 
	##Filter for ratios (control vs. control and control vs. input)
	#selected[,value_A]==0 --> 0 
	#selected[,value_B]==0 --> Inf
	#0/0--> NaN
	#-->Both cases are caught by the imbalance commands!		
	#####################################################
	
	##Control up:
	if(control==T){		
		##Filter for MeDIP control vs. MeDIP treatment
		frames=frames[as.numeric(frames[,value_ratio])>=up,]
		print(paste("Remaining number of frames where control/treatment ratio >= ", up, ": ", length(frames[,1]), sep=""), quote=F)
				
		##Filter for MeDIP control rpm signals higher than given quantile of the Input distribution
		if(input==T){
			print(paste("Estimated rpm threshold for input quantile ", quant, " is: ", t, sep=""), quote=F)
		}
		else{
			print(paste("Given rpm threshold is: ", t, sep=""), quote=F)
		}
		frames=frames[as.numeric(frames[,value_rpm_A])>=t, ]
		print(paste("Remaining number of frames with control rpm >=", t, ": ", length(frames[,1]), sep=""), quote=F)		
		
		##Filter for rpm MeDIP control vs. Input
		if(input==T){
			frames=frames[((as.numeric(frames[,value_rpm_A])/as.numeric(frames[,value_input]))>=up), ]
			print(paste("Remaining number of frames with control/input ratio>=", up, ": ", length(frames[,1]), sep=""), quote=F)
		}
		print(paste("[Note: There are ", length(frames[,1][frames[,value_pvalue.ttest]==0 & frames[,value_pvalue.wilcox]==0]), " frames associated to a p-value==0.]", sep=""), quote=F)
		print(paste("[Note: There are ", length(frames[,1][frames[,value_ratio]=="Inf"]), " frames, where control/treatment ratio = Inf (i.e. treatment==0).]", sep=""), quote=F)
	}
	
	##Treatment up
	else if(control==F){
		##Filter for MeDIP treatment vs. MeDIP control
		frames=frames[as.numeric(frames[,value_ratio])<=down,]
		print(paste("Remaining number of frames where treatment/control ratio <= ", down, ": ", length(frames[,1]), sep=""), quote=F)
				
		##Filter for MeDIP treatment rpm signals higher then given quantile of the Input distribution
		if(input==T){
			print(paste("Estimated rpm threshold for input quantile ", quant, " is: ", t, sep=""), quote=F)
		}
		else{
			print(paste("Given rpm threshold is: ", t, sep=""), quote=F)
		}		
		frames=frames[as.numeric(frames[,value_rpm_B])>=t, ]
		print(paste("Remaining number of frames with treatment rpm >=", t, ": ", length(frames[,1]), sep=""), quote=F)		
		
		##Filter for rpm MeDIP treatment vs. Input
		if(input==T){
			frames=frames[((as.numeric(frames[,value_rpm_B])/as.numeric(frames[,value_input]))>=up), ]
			print(paste("Remaining number of frames with treatment vs. input ratio>=", up, ": ", length(frames[,1]), sep=""), quote=F)
		}
		print(paste("[Note: There are ", length(frames[,1][frames[,value_pvalue.ttest]==0 & frames[,value_pvalue.wilcox]==0]), " frames associated to a p-value==0.]", sep=""), quote=F)
		print(paste("[Note: There are ", length(frames[,1][frames[,value_ratio]==0]), " frames, where control/treatment ratio = 0 (i.e. control=0).]", sep=""), quote=F)
	}	
	
	return(as.data.frame(frames,stringsAsFactors=F))
}


