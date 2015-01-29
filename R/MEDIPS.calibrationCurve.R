##########################################################################
##Function takes a MEDIPSset and COUPLINGSet object, and performs linear regression
##########################################################################
##Input:	MEDIPS SET and COUPLING SET
##Param:	MSet, CSet, input [T | F]
##Output:	List of calibration curve and linear regression results
##Modified:	1/09/2015
##Author:	Lukas Chavez

###############
##Function takes a MEDIPSset object, calculates the calibration curve and the parameters intercept and slope.
###############
MEDIPS.calibrationCurve <- function(MSet=NULL, CSet=NULL, input=FALSE){
	
	signal=		genome_count(MSet)
	coupling=	genome_CF(CSet)	
	maxCoup = floor(max(coupling)*0.8)
				
	cat("Calculating calibration curve...\n")			
	mean_signal = NULL
	coupling_level = NULL
	
	count_decrease_steps = 0
	max_signal_index = NULL

	first = TRUE
	n_coupling = 0

	for(i in 1:(maxCoup+1)){		

		#Test if coupling level is non-empty (can happen for ROIs)
		if(length(signal[coupling==i-1])>0){

			n_coupling = n_coupling + 1
			mean_signal = c(mean_signal, mean(signal[coupling==i-1], trim=0.1))
			coupling_level = c(coupling_level, i-1)
		
			##Test if mean_signal decreases
			if(!first){

				if(is.null(max_signal_index) & mean_signal[n_coupling]<mean_signal[n_coupling-1]){				
					count_decrease_steps=count_decrease_steps+1
					if(count_decrease_steps==3){
						max_signal_index=n_coupling-3		
					}
				}
				else if(is.null(max_signal_index) & mean_signal[n_coupling]>=mean_signal[n_coupling-1]){
					count_decrease_steps=0
				}	
			}
			else{
				first = FALSE
			}
		} 
		else{
			cat(paste("Skipping coupling level ", i, " (no data).\n", sep=""))
		}						
	}

	if(!input & is.null(max_signal_index)){stop("The dependency of coverage signals on sequence pattern (e.g. CpG) densities is different than expected. No linear model can be build, please check the calibration plot by providing the MSet object at ISet.")}
				
 	if(!input){
		##Linear curve by linear least squares regression
		##################
		cat("Performing linear regression...\n")
			
		fit=lm(mean_signal[1:max_signal_index] ~ coupling_level[1:max_signal_index])
		intercept=fit$coefficients[[1]]
		slope=fit$coefficients[[2]]
		cat(paste("Intercept:", intercept, "\n"))
		cat(paste("Slope:", slope, "\n"))
	}
	else{
		intercept=NA
		slope=NA
	}	
	gc()
	return(list(mean_signal=mean_signal, coupling_level=coupling_level, max_index=max_signal_index, intercept=intercept, slope=slope))
}
