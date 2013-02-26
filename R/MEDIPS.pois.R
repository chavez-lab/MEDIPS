##########################################################################
##Function estimates methylation level based on a series of coupling factor dependent poisson distributions.
##########################################################################
##Input:	MEDIPS SET, CS
##Output:	vector of methylation probabilities
##Modified:	11/17/2011
##Author:	Lukas Chavez

MEDIPS.pois <- function(MSet=NULL, CSet=NULL, ccObj=NULL){
	
	counts = genome_count(MSet)
	coupling = genome_CF(CSet)
	#win_size = window_size(MSet)
		
	if(is.null(ccObj)){
		##Calculate calibration curve
		#####################################
		ccObj = MEDIPS.calibrationCurve(MSet=MSet, CSet=CSet, input=F)
	}
	
	mean_count = ccObj$mean_signal
	coupling_level = ccObj$coupling_level
	intercept = ccObj$intercept
	slope = ccObj$slope
	max_index =  ccObj$max_index
	rm(ccObj)
	gc()
	
	##Estimate methylation w.r.t series of poisson ditributions
	############################################################
	cat("Estimate methylation...\n")
	count.pois = rep(1, length(counts))
	count.pois[coupling==0] = 0
	for(i in 2:(max(coupling)+1)){		
		
		##Mean given by calibration curve
		if(i<=max_index){
			count.pois[coupling==coupling_level[i]] = count.pois[coupling==coupling_level[i]] - ppois(q=counts[coupling==coupling_level[i]], lambda=mean_count[i], lower.tail=F)
		}
		##Mean estimated by linear regression obtained parameters
		else{
			est.mean = intercept + (slope * (i-1))
			count.pois[coupling==(i-1)] = count.pois[coupling==(i-1)] - ppois(q=counts[coupling==(i-1)], lambda=est.mean, lower.tail=F)
		}
	}
	
	gc()
	return(count.pois)	
}
