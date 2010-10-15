#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <Rmath.h>

SEXP calculate (SEXP Input,SEXP select,SEXP matrix,R_len_t ii,SEXP MEDIPS,SEXP MEDIPS2,int start_pos, int stop_pos,SEXP rho,SEXP fn,SEXP var_env,SEXP vari,SEXP math_env,SEXP math,SEXP ttest,SEXP t_env){
    
    if(!isFunction(math)){error("math should be function");} 
    R_len_t nrow;
    nrow = nrows(matrix);
    SEXP R_fcall,Val,R_fcall_var,R_fcall_math,R_fcall_ttest,X,nX,CF;
   
    int items=(stop_pos-start_pos)+1;
    REAL(matrix)[ii + nrow * 3] = items;
     
    int protected=0;
    PROTECT(X=allocVector(REALSXP,items)); protected++;
    PROTECT(nX=allocVector(REALSXP,items)); protected++;
    PROTECT(CF=allocVector(REALSXP,items)); protected++;
    
    
    double *ptrgenome_raw,*ptrgenome_pos,*ptrgenome_norm,*ptrgenome_CF;
    ptrgenome_raw=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_raw")));
    ptrgenome_norm=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_norm")));
    ptrgenome_CF=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_CF")));
    ptrgenome_pos=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_pos")));
    
    SEXP IN;double *ptrinput;
    if(Input != R_NilValue){
        PROTECT(IN=allocVector(REALSXP,items)); protected++;
 	ptrinput=NUMERIC_POINTER(Input);
     }
    
    int x=0;
    int cf_start=start_pos;
    while(cf_start<=stop_pos){
 	REAL(CF)[x]=ptrgenome_CF[cf_start];
        if(Input != R_NilValue) REAL(IN)[x]=ptrinput[cf_start];
	cf_start++;
        x++;
    }
    
    PROTECT(R_fcall_math=lang2(math,CF));
    double meanCF =asReal(eval(R_fcall_math, math_env));
    REAL(matrix)[ii+nrow*4]=meanCF;
    UNPROTECT(1);
    
    if (Input != R_NilValue) {
	PROTECT(R_fcall_math=lang2(math,IN));
	double meanIN =asReal(eval(R_fcall_math, math_env));
 	REAL(matrix)[ii+nrow*5]=meanIN;
 	UNPROTECT(1);	
    }
    else REAL(matrix)[ii+nrow*5]=NA_REAL;
		    
    x=0;
    int start_pos_x=start_pos;
    int uniX=0;double diffX=ptrgenome_raw[start_pos_x];int NullX=0;
    int uninX=0;double diffnX=ptrgenome_norm[start_pos_x];int NullnX=0;
    while(start_pos_x<=stop_pos){
        if(diffX==ptrgenome_raw[start_pos_x]){uniX++;}
	if(0==ptrgenome_raw[start_pos_x]){NullX++;}
	if(diffnX==ptrgenome_norm[start_pos_x]){uninX++;}
	if(0==ptrgenome_norm[start_pos_x]){NullnX++;}
        REAL(X)[x]=ptrgenome_raw[start_pos_x];
        REAL(nX)[x]=ptrgenome_norm[start_pos_x];
        REAL(CF)[x]=ptrgenome_CF[start_pos_x];
 	if(Input != R_NilValue) REAL(IN)[x]=ptrinput[start_pos_x];
	start_pos_x++;
        x++;
    }
     
    if(!isS4(MEDIPS2)){
        REAL(matrix)[ii+nrow*7]=NA_REAL;
        REAL(matrix)[ii+nrow*9]=NA_REAL;
	REAL(matrix)[ii+nrow*11]=NA_REAL;
	REAL(matrix)[ii+nrow*13]=NA_REAL;
	REAL(matrix)[ii+nrow*15]=NA_REAL;
	REAL(matrix)[ii+nrow*16]=NA_REAL;
	REAL(matrix)[ii+nrow*17]=NA_REAL;
	REAL(matrix)[ii+nrow*18]=NA_REAL;
    }
	   
    PROTECT(R_fcall_math=lang2(math,X));
    double meanX =asReal(eval(R_fcall_math, math_env));
    REAL(matrix)[ii+nrow*6]=meanX;
    UNPROTECT(1);
	    
    PROTECT(R_fcall_math=lang2(math,nX));
    double meannX =asReal(eval(R_fcall_math, math_env));
    REAL(matrix)[ii+nrow*8]=meannX;
    UNPROTECT(1);
	
    //mean ams1
    if(meanCF!=0) REAL(matrix)[ii+nrow*10]=meannX/meanCF;
    else REAL(matrix)[ii+nrow*10]=-1;
    
	
    if(NUMERIC_VALUE(select)==1){
       //varX
        PROTECT(R_fcall_var=lang2(vari,X));
        double varX =asReal(eval(R_fcall_var, var_env));
        REAL(matrix)[ii+nrow*12]=varX;
        UNPROTECT(1);
        if(meanX!=0){
            //varcX
            double varcX =(sqrt(varX))/meanX;
	    REAL(matrix)[ii+nrow*14]=varcX;
	}
	else REAL(matrix)[ii+nrow*14]=-1;
    }
    else{
        PROTECT(R_fcall_var=lang2(vari,nX));
	double varnX =asReal(eval(R_fcall_var, var_env));
	REAL(matrix)[ii+nrow*12]=varnX;
 	UNPROTECT(1);
	    
	if(meannX!=0){
	    double varcXn =(sqrt(varnX))/meannX;
	    REAL(matrix)[ii+nrow*14]=varcXn;
	}   
	else REAL(matrix)[ii+nrow*14]=-1;
    }
 	
    
    if(isS4(MEDIPS2)){
        SEXP Y,nY; double *ptrgenome_raw2,*ptrgenome_norm2;
        ptrgenome_raw2=NUMERIC_POINTER(GET_SLOT(MEDIPS2,install("genome_raw")));
  	ptrgenome_norm2=NUMERIC_POINTER(GET_SLOT(MEDIPS2,install("genome_norm")));
	PROTECT(Y=allocVector(REALSXP,items)); protected++;
  	PROTECT(nY=allocVector(REALSXP,items)); protected++;
  	x=0;
	int uniY=0;double diffY=ptrgenome_raw2[start_pos];int NullY=0;
	int uninY=0;double diffnY=ptrgenome_norm2[start_pos];int NullnY=0;
	while(start_pos<=stop_pos){
	    if(diffY==ptrgenome_raw2[start_pos]){uniY++;}
	    if(0==ptrgenome_raw2[start_pos]){NullY++;}
	    
	    if(diffnY==ptrgenome_norm2[start_pos]){uninY++;}
	    if(0==ptrgenome_norm2[start_pos]){NullnY++;}
            
	    REAL(Y)[x]=ptrgenome_raw2[start_pos];
            REAL(nY)[x]=ptrgenome_norm2[start_pos];
            start_pos++;
            x++;
 	}
	//meanY
 	PROTECT(R_fcall_math=lang2(math,Y));
 	double meanY =asReal(eval(R_fcall_math, math_env));
 	REAL(matrix)[ii+nrow*7]=meanY;
 	UNPROTECT(1);
		    
	//mean nY
 	PROTECT(R_fcall_math=lang2(math,nY));
 	double meannY =asReal(eval(R_fcall_math, math_env));
 	REAL(matrix)[ii+nrow*9]=meannY;
  	UNPROTECT(1);
   
	//mean SW2
  	if(meanCF!=0) REAL(matrix)[ii+nrow*11]=meannY/meanCF;
        else REAL(matrix)[ii+nrow*11]=-1;
	
	if(NUMERIC_VALUE(select)==1){
 	    PROTECT(R_fcall_var=lang2(vari,Y));
 	    double varY =asReal(eval(R_fcall_var, var_env));
	    REAL(matrix)[ii+nrow*13]=varY;
  	    UNPROTECT(1);
	
 	    if(meanY!=0){
                double ratio=meanX/meanY;
		REAL(matrix)[ii+nrow*16]=ratio;
                
		double varcY =(sqrt(varY))/meanY;
 	        REAL(matrix)[ii+nrow*15]=varcY;
 	    }
 	    else { 
	        REAL(matrix)[ii+nrow*15]=-1;
	        REAL(matrix)[ii+nrow*16]=-1;
	    }
	    if(uniX==items & uniY==items & NullX==items & NullY==items & items<=4){
	        REAL(matrix)[ii+nrow*17]=0;
	        REAL(matrix)[ii+nrow*18]=0;
	    }
	    else{
 	        if(((uniX==items & uniY!=items)  | (uniY==items & uniX!=items) | (uniY!=items & uniX!=items))&items>=5){
 	            PROTECT(R_fcall=lang3(fn,X,Y));
  	            double p=asReal( VECTOR_ELT(eval(R_fcall, rho),2));
  	            REAL(matrix)[ii+nrow*17]=p;
   	            UNPROTECT(1);
		    PROTECT(R_fcall_ttest=lang3(ttest,X,Y));
                    double tp =asReal( VECTOR_ELT(eval(R_fcall_ttest, t_env),2));
                    REAL(matrix)[ii+nrow*18]=tp;
                    UNPROTECT(1); 
		 }
		 else{
		    REAL(matrix)[ii+nrow*17]=0;
    		    REAL(matrix)[ii+nrow*18]=0;
		 }
 	    } 
 	}
 	else{
            PROTECT(R_fcall_var=lang2(vari,nY));
 	    double varnY =asReal(eval(R_fcall_var, var_env));
 	    REAL(matrix)[ii+nrow*13]=varnY;
  	    UNPROTECT(1);
	    
 	    if(meannY!=0){
	        double ratio=meannX/meannY;
 		REAL(matrix)[ii+nrow*16]=ratio;
	        
		//varcX
 	        double varcYn =(sqrt(varnY))/meannY;
 	       	REAL(matrix)[ii+nrow*15]=varcYn; 	        
	    }   
 	    else{
	        REAL(matrix)[ii+nrow*15]=-1;
	        REAL(matrix)[ii+nrow*16]=-1;
	    }
	    
 	    if(uninX==items & uninY==items & NullnX==items & NullnY==items & items<=4){
		REAL(matrix)[ii+nrow*17]=0;
	        REAL(matrix)[ii+nrow*18]=0;
	    }
	    else{
 	        if(((uniX==items & uniY!=items)  | (uniY==items & uniX!=items) | (uniY!=items & uniX!=items))&items>=5){
 	            PROTECT(R_fcall=lang3(fn,X,Y));
  	            double p=asReal( VECTOR_ELT(eval(R_fcall, rho),2));
  	            REAL(matrix)[ii+nrow*17]=p;
   	            UNPROTECT(1);
		    PROTECT(R_fcall_ttest=lang3(ttest,X,Y));
                    double tp =asReal( VECTOR_ELT(eval(R_fcall_ttest, t_env),2));
                    REAL(matrix)[ii+nrow*18]=tp;
                    UNPROTECT(1); 
		 }
		 else{
		    REAL(matrix)[ii+nrow*17]=0;
    		    REAL(matrix)[ii+nrow*18]=0;
		 }
 	    }
	}	
    }
    UNPROTECT(protected);
    return(R_NilValue);
}

SEXP roiprofile (SEXP Input, SEXP select,SEXP ROI,SEXP bin_pos, SEXP MEDIPS,SEXP MEDIPS2,SEXP rho,SEXP fn,SEXP var_env, SEXP vari,SEXP math_env,SEXP math,SEXP ttest,SEXP t_env,SEXP factor){
    int bin_size_c,start,stop,start_pos,stop_pos,*ptrbin_pos,nProtected = 0;
    double *ptrgenome_pos,chr,*ptrfactor;
    bin_size_c=INTEGER_VALUE(GET_SLOT(MEDIPS,install("bin_size")));
    ptrbin_pos=INTEGER_POINTER(bin_pos);
    ptrfactor=NUMERIC_POINTER(factor);
    ptrgenome_pos=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_pos")));
    R_len_t ii, jj, nrow, ncol;
    SEXP matrix;
    nrow = nrows(ROI);ncol = ncols(ROI);
    PROTECT(matrix = allocMatrix(REALSXP, nrow, 19));

    for (ii = 0; ii < nrow; ++ii) {
        int wchr=0;
        chr = REAL(ROI)[ii + nrow * 0] ;
	while(wchr< LENGTH(factor)){
	    if(chr==ptrfactor[wchr]){
		break;
	    }
	    wchr++;
	}
        Rprintf("Analysed %i / %i \r",ii,nrow);

	start = REAL(ROI)[ii + nrow * 1] ;
	stop= REAL(ROI)[ii + nrow * 2] ;	
	REAL(matrix)[ii + nrow * 0] = chr ;
        REAL(matrix)[ii + nrow * 1] = start;	
        REAL(matrix)[ii + nrow * 2] = stop;
	if(wchr == 0){
	    start_pos=ceil(start/bin_size_c);
	    stop_pos=ceil(stop/bin_size_c);
	    if(start != ptrgenome_pos[start_pos]) start_pos=start_pos+1;
	    else start_pos=start_pos;
 	    if(stop != ptrgenome_pos[stop_pos] & (stop_pos-start_pos)>1) stop_pos=stop_pos-1;
 	    else stop_pos=stop_pos;
            if (start < stop & !ISNA(REAL(ROI)[ii + nrow * 0])) {
   		calculate(Input,select,matrix,ii,MEDIPS,MEDIPS2,start_pos,stop_pos,rho,fn, var_env, vari,math_env,math,ttest,t_env);
            }
	    else{
	        for (jj = 3; jj <= 18; ++jj) {REAL(matrix)[ii + nrow * jj]=NA_REAL;}
	    }
	}
	else{
	    start_pos=ptrbin_pos[wchr-1]+(int)ceil(start/bin_size_c);
	    stop_pos=ptrbin_pos[wchr-1]+(int)ceil(stop/bin_size_c);
	    if(start != ptrgenome_pos[start_pos]) start_pos=start_pos+1;
	    else start_pos=start_pos;
	    if(stop != ptrgenome_pos[stop_pos] & (stop_pos-start_pos)>1) stop_pos=stop_pos-1;
 	    else stop_pos=stop_pos;
	    if (start < stop & !ISNA(REAL(ROI)[ii + nrow * 0])) {
	        calculate(Input,select,matrix,ii,MEDIPS,MEDIPS2,start_pos,stop_pos,rho,fn, var_env, vari,math_env,math,ttest,t_env);
            }
	    else{
	        for (jj = 3; jj <= 18; ++jj) {REAL(matrix)[ii + nrow * jj]=NA_REAL;}}
	}
    }
    Rprintf("Analysed %i / %i \r",ii,nrow);
    Rprintf("\n");
    UNPROTECT(1);
    return(matrix);
}
