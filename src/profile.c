#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <Rmath.h>

SEXP profile(SEXP select,SEXP input,SEXP step,SEXP frame, SEXP MEDIPS,SEXP MEDIPS2, SEXP rho,SEXP fn,SEXP var_env, SEXP vari,SEXP math_env,SEXP math,SEXP ttest,SEXP t_env){
 
    if(!isFunction(math)){error("math should be function");}
       
    SEXP R_fcall,Val,pvalue,new_chr,X,Y,nX,nY,matrix,CF,IN;
    SEXP R_fcall_var,R_fcall_math,R_fcall_ttest;
    double *ptrgenome,step_c,frame_c,bin_size_c,*ptrgenome_raw,*ptrgenome_norm,*ptrgenome_CF, *ptr_input,*ptrchr_l;
    step_c=INTEGER_VALUE(step);
    frame_c=INTEGER_VALUE(frame);
    if(step_c==NA_INTEGER){step_c=frame_c;}
    bin_size_c=INTEGER_VALUE(GET_SLOT(MEDIPS,install("bin_size")));
    ptrchr_l=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("chr_lengths")));
    ptrgenome=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_pos")));
    ptrgenome_raw=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_raw")));
    ptrgenome_norm=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_norm")));
    ptrgenome_CF=NUMERIC_POINTER(GET_SLOT(MEDIPS,install("genome_CF")));
 
    
    double *ptrgenome_raw2,*ptrgenome_norm2;
    if(isS4(MEDIPS2)){
        ptrgenome_raw2=NUMERIC_POINTER(GET_SLOT(MEDIPS2,install("genome_raw")));
        ptrgenome_norm2=NUMERIC_POINTER(GET_SLOT(MEDIPS2,install("genome_norm")));
    }
    if(input != R_NilValue){ptr_input=NUMERIC_POINTER(input);}
    
    int l=LENGTH(GET_SLOT(MEDIPS,install("chr_lengths")));
    int position=0;
    int items=0; 
    char *chromosome[l];
    int nrow=0;
    int row_index=0;
    int ncol=19;
    for(int i=0;i<l;i++){
       double step_items = step_c/bin_size_c;
       double stop_index =ceil((ptrchr_l[i]/bin_size_c)/step_items);    
       double wie=(ptrchr_l[i]/bin_size_c)/step_items;
       nrow+=(int)stop_index;
    }
    PROTECT(new_chr=allocVector(STRSXP,nrow));
    PROTECT(matrix = allocMatrix(REALSXP, nrow, ncol));
    
    for(int i=0;i<l;i++){
        chromosome[i] = R_alloc(strlen(CHAR(STRING_ELT(GET_SLOT(MEDIPS,install("chr_names")), i))), sizeof(char)); 
  	strcpy(chromosome[i], CHAR(STRING_ELT(GET_SLOT(MEDIPS,install("chr_names")), i))); 
	int items = frame_c/bin_size_c;
	int start=1;
	int step_items = step_c/bin_size_c;
	int overlap_item=items-step_items;
 	int stop_index =trunc((ptrchr_l[i]/bin_size_c)/step_items);
 	int t=0;
	int add;
	
	while(t<=stop_index){	
	   int prot=0;
	   int stop=start+(items*bin_size_c)-1;
	   int x =0;
	   if(start>=ptrchr_l[i]){break;}
	   if(stop>=ptrchr_l[i]){
	       items=(ptrchr_l[i]-start)/bin_size_c+1;
	       stop=ptrchr_l[i];
	    }
	    row_index++;
            Rprintf("%s start: %i stop: %i (chromosome length: %i)     \r ",chromosome[i],start,stop,(int)ptrchr_l[i]);
  	    
	    REAL(matrix)[row_index-1+nrow*0]=i+1;
	    REAL(matrix)[row_index-1+nrow*1]=start;
	    REAL(matrix)[row_index-1+nrow*2]=stop;
	    REAL(matrix)[row_index-1+nrow*3]=items;
	    PROTECT(X=allocVector(REALSXP,items));
	    PROTECT(nX=allocVector(REALSXP,items));
  	    PROTECT(CF=allocVector(REALSXP,items));
	    PROTECT(IN=allocVector(REALSXP,items));
	    prot+=4;
	    
	    
	    if(isS4(MEDIPS2)){
		PROTECT(Y=allocVector(REALSXP,items));
		PROTECT(nY=allocVector(REALSXP,items));
		prot+=2;
	    }
	    else{
	        REAL(matrix)[row_index-1+nrow*7]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*9]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*11]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*13]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*15]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*16]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*17]=NA_REAL;
		REAL(matrix)[row_index-1+nrow*18]=NA_REAL;
	    }
		
	    int uniX=0;double diffX=ptrgenome_raw[position];int NullX=0;
	    int uniY=0;double diffY= 0; int NullY=0;
	    int uninX=0;double diffnX=ptrgenome_norm[position];int NullnX=0;
	    int uniYn=0;double diffnY=0;int NullnY=0;
	    if(isS4(MEDIPS2)){diffY= ptrgenome_raw2[position];
	        diffnY=ptrgenome_norm2[position];
	    }
	    
	    while( x<items){
	       if(diffX==ptrgenome_raw[position]){uniX++;}
               if(0==ptrgenome_raw[position]){NullX++;}
	       
	       if(diffnX==ptrgenome_norm[position]){uninX++;}
               if(0==ptrgenome_norm[position]){NullnX++;}
	       
	       REAL(X)[x]=ptrgenome_raw[position];
	       REAL(nX)[x]=ptrgenome_norm[position];
	       REAL(CF)[x]=ptrgenome_CF[position];
	          
	       if(isS4(MEDIPS2)){
	           if(diffY==ptrgenome_raw2[position]){uniY++;}
	           if(0==ptrgenome_raw2[position]){NullY++;}
	           
		   if(diffnY==ptrgenome_norm2[position]){uniYn++;}
	           if(0==ptrgenome_norm2[position]){NullnY++;}
		   
		   REAL(Y)[x]=ptrgenome_raw2[position];
		   REAL(nY)[x]=ptrgenome_norm2[position];
	       }
	       if (input != R_NilValue) {REAL(IN)[x]=ptr_input[position];}
	       x++;
  	       position+=1;

	    }
	    
            //mean CF
	    PROTECT(R_fcall_math=lang2(math,CF));
	    double meanCF =asReal(eval(R_fcall_math, math_env));
	    REAL(matrix)[row_index-1+nrow*4]=meanCF;
 	    UNPROTECT(1);
	    
	    if (input != R_NilValue) {
			PROTECT(R_fcall_math=lang2(math,IN));
			double meanIN =asReal(eval(R_fcall_math, math_env));
			REAL(matrix)[row_index-1+nrow*5]=meanIN;
			UNPROTECT(1);	
	    }
	    else{REAL(matrix)[row_index-1+nrow*5]=NA_REAL;}
	    //meanX
	    PROTECT(R_fcall_math=lang2(math,X));
	    double meanX =asReal(eval(R_fcall_math, math_env));
	    REAL(matrix)[row_index-1+nrow*6]=meanX;
 	    UNPROTECT(1);	    
	   
	    //mean nX
	    PROTECT(R_fcall_math=lang2(math,nX));
	    double meannX =asReal(eval(R_fcall_math, math_env));
 	    REAL(matrix)[row_index-1+nrow*8]=meannX;
	    UNPROTECT(1);
	    
	    //mean ams
	    if(meanCF!=0) REAL(matrix)[row_index-1+nrow*10]=meannX/meanCF;
            else REAL(matrix)[row_index-1+nrow*10]=-1;
	    
	    
	    if(NUMERIC_VALUE(select)==1){
                //varX
                PROTECT(R_fcall_var=lang2(vari,X));
                double varX =asReal(eval(R_fcall_var, var_env));
                REAL(matrix)[row_index-1+nrow*12]=varX;
                UNPROTECT(1);
                if(meanX!=0){
		    //varcX
                    double varcX =(sqrt(varX))/meanX;
	            REAL(matrix)[row_index-1+nrow*14]=varcX;
	        }
    	        else REAL(matrix)[row_index-1+nrow*14]=-1;
            }
            else{
                PROTECT(R_fcall_var=lang2(vari,nX));
	        double varnX =asReal(eval(R_fcall_var, var_env));
	        REAL(matrix)[row_index-1+nrow*12]=varnX;
 	        UNPROTECT(1);
                if(meannX!=0){
	            double varcXn =(sqrt(varnX))/meannX;
	            REAL(matrix)[row_index-1+nrow*14]=varcXn;
	        }   
	        else REAL(matrix)[row_index-1+nrow*14]=-1;
            }

	    if(isS4(MEDIPS2)){
		//meanY
       	        PROTECT(R_fcall_math=lang2(math,Y));
	        double meanY =asReal(eval(R_fcall_math, math_env));
	        REAL(matrix)[row_index-1+nrow*7]=meanY;
	        UNPROTECT(1);
	        
	        //mean nY
    	        PROTECT(R_fcall_math=lang2(math,nY));
	        double meannY =asReal(eval(R_fcall_math, math_env));
		REAL(matrix)[row_index-1+nrow*9]=meannY;
 		UNPROTECT(1);
		
		//mean ams
	        if(meanCF!=0) REAL(matrix)[row_index-1+nrow*11]=meannY/meanCF;
                else REAL(matrix)[row_index-1+nrow*11]=-1;
		 
	        
	        if(NUMERIC_VALUE(select)==1){
  	            PROTECT(R_fcall_var=lang2(vari,Y));
  	            double varY =asReal(eval(R_fcall_var, var_env));
 	            REAL(matrix)[row_index-1+nrow*13]=varY;
   	            UNPROTECT(1);
 		    if(meanY!=0){
                        double ratio=meanX/meanY;
 		        REAL(matrix)[row_index-1+nrow*16]=ratio;
                        double varcY =(sqrt(varY))/meanY;
  	                REAL(matrix)[row_index-1+nrow*15]=varcY;
  	            }
		    else {
		        REAL(matrix)[row_index-1+nrow*15]=-1;
			REAL(matrix)[row_index-1+nrow*16]=-1;
			}
	            if(uniX==items & uniY==items & NullX==items & NullY==items & items<=4){
                        REAL(matrix)[row_index-1+nrow*17]=0;
                        REAL(matrix)[row_index-1+nrow*18]=0;
   	            }
 	            else{
		        
  	                if(((uniX==items & uniY!=items)  | (uniY==items & uniX!=items) | (uniY!=items & uniX!=items))&items>=5){
 			    PROTECT(R_fcall=lang3(fn,X,Y));
  	                    double p=asReal( VECTOR_ELT(eval(R_fcall, rho),2));
  	                    REAL(matrix)[row_index-1+nrow*17]=p;
   	                    UNPROTECT(1);
			    PROTECT(R_fcall_ttest=lang3(ttest,X,Y));
                            double tp =asReal( VECTOR_ELT(eval(R_fcall_ttest, t_env),2));
                            REAL(matrix)[row_index-1+nrow*18]=tp;
                            UNPROTECT(1); 
			}
			else{
			    REAL(matrix)[row_index-1+nrow*17]=0;
    			    REAL(matrix)[row_index-1+nrow*18]=0;
			}
			
 	            } 
 	        }
	    	
 	        else{
                    PROTECT(R_fcall_var=lang2(vari,nY));
 	            double varnY =asReal(eval(R_fcall_var, var_env));
 	            REAL(matrix)[row_index-1+nrow*13]=varnY;
  	            UNPROTECT(1);
	        
			    
 		    if(meannY!=0){
 	                 double ratio=meannX/meannY;
 		         REAL(matrix)[row_index-1+nrow*16]=ratio;
	                 //varcX
  	                 double varcYn =(sqrt(varnY))/meannY;
 	                 REAL(matrix)[row_index-1+nrow*15]=varcYn;
  	            }   
  	            else {
 		        REAL(matrix)[row_index-1+nrow*15]=-1;
 			REAL(matrix)[row_index-1+nrow*16]=-1;
 		    }
   	            if(uninX==items & uniYn==items & NullnX==items & NullnY==items & items<=4){
   		        REAL(matrix)[row_index-1+nrow*17]=0;
  		        REAL(matrix)[row_index-1+nrow*18]=0;
 	            }
		
		    else{
		        if(((uniX==items & uniY!=items)  | (uniY==items & uniX!=items) | (uniY!=items & uniX!=items))&items>=5){
 			    PROTECT(R_fcall=lang3(fn,X,Y));
  	                    double p=asReal( VECTOR_ELT(eval(R_fcall, rho),2));
  	                    REAL(matrix)[row_index-1+nrow*17]=p;
   	                    UNPROTECT(1);
			    PROTECT(R_fcall_ttest=lang3(ttest,X,Y));
                            double tp =asReal( VECTOR_ELT(eval(R_fcall_ttest, t_env),2));
                            REAL(matrix)[row_index-1+nrow*18]=tp;
                            UNPROTECT(1); 
			}
			else{
			    REAL(matrix)[row_index-1+nrow*17]=0;
    			    REAL(matrix)[row_index-1+nrow*18]=0;
			}
 		    }
	        }
	    
	    }	
	    
	    UNPROTECT(prot);
	    position+=1;
  	    start+=step_c;
 	    if(t!=stop_index){
		add=(frame_c/bin_size_c)-items;
		position-=overlap_item+1-add;
 	    }
  	    else{
		position-=1;
	    }
 	    t++;
	}
    Rprintf("\n");
    }
    UNPROTECT(2);
//     return(R_NilValue);
    return(matrix);
}

