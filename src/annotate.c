#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <Rmath.h>

SEXP annotate(SEXP ANNO, SEXP REGION){
    R_len_t ai, aj, ar, ac,ri, rj, rr, rc;
    ar = nrows(ANNO);ac = ncols(ANNO);
    rr = nrows(REGION);rc = ncols(REGION);
    SEXP ans;
    int *ID,n=0,*chr,*start,*stop;
    ID = (int *)malloc(n*sizeof(int));	
    chr = (int *)malloc(n*sizeof(int));
    start=(int *)malloc(n*sizeof(int));
    stop=(int *)malloc(n*sizeof(int));
    
    for (ri = 0; ri < rr; ++ri) {
         Rprintf("analysed: %i/%i\r",ri,rr);
         for (ai = 0; ai < ar; ++ai) {
	        if(REAL(REGION)[ri + rr * 1] < REAL(ANNO)[ai + ar * 1] & REAL(REGION)[ri + rr * 2] > REAL(ANNO)[ai + ar * 1]){
                    n++;
	            ID = (int*)realloc(ID, n * sizeof(int));
 		    chr = (int*)realloc(chr,n * sizeof(int) );
		    start = (int*)realloc(start,n * sizeof(int) );
		    stop = (int*)realloc(stop,n * sizeof(int) );
	            chr[n-1]=REAL(REGION)[ri + rr * 0];
                    ID[n-1]=REAL(ANNO)[ai + ar * 3];
	            start[n-1]=REAL(REGION)[ri + rr * 1];
		    stop[n-1]=REAL(REGION)[ri + rr * 2];
		}
 		else if (REAL(REGION)[ri + rr * 1] > REAL(ANNO)[ai + ar * 1] & REAL(REGION)[ri + rr * 2] < REAL(ANNO)[ai + ar * 2]){
                    n++;
	            ID = (int*)realloc(ID, n * sizeof(int));
 		    chr = (int*)realloc(chr,n * sizeof(int) );
		    start = (int*)realloc(start,n * sizeof(int) );
		    stop = (int*)realloc(stop,n * sizeof(int) );
	            chr[n-1]=REAL(REGION)[ri + rr * 0];
                    ID[n-1]=REAL(ANNO)[ai + ar * 3];
	            start[n-1]=REAL(REGION)[ri + rr * 1];
		    stop[n-1]=REAL(REGION)[ri + rr * 2];
		}
 		else if (REAL(REGION)[ri + rr * 1] < REAL(ANNO)[ai + ar * 2] & REAL(REGION)[ri + rr * 2] > REAL(ANNO)[ai + ar * 2]){
                    n++;
	            ID = (int*)realloc(ID, n * sizeof(int));
 		    chr = (int*)realloc(chr,n * sizeof(int) );
		    start = (int*)realloc(start,n * sizeof(int) );
		    stop = (int*)realloc(stop,n * sizeof(int) );
	            chr[n-1]=REAL(REGION)[ri + rr * 0];
                    ID[n-1]=REAL(ANNO)[ai + ar * 3];
	            start[n-1]=REAL(REGION)[ri + rr * 1];
		    stop[n-1]=REAL(REGION)[ri + rr * 2];
		} 
	 }
    }
    Rprintf("analysed: %i/%i found: %i\n",ri,rr,n);
    PROTECT(ans = allocMatrix(REALSXP, n, 4)); //nrow,ncol
    for(int x=0;x<n;x++){
      REAL(ans)[x+n*0]=chr[x];
      REAL(ans)[x+n*1]=start[x];
      REAL(ans)[x+n*2]=stop[x];
      REAL(ans)[x+n*3]=ID[x];
    }
    
    UNPROTECT(1);
    return(ans);
}
