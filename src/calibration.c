#include <R.h>
#include <stdio.h>
#include <stdlib.h>

	
void calibration(int *maxCoup, int *coupling, int *lc, int *signal, int *ls, double *mean_signal,int *mean_coupling){
    for(int i=0; i<= *maxCoup;i++){
        int sumc=0;int sums=0; int c=0; 
        for(int x=0;x<*ls;x++){
	    if(coupling[x]>=i & coupling[x]<i+1){sumc+=coupling[x];sums+=signal[x];c++;}
	}
	if(c!=0){
	  mean_signal[i]=(float)sums/(float)c;
	  mean_coupling[i]=(float)sumc/(float)c;
	}
	else{
	    mean_signal[i]=(float)*maxCoup+5;
	    mean_coupling[i]=(float)*maxCoup+5;
	}
    }
}
