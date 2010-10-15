#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <Rmath.h>

SEXP coupling(SEXP temp_posPos,SEXP bin_size,SEXP temp_genomePos,SEXP temp_genomeCoup,SEXP noInfluenceFlankWin,SEXP fragmentLength,SEXP distanceVector){
    int score=0;
    int bin_size_c=INTEGER_VALUE(bin_size);
    int *ptrtemp_posPos,*ptrtemp_genomePos;
    double *ptrtemp_genomeCoup,*ptrdistanceVector;
    int noInfluenceFlankWin_c=INTEGER_VALUE(noInfluenceFlankWin);
    ptrtemp_posPos=INTEGER_POINTER(temp_posPos);
    ptrtemp_genomePos=INTEGER_POINTER(temp_genomePos);
    ptrdistanceVector=NUMERIC_POINTER(distanceVector);
    ptrtemp_genomeCoup=NUMERIC_POINTER(temp_genomeCoup);
    int genomePos_length=LENGTH(temp_genomePos);
    int fragmentLength_c=INTEGER_VALUE(fragmentLength);
    for (int j=0;j<LENGTH(temp_posPos);j++){
        int relative_posPos=floor(ptrtemp_posPos[j]/bin_size_c);
        int start = 1>(relative_posPos-noInfluenceFlankWin_c)?1:(relative_posPos-noInfluenceFlankWin_c);
        start=start< genomePos_length?start:genomePos_length;
        int stop=1> (relative_posPos+noInfluenceFlankWin_c+1)?1:(relative_posPos+noInfluenceFlankWin_c+1);
        stop=stop<genomePos_length?stop:genomePos_length;
        int k=start-1;
        for(k;k<stop;k++){
            int k_dist=fabs(ptrtemp_genomePos[k]-ptrtemp_posPos[j]);
            if(k_dist<(fragmentLength_c)){
                ptrtemp_genomeCoup[k]=ptrtemp_genomeCoup[k]+ptrdistanceVector[k_dist];
            }
            score=k;
        }
    }
    return(temp_genomeCoup);
}
