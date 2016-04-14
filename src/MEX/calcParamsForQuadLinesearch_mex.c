#include "mex.h"
#include "bcdic.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input params:
    double *ASi, *Delta_i, *M_Delta_i, *Gi_before_shrink, 
            *Free_Neighbors, *Free_idxs;
    int size_Neighbors, size_idxs; 
    
    // output params:
    double *ans_a1,*ans_a2,*ans_ASi,*ans_Delta_i;
    
    
    ASi                 = (double *)mxGetPr(prhs[0]);
    Delta_i             = (double *)mxGetPr(prhs[1]);
    M_Delta_i           = (double *)mxGetPr(prhs[2]);
    Gi_before_shrink    = (double *)mxGetPr(prhs[3]);
    Free_Neighbors      = (double *)mxGetPr(prhs[4]);
    Free_idxs           = (double *)mxGetPr(prhs[5]);
    size_Neighbors      = (int)max(mxGetM(prhs[4]),mxGetM(prhs[4]));
    size_idxs           = (int)max(mxGetM(prhs[5]),mxGetM(prhs[5]));
    

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(size_Neighbors + size_idxs, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(size_Neighbors + size_idxs, 1, mxREAL);
    ans_a2      = mxGetPr(plhs[0]);
    ans_a1      = mxGetPr(plhs[1]);    
    ans_ASi     = mxGetPr(plhs[2]);
    ans_Delta_i = mxGetPr(plhs[3]);
    
    calcParamsForQuadLinesearch(ASi, Delta_i, M_Delta_i, Gi_before_shrink, 
            Free_Neighbors, Free_idxs, size_Neighbors, size_idxs,
            ans_a1, ans_a2, ans_ASi, ans_Delta_i);
}
