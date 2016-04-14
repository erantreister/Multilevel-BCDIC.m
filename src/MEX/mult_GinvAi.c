#include "mex.h"
#include "bcdic.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input params:
    double *in1, *in2, *in3;
    double *set1, *row_translation1;
    double *idxs_translation_inv, *idxs_translation;
    int setsize1, colsize, rows, cols; 
    
    // output params:
    double *result;
    
    
    in1 = (double *)mxGetPr(prhs[0]);  // Gi (nnz's)
    in2 = (double *)mxGetPr(prhs[1]);  // invAi
    set1 = (double *)mxGetPr(prhs[2]); // all_relevant_columns
    idxs_translation_inv = (double *)mxGetPr(prhs[3]); // idxs_in_relevant
    idxs_translation = (double *)mxGetPr(prhs[4]);

    setsize1 = mxGetM(prhs[2]);
    colsize = mxGetM(prhs[1]);
    
    cols = mxGetN(prhs[1]); 
    rows = mxGetM(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    result = mxGetPr(plhs[0]);
    
    mult_AjWj(result, in1, in2, colsize, set1, setsize1, idxs_translation_inv, idxs_translation, rows, cols);
}
