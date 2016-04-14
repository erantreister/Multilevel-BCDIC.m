#include "mex.h"
#include "bcdic.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input params:
    double *in1, *in2, *in3;
    double *set1, *row_translation1;
    double *set2, *row_translation2;
	double *idxs_translation_inv;
    int setsize1, setsize2, colsize, rows, cols; 
    
    // output params:
    double *result;
    
    
    in1 = (double *)mxGetPr(prhs[0]);
    in2 = (double *)mxGetPr(prhs[1]);
    in3 = (double *)mxGetPr(prhs[2]);
    set1 = (double *)mxGetPr(prhs[3]);
    row_translation1 = (double *)mxGetPr(prhs[4]);
    set2 = (double *)mxGetPr(prhs[5]);
    row_translation2 = (double *)mxGetPr(prhs[6]);
    rows = (int)*mxGetPr(prhs[7]); // the number of rows in the result matrix is a parameter
	idxs_translation_inv = (double *)mxGetPr(prhs[8]); // idxs_in_relevant

    setsize1 = mxGetM(prhs[3]);
    setsize2 = mxGetM(prhs[5]);
    colsize = mxGetM(prhs[0]);
    cols = mxGetN(prhs[1]); // in2 is colsize by cols
    
    if (colsize != mxGetM(prhs[1])){
        printf("Error: Input matrices should have the number of rows");
        return;
    }

    
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    result = mxGetPr(plhs[0]);
    
    mult_sparse_AtB(result, in1, in2, colsize, set1, setsize1, row_translation1, NULL, rows, cols, 0);
    mult_sparse_AtB(result, in3, in2, colsize, set2, setsize2, row_translation2, idxs_translation_inv, rows, cols, 1);

}