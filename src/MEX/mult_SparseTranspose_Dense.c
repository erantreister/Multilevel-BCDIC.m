#include "mex.h"
#include "bcdic.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input params:
    double *in_sparse, *in_dense, *in3;
    double *set_sparse;
    double *idxs_translation_inv;
    int set_sparse_size, rows_dense, cols_dense, cols_sparse; 
    
    // output params:
    double *result;
    
    
    in_sparse = (double *)mxGetPr(prhs[0]);  // Gi
    set_sparse = (double *)mxGetPr(prhs[1]); // all_relevant_columns
    set_sparse_size = mxGetM(prhs[1]);
    cols_sparse = (int)*((double*)mxGetPr(prhs[2]));
    
    in_dense = (double *)mxGetPr(prhs[3]);  // invAi
    cols_dense = mxGetN(prhs[3]); 
    rows_dense = mxGetM(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(cols_sparse, cols_dense, mxREAL);
    result = mxGetPr(plhs[0]);
    
	mult_SparseMatTranspose_DenseMat(result,in_sparse,set_sparse,set_sparse_size,cols_sparse,in_dense,cols_dense,rows_dense);
}
