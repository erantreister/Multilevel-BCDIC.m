#include "mex.h"
#include "bcdic.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input params:
    double *values, *Grad, *Prec, *FreeSet;
    double epsilon;

    int rows, cols,elems; 
    
    // output params:
    double *result;
   
	double one_lambda = mxGetScalar(prhs[0]);
    values  = (double *)mxGetPr(prhs[1]);
	rows = (int)mxGetM(prhs[1]);
    cols = (int)mxGetN(prhs[1]);

    epsilon = (double)*mxGetPr(prhs[2]);
	Grad = (double *)mxGetPr(prhs[3]);
	Prec = (double *)mxGetPr(prhs[4]);
    FreeSet = (double*)mxGetPr(prhs[5]);
	elems = max((int)mxGetM(prhs[5]),(int)mxGetN(prhs[5]));
	
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    result = mxGetPr(plhs[0]);
	
	soft_shrinkage_Free(result, one_lambda, Prec, values, Grad, epsilon, FreeSet, elems);
}
