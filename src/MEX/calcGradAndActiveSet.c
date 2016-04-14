
#if !defined(_WIN32)
#define ddot ddot_
#endif

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"
#include "omp.h"

#define BLOCK_SIZE 8192

void calcInnerXXT(double* X,double* Y, double* temp_block,size_t n_rows, size_t n_cols, size_t m, size_t offset){
	char *trA = "T";
	char *trB = "N";
	double one = 1.0, zero = 0.0;
	dgemm(trA, trB, &n_rows, &n_cols, &m, &one, X+offset*m , &m, Y, &m, &zero, temp_block, &n_rows);
}

size_t calcGradAndActive(double* X,double* Y, size_t n, size_t m, double* invAi,
	mwIndex *Ai_rows, mwIndex *Ai_starts, double* Ai_vals, double lambda, double* idxs ,size_t num_idxs, size_t NNZ, 
	unsigned** I, double** GFi, double** XXt, int type){
    size_t nnz = 0;
	size_t num_blocks = n/BLOCK_SIZE + (n%BLOCK_SIZE > 0);
	
	size_t k_block;
	double* block = malloc(BLOCK_SIZE*num_idxs*sizeof(double));
	
	
	for (k_block = 0 ; k_block < num_blocks ; k_block++){
		size_t num_rows = BLOCK_SIZE;
		int j;
		size_t offset = 0;
		size_t i_block = k_block % num_blocks;
		offset = i_block*BLOCK_SIZE;
		if (i_block==num_blocks - 1){
			num_rows = n - offset;
		}
		calcInnerXXT(X,Y, block, num_rows, num_idxs, m, offset);
		for (j=0 ; j < num_idxs ; j++){
			size_t k_inner,i,i_start = 0, i_A;
			if (type==0){
				if ((offset+num_rows)<((size_t)idxs[j])-1){
					continue;
				}
				if (offset < (size_t)idxs[j]-1){
					i_start = (size_t)idxs[j] - 1 - offset;
				}else{
					i_start = 0;
				}
			}
			for (i = i_start ; i < num_rows ; i++){
				double XX_val = 0;
				double A_val = 0;
				double GFi_val;
				int flag_in_active = 0;
				k_inner = i + offset + j*n;
				XX_val = block[i + j*num_rows];
				GFi_val = XX_val - invAi[k_inner];
				if (fabs(GFi_val) > lambda){
					flag_in_active = 1;
				}else{
					for (i_A = Ai_starts[j]; i_A < Ai_starts[j+1] ; i_A++){
						if (Ai_rows[i_A] == i + offset){
							flag_in_active = 1;
							break;
						}
					}
				}
				if (flag_in_active){
					size_t nnz_local;
					{
						if (nnz == NNZ){
							unsigned* big_I = malloc(2*NNZ*sizeof(unsigned));
							double* big_GFi = malloc(2*NNZ*sizeof(double));
							double* big_XXt = malloc(2*NNZ*sizeof(double));
							memcpy ( big_I, *I, NNZ*sizeof(unsigned) );
							memcpy ( big_GFi, *GFi, NNZ*sizeof(double) );
							memcpy ( big_XXt, *XXt, NNZ*sizeof(double) );
							free(*GFi); free(*I); free(*XXt);
                            NNZ *= 2;
							*GFi = big_GFi;  *I = big_I; *XXt = big_XXt;
							//printf("Increased memory in mex calcGradAndActiveSet \n");
						}
						nnz_local = nnz++;
					}
					(*I)[nnz_local] = k_inner;
					(*GFi)[nnz_local] = GFi_val;
					(*XXt)[nnz_local] = XX_val;
					
				}
			}
		}
	}
	free(block);
	return nnz;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
 
    double *X = (double *)mxGetPr(prhs[0]); 
	size_t n = mxGetN(prhs[0]);
	size_t m = mxGetM(prhs[0]);
	
	double *Y = (double *)mxGetPr(prhs[1]); 

	mwIndex *Ai_rows = mxGetIr(prhs[2]);
    mwIndex *Ai_starts = mxGetJc(prhs[2]);
    double* Ai_vals = mxGetPr(prhs[2]); 
	
	double* invAi = (double *)mxGetPr(prhs[3]);

	
	double lambda = *mxGetPr(prhs[4]);
	double* idxs = (double *)mxGetPr(prhs[5]);
	size_t size_idxs = max(mxGetM(prhs[5]),mxGetN(prhs[5]));
	size_t total_nnz,i;
	unsigned NNZ = (unsigned)*mxGetPr(prhs[6]);
	
	int type = (int)*mxGetPr(prhs[7]);
	
	
	
	unsigned* I = malloc(NNZ*sizeof(unsigned));
	double* GFi = malloc(NNZ*sizeof(double));
	double* XXt = malloc(NNZ*sizeof(double));
	
	double* GF_ans;
	double* XXt_ans;
	double* I_ans;
	double* J_ans;
	
	if ((I==0)||(GFi==0)||(XXt==0))
	{
		printf("ERROR: Out of memory in mulXXT\n");
		return;
	}
	total_nnz = calcGradAndActive(X,Y, n, m, invAi, Ai_rows, Ai_starts, Ai_vals, lambda, idxs , size_idxs, NNZ, &I, &GFi, &XXt,type);
	
	plhs[0] = mxCreateDoubleMatrix(total_nnz, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(total_nnz, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(total_nnz, 1, mxREAL);
	I_ans = (double*)mxGetPr(plhs[0]); 
    GF_ans = (double*)mxGetPr(plhs[1]); 
    XXt_ans = (double*)mxGetPr(plhs[2]);
	
	for (i=0 ; i < total_nnz ; i++){
		I_ans[i] = (double)I[i]+1;
	}
	memcpy (GF_ans, GFi, sizeof(double)*total_nnz);
	memcpy (XXt_ans, XXt, sizeof(double)*total_nnz);
	free(GFi); free(I); free(XXt); 
}
