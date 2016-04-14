
#if !defined(_WIN32)
#define ddot ddot_
#endif

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"


int multOffDiagonalBlock(double* X, double* temp_block, int n_rows, int n_cols, int m, int** I, int** J, double** V, double lambda, int maxNNZ, int* nnz, int offset_row, int offset_col){
	int k;
	char *trA = "T";
	char *trB = "N";
	double one = 1.0, zero = 0.0;
	
	size_t u_n_rows = n_rows,u_n_cols = n_cols, u_m = m;
	
	
	dgemm(trA, trB, &u_n_rows, &u_n_cols, &u_m, &one, X +offset_row*m , &u_m, X+offset_col*m, &u_m, &zero, temp_block, &u_n_rows);
	for (k = 0; k < n_rows*n_cols; k++){
		double v = 0;
		int i,j;
		int one = 1;
		
		
		v = temp_block[k];
		i = k % n_rows;
		j = k / n_rows;

		i = i + offset_row;
		j = j + offset_col;
		if (*nnz >= maxNNZ){
			int*    big_I = malloc(2*maxNNZ*sizeof(int));
			int*    big_J = malloc(2*maxNNZ*sizeof(int));
			double* big_V = malloc(2*maxNNZ*sizeof(double));
			memcpy (big_I, *I, maxNNZ*sizeof(int) );
			memcpy (big_J, *J, maxNNZ*sizeof(int) );
			memcpy (big_V, *V, maxNNZ*sizeof(double) );
			free(*I); free(*J); free(*V);
            maxNNZ *= 2;
			*I = big_I;  *J = big_J; *V = big_V;
			//printf("Increased memory in mex mulXXt_sparse_new\n");
		}
		
		if (fabs(v) > lambda){
			int nnz_local;
			nnz_local = (*nnz)++;
			(*I)[nnz_local] = i;
			(*J)[nnz_local] = j;
			(*V)[nnz_local] = v;
		}
	}
    return maxNNZ;
}
int multDiagonalSymmetricSquareBlock(double* X,double* temp_block, int n, int m, int** I, int** J, double** V, double lambda, int maxNNZ, int* nnz, int offset){
	int k;
	char *trA = "T";
	char *trB = "N";
	double one = 1.0, zero = 0.0;
	
	size_t u_n = n, u_m = m;
	
	
	dgemm(trA, trB, &u_n, &u_n, &u_m, &one, X +offset*m , &u_m, X+offset*m, &u_m, &zero, temp_block, &u_n);
	for (k = 0; k < n*n; k++){
		double v = 0;
		int l, i,j,im,jm;
		int one = 1;
		
		
		v = temp_block[k];
		i = k % n;
		j = k / n;
		if (i>j){
			continue;
		}
		i = i + offset;
		j = j + offset;
		if (*nnz >= maxNNZ){
			int*    big_I = malloc(2*maxNNZ*sizeof(int));
			int*    big_J = malloc(2*maxNNZ*sizeof(int));
			double* big_V = malloc(2*maxNNZ*sizeof(double));
			memcpy (big_I, *I, maxNNZ*sizeof(int) );
			memcpy (big_J, *J, maxNNZ*sizeof(int) );
			memcpy (big_V, *V, maxNNZ*sizeof(double) );
			free(*I); free(*J); free(*V);
            maxNNZ *= 2;
			*I = big_I;  *J = big_J; *V = big_V;
			//printf("Increased memory in mex mulXXt_sparse_new\n");
		}
		
		if (fabs(v) > lambda){
			int nnz_local;
			if (i==j) // Elements in the diagonal are multiplied by 0.5 for the future creation of the Symmetric sparse matrix.
				v /=2;
			nnz_local = (*nnz)++;
			(*I)[nnz_local] = i;
			(*J)[nnz_local] = j;
			(*V)[nnz_local] = v;
		}
	}
    return maxNNZ;
}


int mulXXT_new(double* X, int n, int m, int** I, int** J, double** V, double lambda, int maxNNZ)
{
	int nnz = 0;
	int blockSize =  2048;
	int num_blocks = n/blockSize + (n%blockSize > 0);
	int k_block;
	
	double* block = malloc(blockSize*blockSize*sizeof(double));
	for (k_block = 0; k_block < (num_blocks*(num_blocks+1))/2; k_block ++){
		int tmp;
		int num_cols = blockSize, num_rows = blockSize;
		int l, i;
		
		
		int low = 0, high = num_blocks-1, j = (low+high)/2;
		tmp = ((num_blocks+(num_blocks-j))*(j+1)) / 2;
		while(low<=high)
		{
			if (tmp>k_block && k_block >=tmp-(num_blocks-j) )
				break;
			else if (tmp<=k_block)
				low = j+1;
			else
				high = j-1;
			j = (low+high)/2;
			tmp = ((num_blocks+(num_blocks-j))*(j+1)) / 2;
		}
		
		i = num_blocks - (tmp - k_block); // i,j are the block numbers.
		
		if (i==num_blocks - 1){
			num_rows = n - i*blockSize;
		}
		if (j==num_blocks - 1){
			num_cols = n - j*blockSize;
		}
		
		if (i==j){
			maxNNZ = multDiagonalSymmetricSquareBlock( X,block, num_rows,  m,  I,  J, V, lambda, maxNNZ, &nnz, i*blockSize);
		}else{
			maxNNZ = multOffDiagonalBlock(X,block, num_rows, num_cols, m, I, J, V, lambda, maxNNZ, &nnz, i*blockSize, j*blockSize);
		}
	}
	
	free(block);
	return nnz;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input params:
    double *X;
	// output params:
    double *Ians, *Jans, *Vans, lambda,NNZ_factor;
    int n,k;
	int i;
	int *I, *J;
	double *V;
	int NNZ = 0, total_nnz = 0;
	
	
    X = (double *)mxGetPr(prhs[0]);
	lambda = *mxGetPr(prhs[1]);
	NNZ_factor = *mxGetPr(prhs[2]);

	n = mxGetN(prhs[0]);
	k = mxGetM(prhs[0]);
	NNZ = ((int)NNZ_factor)*n;

	I = malloc(NNZ*sizeof(int));
	J = malloc(NNZ*sizeof(int));
	V = malloc(NNZ*sizeof(double));
	if ((I==0)||(J==0)||(V==0))
	{
		printf("ERROR: Out of memory in mulXXT\n");
		return;
	}	
	
	total_nnz = mulXXT_new(X, n, k, &I, &J, &V, lambda, NNZ);
	
	plhs[0] = mxCreateDoubleMatrix(total_nnz, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(total_nnz, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(total_nnz, 1, mxREAL);
	Ians = (double *)mxGetPr(plhs[0]); 
    Jans = (double *)mxGetPr(plhs[1]); 
    Vans = (double *)mxGetPr(plhs[2]);
    
    
	for (i=0 ; i < total_nnz ; i++){
		Ians[i] = (double)I[i]+1;
		Jans[i] = (double)J[i]+1;
	}
	memcpy (Vans, V, sizeof(double)*total_nnz);
	free(V); free(I); free(J); 
}
