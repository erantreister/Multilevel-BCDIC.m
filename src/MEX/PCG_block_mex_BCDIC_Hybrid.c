#include "mex.h"
#include <omp.h>
#include <math.h>
#include <time.h>
#include "MAT_hybrid.h"

#define PREC_ITER 1
#define OMEGA_SOR 1.0


void local_SSOR(Args *args, B_matrix *B, mwIndex cb, double* x, const double* r){
	mwIndex j;
	int i,iter;
	double tmp,r_norm;
	mwIndex offset = args->coreBlocks[cb];
	r_norm = 0;
	for (i = args->coreBlocks[cb] ; i < args->coreBlocks[cb+1] ; i++) {
		r_norm += r[i]*r[i];
	}
	if (r_norm < 1e-16){
		return;
	}
	for (iter = 0 ; iter < PREC_ITER ; iter++){
		for (i = args->coreBlocks[cb] ; i < args->coreBlocks[cb+1] ; i++) {
			tmp = r[i];
			for (j = B->starts[i] ; j < B->starts[i+1] ; ++j) {					
				tmp -= B->vals[j]*x[B->cols[j]+offset];
			}
			
			x[i] += tmp*B->invDiag[i]; // alternative if OMEGA_SOR is not 1: x[i] += OMEGA_SOR*tmp*B->invDiag[i];
		}
		for (i = args->coreBlocks[cb+1]-1 ; i >= (int)args->coreBlocks[cb] ; i--) {
			tmp = r[i];
			for (j = B->starts[i] ; j < B->starts[i+1] ; ++j) {					
				tmp -= B->vals[j]*x[B->cols[j]+offset];
			}
			x[i] += tmp*B->invDiag[i]; //alternative if OMEGA_SOR is not 1: x[i] += OMEGA_SOR*tmp*B->invDiag[i];
		}
	}
}
void local_SSOR_lower(Args *args, B_matrix *B, mwIndex cb, double* x, double* r_local){
	mwIndex j;
	int i,iter;
	double tmp,r_norm;
	mwIndex offset = args->coreBlocks[cb];
	r_norm = 0;
	for (i = 0 ; i < (args->coreBlocks[cb+1]-args->coreBlocks[cb]) ; i++) {
		r_norm += r_local[i]*r_local[i];
	}
	if (r_norm < 1e-16){
		return;
	}
	for (iter = 0 ; iter < PREC_ITER ; iter++){
		for (i = args->coreBlocks[cb] ; i < args->coreBlocks[cb+1] ; i++) {
			tmp = r_local[i-offset];
			for (j = B->starts[i] ; j < B->starts[i+1] ; ++j) {					
				tmp -= B->vals[j]*x[B->cols[j]+offset];
			}
			tmp /= B->vals[B->starts[i]]; // if OMEGA_SOR is not 1 : tmp *= (OMEGA_SOR/B->vals[B->starts[i]]);
			
			x[i] += tmp;
			for (j = B->starts[i]+1 ; j < B->starts[i+1] ; ++j) {					
				r_local[B->cols[j]] -= B->vals[j]*tmp;
			}
		}
		
		for (i = args->coreBlocks[cb+1]-1 ; i >= (int)args->coreBlocks[cb] ; i--) {
			tmp = r_local[i-offset];
			for (j = B->starts[i] ; j < B->starts[i+1] ; ++j) {					
				tmp -= B->vals[j]*x[B->cols[j]+offset];
			}
			tmp /= (B->vals[B->starts[i]]); //if OMEGA_SOR is not 1 : tmp *= (OMEGA_SOR/B->vals[B->starts[i]]);
			x[i] += tmp;
			for (j = B->starts[i]+1 ; j < B->starts[i+1] ; ++j) {					
				r_local[B->cols[j]] -= B->vals[j]*tmp;
			}
		}
	}
}

void applyPrec(Args *args, B_matrix *B, double* z_aux, const double* r_aux){
	int cb;
	memset(z_aux, 0, args->n*sizeof(double));
	if (args->type==1){
		for (cb = 0 ; cb < args->numCoreBlocks ; cb++){
			// see if norm is small enough...
			local_SSOR(args, B, cb, z_aux, r_aux);
		}
	}else{
		int max_block_size = 0;
		double* aux_local_vec;
		for (cb = 0 ; cb < args->numCoreBlocks ; cb++){
			if (args->coreBlocks[cb+1] - args->coreBlocks[cb] > max_block_size){
				max_block_size = args->coreBlocks[cb+1] - args->coreBlocks[cb];
			}
		}
		aux_local_vec = (double*)malloc(max_block_size*sizeof(double));
		for (cb = 0 ; cb < args->numCoreBlocks ; cb++){
			memcpy(aux_local_vec, r_aux+args->coreBlocks[cb], (args->coreBlocks[cb+1] - args->coreBlocks[cb])*sizeof(double));
			local_SSOR_lower(args, B, cb, z_aux, aux_local_vec);
		}
		free(aux_local_vec);
	}
}

void applyInnerPCG(Args *args, B_matrix *B, C_matrix *C, mwIndex n ,double* x, int use_prec, int max_iter, double tol, double* r_aux, real* p_aux, double* z_aux){
    // this function treats one linear system, and uses the auxiliary vectors.
	
	double rnorm = 0, ztr = 0, ztr_new = 0, alpha, *tmp_pointer,time,acc_time = 0 , acc_time_A = 0;
	double beta, temp = 0;
    mwIndex iter = 0,j;
    
    
	for (j=0 ; j<n ; j++){
        rnorm += r_aux[j]*r_aux[j];
    }
    if (rnorm < 1e-16){
        return; 
    }
	
	if (use_prec == 0){ // No preconditioner:
		for (j=0 ; j<n ; j++){
			p_aux[j] = (real)r_aux[j];
		}
		ztr = rnorm; 
	}else{ 
		ztr = 0;
        applyPrec(args, B, z_aux,r_aux);
		for (j=0 ; j<n ; j++){
			p_aux[j] = (real)z_aux[j];
			ztr += z_aux[j]*r_aux[j];
		}
	}	
    tol = (tol*tol)*rnorm; // we dont apply the sqrt of the norm here...
    
    while (iter < max_iter){
        iter++;
        if (args->type==0){
			memset(z_aux, 0, n*sizeof(double));
			multiplySparseHybridMatLowerVec(args, B, C, z_aux, p_aux, 0, args->numCoreBlocks);
		}else{
			multiplySparseHybridMatVec(args, B, C, z_aux, p_aux, 0, args->numCoreBlocks);
		}
        
        temp = 0;
        for (j=0 ; j<n ; j++){
            temp += p_aux[j]*z_aux[j]; // temp is p'*Ap
        }
        alpha = ztr / temp;
        rnorm = 0;
        for (j=0 ; j<n ; j++){
            x[j] += alpha*p_aux[j];
            r_aux[j] -= alpha*z_aux[j];
            rnorm += r_aux[j]*r_aux[j];
        }
        if (rnorm < tol){
            break;
        }
		
		if (use_prec == 0){
			ztr_new = rnorm;
            tmp_pointer = r_aux;
		}else{
			ztr_new = 0;
			applyPrec(args, B, z_aux,r_aux);
			for (j=0 ; j<n ; j++){
				ztr_new += z_aux[j]*r_aux[j];
			}
			tmp_pointer = z_aux;
		}
		
        
        beta = ztr_new/ztr;
        for (j=0 ; j<n ; j++){
            p_aux[j] = tmp_pointer[j] + beta*p_aux[j];
        }
        ztr = ztr_new;
    }
    
}

    
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int j;
    
	Args args;
    B_matrix B;
    C_matrix C;
	
          
    mwIndex *C_t = mxGetIr(prhs[0]);
    mwIndex *starts_t = mxGetJc(prhs[0]);
    double* vals_t = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    
    double* x_in = mxGetPr(prhs[1]);
	
    double* r = 0;
	int r_canonical = 0;
    mwIndex k = mxGetN(prhs[2]);
    int max_iter = (int)(*mxGetPr(prhs[3]));
    double tol = (double)(*mxGetPr(prhs[4]));
	
	double* rowsToReturn = 0;
	int flag_rowsToReturn = 0;
	int numRowsToReturn = 0;
	
    double* x_out = 0;
    
    double* r_aux = 0;
    real* p_aux = 0;
    double* z_aux = 0, *x_aux = 0;
    
    int num_of_threads = 1;
	int use_prec = 0;

    // Handling the input r: if its a row vector - then its a list of canonical vectors.
	if (mxGetM(prhs[2])==1){
		r_canonical = 1;
		if (x_in != 0){
			printf("ERROR: canonical + x_in not implemented\n");
			return;
		}
	}
	r = mxGetPr(prhs[2]);


	num_of_threads = (int)(*mxGetPr(prhs[5]));
    if (num_of_threads<=0){
        num_of_threads = 1;
    }
	rowsToReturn = mxGetPr(prhs[6]);
	
    /* Output Variables */
	if (rowsToReturn!=0){
		numRowsToReturn = max(mxGetN(prhs[6]),mxGetM(prhs[6]));
		plhs[0] = mxCreateDoubleMatrix(numRowsToReturn, k, mxREAL);
		flag_rowsToReturn = 1;
	}else{
		plhs[0] = mxCreateDoubleMatrix(n, k, mxREAL);
	}
	x_out = mxGetPr(plhs[0]);	
    
	///////////////////////////////////////////////////////////////////////////
	args.cols = C_t;
    args.starts = starts_t;
    args.vals = vals_t;
    args.map = (unsigned int*)mxGetPr(prhs[7]);
    args.type = (int)*mxGetPr(prhs[8]);
    use_prec = (int)*mxGetPr(prhs[9]);
	args.n = n;
    
	setup(&args, &B , &C);
	
	if (use_prec==1 && args.type == 1){
		calc_invDiag_B(&args, &B);
	}
	
    r_aux = (double*)malloc(num_of_threads*n*sizeof(double));
    p_aux = (real*)malloc(num_of_threads*n*sizeof(real));
    z_aux = (double*)malloc(num_of_threads*n*sizeof(double));
	x_aux = (double*)malloc(num_of_threads*n*sizeof(double));
	       
    #pragma omp parallel for schedule(dynamic) num_threads(num_of_threads)
    for (j=0 ; j<k ; j++){
		int i;
		int id = omp_get_thread_num();
		double* x_pointer = x_aux + id*n;
		double* r_pointer = r_aux + id*n;
		if (r_canonical){
			memset(r_pointer, 0, n*sizeof(double));
			r_pointer[args.reversed_perm[(int)r[j]-1]] = 1;
		}else{
			for (i=0 ; i < n ; i++){
				r_pointer[i] = r[j*n + args.perm[i]];
			}
		}
		if (x_in != 0){
			for (i=0 ; i < n ; i++){
				x_pointer[i] = x_in[j*n + args.perm[i]];
			}
		}else{
			memset(x_pointer, 0, n*sizeof(double));
		}
        applyInnerPCG(&args, &B, &C, n , x_pointer,use_prec, max_iter, tol, r_pointer,  p_aux+id*n,  z_aux+id*n);
		if (flag_rowsToReturn){
			for (i = 0; i < numRowsToReturn ; i++){
				x_out[j*numRowsToReturn + i] = x_pointer[args.reversed_perm[(int)rowsToReturn[i]-1]];
			}
		}else{
			for (i = 0; i < n ; i++){
				x_out[j*n + i] = x_pointer[args.reversed_perm[i]];
			}
		}
    }
	free(x_aux);
    free(r_aux);
    free(p_aux);
    free(z_aux);
	free_memory(args,B,C);
}