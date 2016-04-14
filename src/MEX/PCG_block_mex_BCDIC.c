#include "mex.h"
#include <omp.h>
#include <math.h>
#include <time.h>

void applyInvLt(double* x, double* b, mwIndex n, mwIndex *R_L, mwIndex *starts_L, double* vals_L){
	// x - output so that Lt*x = b. We hold L column-wise.
	mwIndex j,k,i;
	double temp, invLii;
	for (k = n ; k > 0 ; --k){
		// k should be k-1 This is the true iterator, we just want to make sure that the iteration ends. MwIndex is unsigned.
        i = k-1;
        temp = b[i]; // the usage of temp allows the situation where b=x.
        /// WARNING: IF NOT INVERTING IN MAIN MEX THEN THIS LINE SHOULD BE (1/vals_L[starts_L[k-1]]);
		invLii = vals_L[starts_L[i]];
        for (j = starts_L[i]+1 ; j < starts_L[k] ; ++j){
			temp -= vals_L[j]*x[R_L[j]];
		}
		//the value in vals_L[starts_L[i]] is Lii;
		x[i] = temp*invLii;
	}    
}

void applyInvUt(double* x, double* b, mwIndex n, mwIndex *R_U, mwIndex *starts_U, double* vals_U){
	// x - output so that Ut*x = b. We hold U column-wise.
	mwIndex j,i;
	double temp;
    
    for (i = 0 ; i < n ; i++){ 
		temp = b[i]; 
		for (j = starts_U[i] ; j < starts_U[i+1]-1 ; ++j){
			temp -= vals_U[j]*x[R_U[j]];
		}
		/// WARNING: IF NOT INVERTING IN MAIN MEX THEN THIS LINE SHOULD BE temp*(1/vals_U[j]);
        x[i] = temp*vals_U[j];// j here is equal to starts_U[i+1]-1
        
	}
}




void applyInvL(double* x, double* b, mwIndex n, mwIndex *R_L, mwIndex *starts_L, double* vals_L){
	// x - output so that L*x = b. We hold L column-wise.
	mwIndex i,j;
	double temp;
	if (x!=b){
		for ( i = 0 ; i < n ; i++){
			x[i] = b[i];
		}
	}
	for ( i = 0 ; i < n ; i++){ // x[0...i-1] - the value of x, x[i...n] - the value of r.
        /// WARNING: IF NOT INVERTING IN MAIN MEX THEN THIS LINE SHOULD BE x[i]/vals_L[starts_L[i]]; ;
		temp = x[i]*vals_L[starts_L[i]]; 
		for (j = starts_L[i]+1 ; j < starts_L[i+1] ; ++j){
			x[R_L[j]] -= vals_L[j]*temp;
		}
		x[i] = temp;
	}
}


void applyInnerPCG(mwIndex* C_t, double* vals_t,  mwIndex* starts_t, mwIndex n ,double* x,  int max_iter, double tol, mwIndex *R_L, mwIndex *starts_L, double* vals_L, 
    mwIndex *R_U, mwIndex *starts_U, double* vals_U, double* r_aux, double* p_aux, double* z_aux){
    // this function treats one linear system, and uses the auxiliary vectors.
	
	double temp = 0, rnorm = 0, ztr = 0, ztr_new = 0, alpha, beta, *tmp_pointer,time,acc_time = 0 , acc_time_A = 0;
    mwIndex iter = 0,j,l;
    
    
	for (j=0 ; j<n ; j++){
        rnorm += r_aux[j]*r_aux[j];
    }
    if (rnorm < 1e-16){
        return; 
    }
	
	if (R_L==NULL){ // No preconditioner:
		for (j=0 ; j<n ; j++){
			p_aux[j] = r_aux[j];
		}
		ztr = rnorm; 
	}else{ // Triangular Preconditioner
		
		ztr = 0;
        if (R_U!=NULL){
            applyInvUt(z_aux, r_aux, n, R_U, starts_U, vals_U);
        }else{
            applyInvL (z_aux, r_aux, n, R_L, starts_L, vals_L);
        }	
		applyInvLt(z_aux, z_aux, n, R_L, starts_L, vals_L);
		for (j=0 ; j<n ; j++){
			p_aux[j] = z_aux[j];
			ztr += z_aux[j]*r_aux[j];
		}
	}	
    tol = (tol*tol)*rnorm; // we dont apply the sqrt of the norm here...
    
    while (iter < max_iter){
        iter++;
        memset(z_aux, 0, n*sizeof(double));
        for (j=0 ; j<n ; j++){
            temp = 0;
            for (l = starts_t[j] ; l < starts_t[j+1]-1 ; ++l)
            {
                temp += vals_t[l]*p_aux[C_t[l]];
                z_aux[C_t[l]] += vals_t[l]*p_aux[j];
            }
            temp += vals_t[l]*p_aux[C_t[l]];
            z_aux[j] += temp; // this is A*p we actually perform A'*p, but A is symmetric
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
		if (R_L==NULL){ // No preconditioner:
			ztr_new = rnorm;
            tmp_pointer = r_aux;
		}else{
			ztr_new = 0;
            if (R_U!=NULL){
                applyInvUt(z_aux, r_aux, n, R_U, starts_U, vals_U);
            }else{
                applyInvL(z_aux, r_aux, n, R_L, starts_L, vals_L);
            }
            applyInvLt(z_aux, z_aux, n, R_L, starts_L, vals_L);
            tmp_pointer = z_aux;
            for (j=0 ; j<n ; j++){
				ztr_new += z_aux[j]*r_aux[j];
			}
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
    mwIndex i;
    long j;
    mwIndex *R_L = 0;
    mwIndex *starts_L = 0;
    double* vals_L = 0;
    int inverted_U = 0, inverted_L = 0; 
    
    mwIndex *R_U = 0;
    mwIndex *starts_U = 0;
    double* vals_U = 0;
          
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
    double* p_aux = 0;
    double* z_aux = 0, *x_aux = 0;
    
    int num_of_threads = 1;
    //////////////////////
    // Handling the input r: if its a row vector - then its a list of canonical vectors.
	if (mxGetM(prhs[2])==1){
		r_canonical = 1;
		if (x_in != 0){
			printf("ERROR: canonical + x_in not implemented\n");
			return;
		}
	}
	r = mxGetPr(prhs[2]);

    /* Output Variables */

	rowsToReturn = mxGetPr(prhs[8]);
	
	if (rowsToReturn!=0){
		numRowsToReturn = max(mxGetN(prhs[8]),mxGetM(prhs[8]));
		plhs[0] = mxCreateDoubleMatrix(numRowsToReturn, k, mxREAL);
		flag_rowsToReturn = 1;
		x_out = mxGetPr(plhs[0]);
	}else{
		plhs[0] = mxCreateDoubleMatrix(n, k, mxREAL);
		x_out = mxGetPr(plhs[0]);
		if (x_in != 0){
			for (i = 0 ; i < (n*k) ; i++){
				x_out[i] = x_in[i];
			}
		}
	}

    if (mxGetPr(prhs[5])!=0){
        R_L = mxGetIr(prhs[5]);
        starts_L = mxGetJc(prhs[5]);
        vals_L = mxGetPr(prhs[5]);
        for (i = 0 ; i < n ; i++){
            vals_L[starts_L[i]] = 1/vals_L[starts_L[i]];
        }
        inverted_L = 1;
    }
    
    if (mxGetPr(prhs[6])!=0){
        R_U = mxGetIr(prhs[6]);
        starts_U = mxGetJc(prhs[6]);
        vals_U = mxGetPr(prhs[6]);
        for (i = 0 ; i < n ; i++){
            vals_U[starts_U[i+1]-1] = 1/vals_U[starts_U[i+1]-1];
        }
        inverted_U = 1;
    }
    num_of_threads = (int)(*mxGetPr(prhs[7]));
    if (num_of_threads<=0){
        num_of_threads = 1;
    }


    r_aux = (double*)malloc(num_of_threads*n*sizeof(double));
    p_aux = (double*)malloc(num_of_threads*n*sizeof(double));
    z_aux = (double*)malloc(num_of_threads*n*sizeof(double));
	if (flag_rowsToReturn){
		x_aux = (double*)malloc(num_of_threads*n*sizeof(double));
	}
	
	
            
    #pragma omp parallel for schedule(dynamic) num_threads(num_of_threads)
    for (j=0 ; j<k ; j++){
		int id = omp_get_thread_num();
		double* x_pointer;
		if (r_canonical){
			memset(r_aux + id*n, 0, n*sizeof(double));
			r_aux[(int)r[j]-1+id*n] = 1;
		}else{
			memcpy(r_aux + id*n, r+j*n, n*sizeof(double));
		}
		if (flag_rowsToReturn){
			if (x_in != 0){
				memcpy(x_aux + id*n, x_in + j*n, n*sizeof(double));
			}else{
				memset(x_aux + id*n, 0, n*sizeof(double));
			}
			x_pointer = x_aux + id*n;
		}else{
			x_pointer = x_out + j*n;
		}
        applyInnerPCG(C_t, vals_t, starts_t, n , x_pointer, max_iter, tol, R_L, starts_L, vals_L, R_U, starts_U, vals_U, r_aux+id*n,  p_aux+id*n,  z_aux+id*n);
		if (flag_rowsToReturn){
			int i = 0;
			for (i = 0; i < numRowsToReturn ; i++){
				x_out[j*numRowsToReturn + i] = x_pointer[(int)rowsToReturn[i]-1];
			}
		}
    }
	if (flag_rowsToReturn){
		free(x_aux);
	}
    free(r_aux);
    free(p_aux);
    free(z_aux);
   if (inverted_U){
       for (i = 0 ; i < n ; i++){
            vals_U[starts_U[i+1]-1] = 1/vals_U[starts_U[i+1]-1];
        }
   }
   if (inverted_L){
       for (i = 0 ; i < n ; i++){
            vals_L[starts_L[i]] = 1/vals_L[starts_L[i]];
        }
   }
}