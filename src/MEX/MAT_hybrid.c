#include "MAT_hybrid.h"

void setup(Args *args, B_matrix *B, C_matrix *c){
    mwIndex i;
	Count count = {0};
    mwIndex numBlocks = 0;
    mwIndex n = args->n;
	mwIndex* coreBlocks2, *new_perm;
    B->invDiag = 0;
	for (i = 0 ; i < n ; i++){
        if (args->map[i] > numBlocks){
			numBlocks = args->map[i];
		}
    }
	numBlocks++;
	coreBlocks2 =  (mwIndex*) malloc((numBlocks+1) * sizeof (mwIndex));
	new_perm =  (mwIndex*) malloc(n * sizeof (mwIndex));
	
	for (i = 0 ; i < (numBlocks+1); i++){
        coreBlocks2[i]=0;
    }
	for (i = 0 ; i < n ; i++){
        coreBlocks2[args->map[i]+1]++ ;
    }
	for (i = 1 ; i < (numBlocks+1); i++){
        coreBlocks2[i]=coreBlocks2[i] + coreBlocks2[i-1];
    }	
	for (i = 0 ; i < n ; i++){
		new_perm[coreBlocks2[args->map[i]]] = i;
		coreBlocks2[args->map[i]]++;
    }
	for (i = 0 ; i < (numBlocks+1); i++){
        coreBlocks2[i]=0;
    }
	for (i = 0 ; i < n ; i++){
        coreBlocks2[args->map[i]+1]++ ;
    }
	for (i = 1 ; i < (numBlocks+1); i++){
        coreBlocks2[i]=coreBlocks2[i] + coreBlocks2[i-1];
    }
	coreBlocks2[numBlocks] = n;
	
	args->coreBlocks = coreBlocks2;
	args->perm = new_perm;
	args->numCoreBlocks = numBlocks;
    
    args->reversed_perm = (mwIndex*) malloc(args->n * sizeof (mwIndex));
    c->blocks = (mwIndex*) malloc(sizeof(mwIndex)*(args->numCoreBlocks*args->numCoreBlocks+1));
    B->vc_idx = 0;
    B->starts = (mwIndex*)malloc((n+1) * sizeof (mwIndex));
 
    for(i = 0; i < n ; i++){
        args->reversed_perm[args->perm[i]] = i;
    }
    
    for(i = 0; i < (args->numCoreBlocks*args->numCoreBlocks+1); i++){
		c->blocks[i] = 0;
    }
    
    count_ans_create_C_blocks(args, c, &count);
    // note that currently, the size of block ij is held at C_blocks[i*args.numCoreBlocks+j+1]
    
    c->cols = (idx_type*)malloc(count.ans_C * sizeof (idx_type));
    c->rows = (idx_type*)malloc(count.ans_C * sizeof (idx_type));
    c->vals = (double*)malloc(count.ans_C * sizeof (double));
    
    
     if(args->type){
        B->cols = (idx_type*) malloc(count.ans_B * sizeof (idx_type));
        B->vals = (double*) malloc(count.ans_B * sizeof (double));
        calculateB(args, B, c, &count);
    } else{        
        B->cols = (idx_type*)malloc(count.ans_B_lower * sizeof (idx_type));
        B->vals = (double*)malloc(count.ans_B_lower * sizeof (double));
        calculateBLower(args, B, c, &count);
    }
}

void calc_invDiag_B(Args *args, B_matrix *B){
mwIndex cb, offset, i,j;
B->invDiag = (double*)malloc(args->n*sizeof (double));
for (cb = 0 ; cb < args->numCoreBlocks ; cb++){
	offset = args->coreBlocks[cb];
	for (i = args->coreBlocks[cb] ; i < args->coreBlocks[cb+1] ; i++) {            
		for (j = B->starts[i] ; j < B->starts[i+1] ; ++j) {
			if (B->cols[j]+offset == i){ // this is a diagonal element
				B->invDiag[i] = 1/B->vals[j];
			}	
		}
	}
}
}

void count_ans_create_C_blocks(Args *args, C_matrix *c, Count *count){
    mwIndex i,j,p,cluster_p, r, l, col;
    mwIndex size_of_C_blocks = args->numCoreBlocks;
    for (i = 0 ; i < args->n ; i++){
        p = args->reversed_perm[i];
        cluster_p = args->map[i];
        r = args->coreBlocks[cluster_p+1];
        l = args->coreBlocks[cluster_p];
        
        for (j = args->starts[i]; j < args->starts[i+1]; j++){
            col = args->reversed_perm[args->cols[j]];
            if (col < r && col >= l){
                (*count).ans_B++;
                (*count).ans_B_lower += (i <= args->cols[j]);
            }else{
                c->blocks[args->map[i]*size_of_C_blocks + args->map[args->cols[j]]+1]++;
                (*count).ans_C++;
            }
        }
    }
    c->blocks[0] = 0;
	for(i = 1; i < (size_of_C_blocks*size_of_C_blocks+1); i++){
		c->blocks[i] = c->blocks[i-1] + c->blocks[i];
	}    
}

void calculateB(Args *args, B_matrix *B, C_matrix *c, Count *count){
mwIndex i,j,p,k,cluster_p,l,r,col, n = args->n;          
mwIndex size_of_C_blocks = args->numCoreBlocks;
for (i = 0 ; i < n ; i++){
        B->starts[i] = B->vc_idx;
        p = args->perm[i];
        cluster_p = args->map[p];
		l = args->coreBlocks[cluster_p];
        r = args->coreBlocks[cluster_p+1];
        for (j = args->starts[p]; j < args->starts[p+1]; j++){
            col = args->reversed_perm[args->cols[j]];
            if (col < r && col >= l){
                B->cols[B->vc_idx] = (idx_type)(col - l);
                B->vals[B->vc_idx] = args->vals[j];
                (B->vc_idx)++;
            }else{
			    mwIndex cluster_col = args->map[args->cols[j]];
                k = c->blocks[cluster_p*size_of_C_blocks + cluster_col]++;
                c->vals[k] = args->vals[j];
                c->rows[k] = (idx_type)(i - l);
                c->cols[k] = (idx_type) (col - args->coreBlocks[cluster_col]);
            }
        }
    }
	B->starts[n] = B->vc_idx;
	
    //Correct C-blocks which now holds "ends" back to "starts".
	for(i = (size_of_C_blocks*size_of_C_blocks); i > 0; i--){
		c->blocks[i] = c->blocks[i-1];
	}
	c->blocks[0] = 0;
	c->blocks[(size_of_C_blocks*size_of_C_blocks)] = count->ans_C;
}

void calculateBLower(Args *args, B_matrix *b_lower, C_matrix *c, Count *count){
mwIndex i,j,p,k,cluster_p,l,r,col, n = args->n;         
mwIndex size_of_C_blocks = args->numCoreBlocks;
    for (i = 0 ; i < n ; i++){
        b_lower->starts[i] = b_lower->vc_idx;
        p = args->perm[i];
        cluster_p = args->map[p];
		l = args->coreBlocks[cluster_p];
        r = args->coreBlocks[cluster_p+1];
        for (j = args->starts[p]; j < args->starts[p+1]; j++){
            col = args->reversed_perm[args->cols[j]];
            if (col < r && col >= l){
                if(i <= col){
                   b_lower->cols[b_lower->vc_idx] = (idx_type)(col) - l;
                   b_lower->vals[b_lower->vc_idx] = args->vals[j];
                   b_lower->vc_idx++;
                }
            }else{
			    mwIndex cluster_col = args->map[args->cols[j]];
                k = c->blocks[cluster_p*size_of_C_blocks + cluster_col]++;
                c->vals[k] = args->vals[j];
                c->rows[k] = (idx_type)(i - l);
                c->cols[k] = (idx_type) (col - args->coreBlocks[cluster_col]);
            }
        }
    }
    b_lower->starts[n] = b_lower->vc_idx;
	
    //Correct C-blocks which now holds "ends" back to "starts"->
	for(i = (size_of_C_blocks*size_of_C_blocks); i > 0; i--){
		c->blocks[i] = c->blocks[i-1];
	}
	c->blocks[0] = 0;
	c->blocks[(size_of_C_blocks*size_of_C_blocks)] = count->ans_C;
}
/*************************************************************************/
void multiplySparseHybridMatVec(Args *args, B_matrix *B, C_matrix *c, real* x_out, real *x_in,int cb_start, int cb_end){
	mwIndex offset,offset2;
	int cb, cb2, i,j;
    real tmp;
    
	// Multiplying B:
    for (cb = cb_start ; cb < cb_end ; cb++){
		offset = args->coreBlocks[cb];
		for (i = args->coreBlocks[cb] ; i < args->coreBlocks[cb+1] ; i++) {
			tmp = 0;
			for (j = B->starts[i] ; j < B->starts[i+1] ; ++j) {					
				tmp += B->vals[j]*x_in[B->cols[j]+offset];
			}
			x_out[i] = tmp;
		}
	}
	// Multiplying C:
    for (cb = cb_start ; cb < cb_end ; cb++){
		offset = args->coreBlocks[cb];
		for (cb2 = 0 ; cb2 < args->numCoreBlocks ; cb2++){
			offset2 = args->coreBlocks[cb2];
			for (i = c->blocks[cb*args->numCoreBlocks+cb2] ; i < c->blocks[cb*args->numCoreBlocks + cb2 + 1] ; i++) {
                x_out[c->rows[i] + offset] += c->vals[i]*x_in[c->cols[i]+offset2];
			}
		}
	}	
}
     
void multiplySparseHybridMatLowerVec(Args *args, B_matrix *B, C_matrix *c, double* x_out, real *x_in, int cb_start, int cb_end){
	mwIndex i,j,offset,offset2;
    int cb, cb2;
    double tmp;
    for (cb = cb_start ; cb < cb_end ; cb++){
		offset = args->coreBlocks[cb];
		for (i = args->coreBlocks[cb] ; i < args->coreBlocks[cb+1] ; i++) {
			tmp = 0;
            
			for (j = B->starts[i]+1 ; j < B->starts[i+1] ; ++j) {
				tmp += B->vals[j]*x_in[B->cols[j]+offset];
                x_out[B->cols[j]+offset] += B->vals[j]*x_in[i];
			}
            j = B->starts[i];
            x_out[i] += tmp;
            x_out[i] += B->vals[j]*x_in[B->cols[j]+offset];	
		}
	}
    
    // Multiplying C:
	for (cb = cb_start ; cb < cb_end ; cb++){
		offset = args->coreBlocks[cb];
		for (cb2 = 0 ; cb2 < args->numCoreBlocks ; cb2++){
			offset2 = args->coreBlocks[cb2];
			for (i = c->blocks[cb*args->numCoreBlocks+cb2] ; i < c->blocks[cb*args->numCoreBlocks + cb2 + 1] ; i++) {
				x_out[c->rows[i] + offset] += c->vals[i]*x_in[c->cols[i]+offset2];
			}
		}
	}
}

void multiplySparseMatVec(mwIndex *C, mwIndex *starts, double* vals, mwIndex istart, mwIndex iend, real* x_in, double* x_out){
    // This function is only for CSR data
    mwIndex i;
	mwIndex j;
	double tmp;
    for (i=istart ; i<iend ; i++) {
        tmp = 0;
        for (j = starts[i] ; j < starts[i+1] ; ++j) {
            tmp += vals[j]*x_in[C[j]];
        }
        x_out[i] = tmp; // this is A*p we actually perform A'*p, but A is symmetric
    }
}

void multiplySparseMatVecSymTriu(mwIndex *C, mwIndex *starts, double* vals, mwIndex n, real* x_in, double* x_out){
    // This function is only for CSR data
    int i;
	mwIndex j;
	double tmp;
	memset(x_out, 0, n*sizeof(double));
    for (i=0 ; i<n ; i++) {
        tmp = 0;
        for (j = starts[i] ; j < starts[i+1]-1 ; ++j) {
			tmp += vals[j]*x_in[C[j]];
			x_out[C[j]] += vals[j]*x_in[i];
        }
		tmp += vals[j]*x_in[C[j]];
        x_out[i] += tmp; // this is A*p we actually perform A'*p, but A is symmetric
    }
}


void free_memory(Args args, B_matrix B, C_matrix c){
    free(args.reversed_perm);
	free(args.perm);
	free(args.coreBlocks);
    free(B.starts);
    free(B.cols);
    free(B.vals);
	if (B.invDiag!=0){
		free(B.invDiag);
	}
    free(c.cols);
    free(c.rows);
    free(c.vals);
    free(c.blocks);
}