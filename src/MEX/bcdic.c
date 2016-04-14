#include "bcdic.h"


void mult_sparse_AtB(double *result, double *in1, double* in2, int colsize, double *set, int setsize, double *row_translation, double *row_translation_inverted, int rows, int cols, int flag)
{
    int i;
    
    //iterate every non-zero in the result
    #pragma omp parallel for schedule(dynamic) private(i)
    for (i = 0; i < setsize; i++)
    {
        double solution = 0;
        int j, r = 0, c = 0, si = (int)set[i]-1;
        
        // Translate which is the row of "in1" from the row of result
        r = (int)row_translation[si % rows]-1;
        
        // The column of "in2" is the same as in result
        c = si / rows;

		if (flag && r<c)
			continue;

		for (j = 0; j < colsize; j++)
        {
            solution += in1[j+colsize*r] * in2[j+colsize*c];
        }
        
		if (flag)
		{
			result[(int)row_translation_inverted[c]-1+r*colsize] = solution;
		}
		
        result[si] = solution;
    }
    
}

void soft_shrinkage_Free(double *result, double lambda, double* P ,double *X, double* Grad, double epsilon, double* FreeSet, int elems)
{
    int i, elem;
	double tmp,c;
	#pragma omp parallel for schedule(dynamic) private(elem, tmp, c, i)
    for (i = 0; i < elems; ++i)
    {
		elem = ((int)FreeSet[i])-1;
		tmp = X[elem] - Grad[elem]*P[elem];
        c = lambda*P[elem];
        result[elem] = (tmp > c+epsilon) * (tmp - c) + (tmp < -c-epsilon) * (tmp + c) - X[elem];
    }
}

void mult_AjWj(double *result, double *in1, double *in2, int colsize, double *set, int setsize, double *row_translation_inverted, double *row_translation, int rows, int cols)
{
    int i,j,k;
    
    
    #pragma omp parallel for private(j, k , i)
    for (j = 0; j < cols; ++j)
    {
        int curr_col = 0;
        double acc = 0.0;
		int coljStarts = rows*j;
		int c,r,pos;
		i = 0;
		
		for (k = 0; k < setsize; ++k)
        {
            // Compute the off-diagonal block  = (sparse matrix)*(dense matrix)
			pos = (int)set[k]-1;
            r = pos % rows; // the row in invAi (in2) - the col is "j"
            c = pos / rows;
			
			if ((int)row_translation[r] == 0){ // this means that we're in an off-diagonal block
				result[r + coljStarts] += in1[k] * in2[(int)row_translation_inverted[c]-1 + coljStarts];
			}
			//Compute the diagonal block = (sparse matrix)'*(dense matrix)
			
            if (c != curr_col)
            {
                result[(int)row_translation_inverted[i]-1+coljStarts] = acc;
                acc = 0.0;
                // next row that is part of the block diagonal 
                i++;
            }
            acc += in1[k] * in2[r + coljStarts];
            curr_col = c;
        }
        result[(int)row_translation_inverted[i]-1+coljStarts] = acc; // don't forget to save the last value!
	}
}

void mult_SparseMatTranspose_DenseMat(double *result, double *in_sparse, double* set_sparse ,int set_sparse_size,int cols_sparse, double *in_dense, int cols_dense, int rows_dense){
    int i,j,k;

    // Compute the diagonal block = (sparse matrix)'*(dense matrix)
    #pragma omp parallel for schedule(dynamic) private(j, k, i)
    for (j = 0; j < cols_dense; ++j)
    {
        // iterate over non-zeros in Gi (in_sparse)
        int curr_col = 0;
        double acc = 0.0;
        // find the first row that is part of the block diagonal (TODO: bring the rows as a set)
        for (k = 0; k < set_sparse_size; ++k)
        {
            int pos = (int)set_sparse[k]-1;
            int r = pos % rows_dense; // the row in invAi (in_dense) - the col is "j"
            int c = pos / rows_dense;
            if (c != curr_col)
            {
                result[curr_col+j*cols_sparse] = acc;
                acc = 0.0;
                // next row that is part of the block diagonal 
            }
            acc += in_sparse[k] * in_dense[r + rows_dense*j];
            curr_col = c;
        }
        result[curr_col+j*cols_sparse] = acc; // don't forget to save the last value!
    }
}


void calcParamsForQuadLinesearch(double* ASi, double* Delta_i, double* M_Delta_i,
        double* Gi_before_shrink, double* Free_Neighbors, double* Free_idxs, int size_Neighbors, int size_idxs,
        double* ans_a1, double* ans_a2, double* ans_ASi, double* ans_Delta_i){
    
    int k,index;
    double a1 = 0, a2 = 0;
// preparing the vectors for the l-1:
    for (k = 0; k < size_Neighbors; k++){
        index = (int)Free_Neighbors[k]-1;
        ans_ASi[k] = 2*ASi[index];
        ans_Delta_i[k] = 2*Delta_i[index];
    }
    for (k = 0; k < size_idxs; k++){
        index = (int)Free_idxs[k]-1;
        ans_ASi[k + size_Neighbors] = ASi[index];
        ans_Delta_i[k + size_Neighbors] = Delta_i[index];
    }
// applying inner products with the l-1 data
    for (k = 0; k < size_Neighbors; k++){
        index = (int)Free_Neighbors[k]-1;
        a1 += Gi_before_shrink[index]*ans_Delta_i[k];
    }
    for (k = 0; k < size_idxs; k++){
        index = (int)Free_idxs[k]-1;
        a1 += Gi_before_shrink[index]*ans_Delta_i[k+size_Neighbors];
    }
    
    for (k = 0; k < size_Neighbors; k++){
        index = (int)Free_Neighbors[k]-1;
        a2 += M_Delta_i[index]*ans_Delta_i[k];
    }
    for (k = 0; k < size_idxs; k++){
        index = (int)Free_idxs[k]-1;
        a2 += M_Delta_i[index]*ans_Delta_i[k+size_Neighbors];
    }
    
    a2 *= 0.5;
    *ans_a1 = a1;
    *ans_a2 = a2;
}