#ifndef __BCDIC__H__
#define __BCDIC__H__

#include "mex.h"
#include <omp.h>



void mult_sparse_AtB(double *result, double *in1, double* in2, int colsize, double *set, int setsize, double *row_translation, double *row_translation_inverted, int rows, int cols, int flag);

void soft_shrinkage_Free(double *result, double lambda, double* P ,double *X, double* Grad, double epsilon, double* FreeSet, int elems);

void mult_SparseMatTranspose_DenseMat(double *result, double *in_sparse, double* set_sparse ,int set_sparse_size,int cols_sparse, double *in_dense, int cols_dense, int rows_dense);

void mult_AjWj(double *result, double *in1, double *in2, int colsize, double *set, int setsize, double *row_translation_inverted, double *row_translation, int rows, int cols);

void calcParamsForQuadLinesearch(double* ASi, double* Delta_i, double* M_Delta_i, 
        double* Gi_before_shrink, double* Free_Neighbors, double* Free_idxs, int size_Neighbors, int size_idxs,
        double* ans_a1, double* ans_a2, double* ans_ASi, double* ans_Delta_i);
#endif
