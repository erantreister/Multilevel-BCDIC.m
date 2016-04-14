#ifndef PERMUTEMAT_H
#define PERMUTEMAT_H
#include <math.h>
#include <time.h>
#include "mex.h"

typedef double real; 
typedef unsigned short idx_type;

typedef struct{
    mwIndex ans_B;
    mwIndex ans_B_lower;
    mwIndex ans_C;    
} Count;

typedef struct{
    mwIndex *starts;
    idx_type *cols;
    double *vals;
	double *invDiag; 
    mwIndex vc_idx;
} B_matrix;

typedef struct{
    idx_type *rows;
    idx_type *cols;
    double *vals;
    mwIndex *blocks;
} C_matrix;

typedef struct{
    mwIndex *cols;
    mwIndex *starts;
    double *vals;
    mwIndex nzmax;
    mwIndex* perm;
    mwIndex *reversed_perm;
    mwIndex* map;
    mwIndex* coreBlocks;
    int numCoreBlocks;
    int type;
    mwIndex n;
} Args;


void setup(Args *args, B_matrix *B, C_matrix *c);

void count_ans_create_C_blocks(Args *args, C_matrix *c, Count *count);

void calculateB(Args *args, B_matrix *b, C_matrix *c, Count *count);

void calculateBLower(Args *args, B_matrix *b_lower, C_matrix *c, Count *count);

void print(int _size, void* _array, int _type);


void multiplySparseHybridMatVec(Args *args, B_matrix *b, C_matrix *c, double* x_out, real *x_in,int cb_start, int cb_end);

void multiplySparseHybridMatLowerVec(Args *args, B_matrix *b_lower, C_matrix *c, double* x_out, real *x_in, int cb_start, int cb_end);

void multiplySparseMatVec(mwIndex *C, mwIndex *starts, double* vals, mwIndex istart, mwIndex iend, real* x_in, double* x_out);

void multiplySparseMatVecSymTriu(mwIndex *C, mwIndex *starts, double* vals, mwIndex n, real* x_in, double* x_out);

void free_memory(Args args, B_matrix B, C_matrix c);

#endif /* PERMUTEMAT_H */