/*---------------------------------------------------------------------------
Copyright (2012): Eran Treister and Irad Yavneh. 
This code is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."
---------------------------------------------------------------------------*/
#include "mex.h"
#include <math.h>
#include <omp.h>
#include <time.h>

double sign(double a){
    return 2*(a >=0) - 1; // this returns 1 for 0 as well...
}
double Max(double a, double b);
double Min(double a, double b);
void quickSort(double numbers[], int indices[], int array_size);
void q_sort(double numbers[], int indices[], int left, int right);

void mergeSort(double numbers[], int indices[], int array_size);
void partition(double numbers[], int indices[],int low,int high, double tmp_numbers[], int tmp_indices[]);
void merge(double numbers[], int indices[], int low,int mid,int high, double tmp_numbers[], int tmp_indices[]);

// if subInd!=null, then n is the length of subInd
double ApplyLinesearch(double a0, double a1, double a2, double mu, double a, double b, double* x, double* v, double* subInd, int n,double* J_opt){
	double alpha_opt = 0,alpha_l,alpha_r,alpha,J;
	double* alpha_critic;
	int *indices_critic, num_critic = 0, up_bound_num_critic = 0;
	double sumx = 0, sumv = 0, rhox = 0 , rhov = 0, tmp;
	double sign_t, eps = 1e-15, eps_subsection = 1e-7;
	int i, debug = 0;
	
	double tic;
	
	
	if (subInd==0){
		for (i=0 ; i<n ; i++){
			up_bound_num_critic += (v[i]!=0.0);
		}
    }else{
		up_bound_num_critic = n;
	}
    
  	alpha_critic = (double *)malloc(sizeof(double)*(up_bound_num_critic));
    indices_critic = (int *)malloc(sizeof(int)*(up_bound_num_critic));
	
	if (subInd==0){
		for (i = 0 ; i < n ; i++){
			tmp = (v[i]!=0.0) ? -x[i]/v[i] : b+1;
			if (tmp < a + eps || tmp > b - eps){ 
				sign_t = sign(x[i] + (a+eps)*v[i]);
				sumx += sign_t*x[i];
				sumv += sign_t*v[i];
			}else{ // alpha is critic, in [a,b]
				sign_t = sign(x[i] + (a+eps)*v[i]);
				rhox += sign_t*x[i];
				rhov += sign_t*v[i];
				alpha_critic[num_critic] = tmp;
				indices_critic[num_critic] = i;
				num_critic++;
			}
		}
	}else{
		for (i = 0; i < n; i++) {
			int sub_i = (int)(subInd[i])-1;
			tmp = (v[sub_i] != 0.0) ? -x[sub_i] / v[sub_i] : b + 1;
			if (tmp < a + eps || tmp > b - eps) { 
				sign_t = sign(x[sub_i] + (a + eps) * v[sub_i]);
				sumx += sign_t * x[sub_i];
				sumv += sign_t * v[sub_i];
			} else { // alpha is critic, in [a,b]
				sign_t = sign(x[sub_i] + (a+eps) * v[sub_i]);
				rhox += sign_t * x[sub_i];
				rhov += sign_t * v[sub_i];
				alpha_critic[num_critic] = tmp;
				indices_critic[num_critic] = sub_i;
				num_critic++;
			}
		}
		
	}
	
	
	
	a0 += mu*sumx;
    a1 += mu*sumv;
	if (num_critic==0){
       
		tmp = -a1/(2*a2);
        tmp = Max(Min(tmp,b),a);
        *J_opt = a2*tmp*tmp + a1*tmp + a0;
        alpha_opt = tmp;
	}else{
		mergeSort(alpha_critic, indices_critic, num_critic);
		// quickSort(alpha_critic, indices_critic, num_critic);
		
		alpha_l = a;
		alpha_r = alpha_critic[0];
		tmp = (mu*rhov+a1);
		alpha_opt = -tmp/(2*a2);
		alpha_opt = Max(Min(alpha_r,alpha_opt),alpha_l);
		*J_opt = (a2*alpha_opt + tmp)*alpha_opt + a0 + mu*rhox;
		sign_t = sign(x[indices_critic[0]] + (a+eps)*v[indices_critic[0]]);
		rhox = rhox - 2*x[indices_critic[0]]*sign_t;
        rhov = rhov - 2*v[indices_critic[0]]*sign_t;
		i = 1;
		while (fabs(alpha_critic[i] - alpha_critic[i-1]) < eps_subsection && i<num_critic){
            sign_t = sign(x[indices_critic[i]] + (a+eps)*v[indices_critic[i]]);
			rhox = rhox - 2*x[indices_critic[i]]*sign_t;
			rhov = rhov - 2*v[indices_critic[i]]*sign_t;
            i++;
		}
		while (i<num_critic){
			alpha_l = alpha_critic[i-1];
			alpha_r = alpha_critic[i];
			tmp = (mu*rhov+a1);
			alpha = -tmp/(2*a2);
			alpha = Max(Min(alpha_r,alpha),alpha_l);
			J = (a2*alpha + tmp)*alpha + a0 + mu*rhox;
			if ( J <= *J_opt){
				*J_opt = J;
				alpha_opt = alpha;
			}else{// optional? else break. A shortcut because of the convexity.
                break;
            }
			sign_t = sign(x[indices_critic[i]] + (a+eps)*v[indices_critic[i]]); 
			rhox = rhox - 2*x[indices_critic[i]]*sign_t;
			rhov = rhov - 2*v[indices_critic[i]]*sign_t;
            i++;
			while (fabs(alpha_critic[i] - alpha_critic[i-1])< eps_subsection && i<num_critic){
				sign_t = sign(x[indices_critic[i]] + (a+eps)*v[indices_critic[i]]);
                rhox = rhox - 2*x[indices_critic[i]]*sign_t;
				rhov = rhov - 2*v[indices_critic[i]]*sign_t;
                i++;
			}
		}
		alpha_l = alpha_critic[num_critic-1];
		alpha_r = b;
		tmp = (mu*rhov+a1);
		alpha = -tmp/(2*a2);
		alpha = Max(Min(alpha_r,alpha),alpha_l);
		J = (a2*alpha + tmp)*alpha + a0 + mu*rhox;
		if ( J < *J_opt){
			*J_opt = J;
			alpha_opt = alpha;
		}
	}
	
	free(alpha_critic);
	free(indices_critic);
	return alpha_opt;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
    // input parameters:
    double a2, a1, a0, a, b, mu;
    double *x, *v, *subInd = 0;
    int n_xv,i;
    // output parameters:
    double *alpha_opt,*J_opt;
    // help variables
   
    x = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    
    a2 = *mxGetPr(prhs[2]);
    a1 = *mxGetPr(prhs[3]);
    a0 = *mxGetPr(prhs[4]);
    a = *mxGetPr(prhs[5]);
    b = *mxGetPr(prhs[6]);
    mu = *mxGetPr(prhs[7]);
	if (mxGetM(prhs[0])!=mxGetM(prhs[1])){
        printf("error: x and v are not from the same size! Aborting!");
        return;
    }
    if (mxGetN(prhs[0])!=mxGetN(prhs[1])){
        printf("error: x and v are not from the same size! Aborting!");
        return;
    }
	if (nrhs > 8){
		subInd = mxGetPr(prhs[8]);
		// n_xv is the length of the sub-indices vector. All the rest of the vectors are assumed to be zeros.
		if (mxGetM(prhs[8])==1){
			n_xv = (int)mxGetN(prhs[8]);
		}else{
			n_xv = (int)mxGetM(prhs[8]);
		}	
	}else{
		 n_xv = (int)(mxGetM(prhs[0])*mxGetN(prhs[0]));
	}
	
	/*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	/*Get a pointer to the data space in our newly allocated memory*/
    alpha_opt = mxGetPr(plhs[0]);
    J_opt     = mxGetPr(plhs[1]);

	// END OF MEX INTERFACE ISSUES.
    alpha_opt[0] =  ApplyLinesearch(a0, a1, a2,  mu,  a,  b,  x,  v, subInd, n_xv, J_opt);
}



void quickSort(double numbers[], int indices[], int array_size)
{
  q_sort(numbers, indices, 0, array_size - 1);
}
 
void q_sort(double numbers[], int indices[], int left, int right)
{
  int pivot_idx, l_hold, r_hold;
  double pivot;
  time_t t;
  srand((unsigned) time(NULL));
  l_hold = left;
  r_hold = right;
  pivot_idx = left + rand()%(right-left);
  pivot_idx = left;
  pivot = numbers[pivot_idx];
  pivot_idx = indices[pivot_idx];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      indices[left] = indices[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      indices[right] = indices[left];
      right--;
    }
  }
  numbers[left] = pivot;
  indices[left] = pivot_idx;
  pivot_idx = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot_idx)
    q_sort(numbers, indices, left, pivot_idx-1);
  if (right > pivot_idx)
    q_sort(numbers, indices, pivot_idx+1, right);
}


void mergeSort(double numbers[], int indices[], int array_size){
	double* tmp_numbers = (double *)malloc(sizeof(double)*(array_size));
    int* tmp_indices = (int *)malloc(sizeof(int)*(array_size));
	partition(numbers,indices,0,array_size-1,tmp_numbers,tmp_indices);
	free(tmp_numbers);
	free(tmp_indices);
}


void partition(double numbers[], int indices[],int low,int high, double tmp_numbers[], int tmp_indices[]){

    int mid,k;

    if(low<high){
         mid=(low+high)/2;
         partition(numbers,indices,low,mid,     tmp_numbers, tmp_indices);
         partition(numbers,indices,mid+1, high, tmp_numbers, tmp_indices);
         merge(numbers,indices,low,mid,high,tmp_numbers, tmp_indices);
    }
}

void merge(double numbers[], int indices[], int low,int mid,int high, double tmp_numbers[], int tmp_indices[]){

    int i,m,k,l;

    l=low;
    i=low;
    m=mid+1;

    while((l<=mid)&&(m<=high)){

         if(numbers[l] <= numbers[m]){
             tmp_numbers[i] = numbers[l];
			 tmp_indices[i] = indices[l];
             l++;
         }else{
             tmp_numbers[i] = numbers[m];
			 tmp_indices[i] = indices[m];
             m++;
         }
         i++;
    }

    if(l>mid){
         for(k=m; k <= high; k++){
             tmp_numbers[i] = numbers[k];
			 tmp_indices[i] = indices[k];
             i++;
         }
    }else{
         for(k=l; k <= mid; k++){
             tmp_numbers[i] = numbers[k];
			 tmp_indices[i] = indices[k];
             i++;
         }
    }
   
    for(k=low;k<=high;k++){
         numbers[k] = tmp_numbers[k];
		 indices[k] = tmp_indices[k];
    }
}








double Max(double a, double b){
    if (a < b){
        return b;
    }
    return a;      
}

double Min(double a, double b){
    if (a < b){
        return a;
    }
    return b;      
}
