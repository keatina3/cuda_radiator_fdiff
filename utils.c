#include <stdio.h>
#include <math.h>
#include "utils.h"

void init_mat(float **u_ptr, float *u_vals, int n, int m){
	int i;
	for(i=0;i<n;i++)
		u_ptr[i]=&u_vals[i*m];
}

void sum_rows(float **A, float *b, int n, int m){
	int i,j;
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			b[i] += fabs(A[i][j]);
}

/*
// sum all values on each column //
void sum_cols(float **A, float *b, int n, int m){
	int i,j;
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			b[j] += f_abs(A[i][j]);
}
*/
// reduce vector to its sum //
float vec_reduce(float *vec, int n){
	float sum = 0.0;
	int i;
	for(i=0;i<n;i++)
		sum += vec[i];

	return sum;
}
