#include <stdio.h>
#include <math.h>
#include "utils.h"

void init_mat(float **u_ptr, float *u_vals, int n, int m){
	int i;
	for(i=0;i<n;i++)
		u_ptr[i]=&u_vals[i*m];
}

void red_rows(float **A, float *b, int n, int m){
	int i,j;
	for(i=0;i<n;i++)
		for(j=0;j<m;j++)
			b[i] += fabs(A[i][j]);
}

float vec_reduce(float *vec, int n){
	float sum = 0.0;
	int i;
	for(i=0;i<n;i++)
		sum += vec[i];

	return sum;
}

float sse(float *a, float *b, int n){
    int i;
    float sse = 0.0;
    for(i=0;i<n;i++)
        sse += (a[i]-b[i])*(a[i]-b[i]);
    return sse;
}
