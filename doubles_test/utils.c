#include <stdio.h>
#include "utils.h"

void init_mat(double **u_ptr, double *u_vals, int n, int m){
    int i;
    for(i=0;i<n;i++)
        u_ptr[i]=&u_vals[i*m];
}

void red_rows(double *A, double *b, int n, int m){
    int i,j;
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
            b[i] += A[j+i*m];
}

double vec_reduce(double *vec, int n){
    double sum = 0.0;
    int i;
    for(i=0;i<n;i++)
        sum += vec[i];

    return sum;
}

double sse(double *a, double *b, int n){
    int i;
    double sse = 0.0;
    for(i=0;i<n;i++)
        sse += (a[i]-b[i])*(a[i]-b[i]);
    return sse;
}
