#ifndef _UTILS_H
#define _UTILS_H

void init_mat(float **u_ptr, float *u_vals, int m, int n);
void sum_rows(float **u, float *t, int n, int m);
//void sum_cols(float **u, float *t, int n, int m);
float vec_reduce(float *vec, int n);

#endif
