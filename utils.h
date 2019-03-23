#ifndef _UTILS_H
#define _UTILS_H

typedef struct {
	float calc_CPU;
	float alloc_GPU;
	float transf_GPU;
	float calc_GPU;
	float calc_avg;
	float transf_RAM;
} Tau;

void init_mat(float **u_ptr, float *u_vals, int m, int n);
void red_rows(float *u, float *t, int n, int m);
float vec_reduce(float *vec, int n);
float sse(float *a, float *b, int n);

#endif
