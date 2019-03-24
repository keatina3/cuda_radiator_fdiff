#ifndef _UTILS_H
#define _UTILS_H

#define RED  "\x1B[91m"     // for colour setting in print
#define GREEN  "\x1B[92m"     // for colour setting in print
#define BLUE  "\x1B[94m"     // for colour setting in print
#define MAG  "\x1B[95m"     // for colour setting in print
#define RESET "\x1B[0m"

typedef struct {
	float calc_CPU;
	float alloc_GPU;
	float transf_GPU;
	float calc_GPU;
	float calc_avgGPU;
    float calc_avgCPU;
	float transf_RAM;
    float tot_GPU;
    float tot_CPU;
} Tau;

void init_mat(float **u_ptr, float *u_vals, int m, int n);
void red_rows(float *u, float *t, int n, int m);
float vec_reduce(float *vec, int n);
float sse(float *a, float *b, int n);

#endif
