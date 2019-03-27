#ifndef _UTILS_H
#define _UTILS_H

#define RED  "\x1B[91m"
#define GREEN  "\x1B[92m"
#define BLUE  "\x1B[94m"
#define MAG  "\x1B[95m"
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

void init_mat(double **u_ptr, double *u_vals, int m, int n);
void red_rows(double *u, double *t, int n, int m);
double vec_reduce(double *vec, int n);
double sse(double *a, double *b, int n);

#endif
