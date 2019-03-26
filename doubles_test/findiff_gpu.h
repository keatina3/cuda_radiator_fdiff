#ifndef _FINDIFF_GPU_H
#ifdef __cplusplus

extern "C" {
    void fdiff_gpu(double *u_vals, double *temps, int n, int m, int p, int block_size_Y, Tau* tau, int mallocPitch, int red);
    void fdiff_gpu_slow(double* u_vals, double* temps, int n, int m, int p, int block_size, Tau* tau, int red);
}

#endif

#define _FINDIFF_GPU_H

__device__ void calc_iterate(double *unew, double *uold, 
                        int n, int m, int idx, int idy, int ind);
__device__ void glob_shared_cpy(double *u_glob, double *unew, double *uold, 
                        int n, int m, int idx, int idy, int ind);
__device__ void shared_glob_cpy(double *uglob, double *unew, 
                        int n, int m, int idx, int idy, int ind);

__global__ void red_rows(double* u_glob, double* u_glob_out, int pitch, int n, int m);
__global__ void red_rows_glob(double* u_glob, double* u_glob_out, int pitch, int n, int m);
__global__ void iterate_gpu(double *u_glob, int n, int m);
__global__ void iterate_gpu_(double* unew_glob, double* uold_glob, int n, int m);

#endif
