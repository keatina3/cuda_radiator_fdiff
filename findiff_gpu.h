#ifndef _FINDIFF_GPU_H
#ifdef __cplusplus

extern "C" {
    void fdiff_gpu(float* u_vals, int n, int m, int p, Tau* tau);
    void fdiff_gpu_slow(float* u_vals, int n, int m, int p, Tau* tau);
}

#endif

#define _FINDIFF_GPU_H

#define BLOCK_SIZE 4

extern const int block_size_Y = BLOCK_SIZE;
extern const int block_size_X = BLOCK_SIZE;

__device__ void calc_iterate(float *unew, float *uold, 
                        int n, int m, int idx, int idy, int ind);
__device__ void glob_shared_cpy(float *u_glob, float *unew, float *uold, 
                        int n, int m, int idx, int idy, int ind, int block_size_Y);
__device__ void shared_glob_cpy(float *uglob, float *unew, 
                        int n, int m, int idx, int idy, int ind);
__global__ void iterate_gpu(float *u_glob, int n, int m);
__global__ void iterate_gpu_slow(float* unew_glob, float* uold_glob, int n, int m);

#endif
