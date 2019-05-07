#ifndef _FINDIFF_GPU_H
#ifdef __cplusplus

extern "C" {
    void fdiff_gpu(float *u_vals, float *temps, int n, int m, int p, int block_size_Y, Tau* tau, int mallocPitch, int red);
    void fdiff_gpu_slow(float* u_vals, float* temps, int n, int m, int p, int block_size, Tau* tau, int red);
}

#endif

#define _FINDIFF_GPU_H

__device__ void calc_iterate(float *unew, float *uold, 
                        int n, int m, int idx, int idy, int ind);
__device__ void glob_shared_cpy(float *u_glob, float *unew, float *uold, 
                        int n, int m, int idx, int idy, int ind);
__device__ void shared_glob_cpy(float *uglob, float *unew, 
                        int n, int m, int idx, int idy, int ind);

__global__ void red_rows(float* u_glob, float* u_glob_out, int pitch, int n, int m, int m_tot);
__global__ void red_rows_glob(float* u_glob, float* u_glob_out, int pitch, int n, int m);
__global__ void iterate_gpu(float *u_glob, int n, int m);
__global__ void iterate_gpu_(float* unew_glob, float* uold_glob, int n, int m);

#endif
