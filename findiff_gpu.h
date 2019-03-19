#ifndef _FINDIFF_GPU_H
#ifdef __cplusplus

extern "C" {
	void fdiff_gpu(float* u_vals, int block_size_X, int block_size_Y, int n, int m, int p, Tau* tau);
}

#endif

#define _FINDIFF_GPU_H

__device__ void calc_iterate(float *u_glob, float *unew, float *uold, int n, int m, int block_size_X);
__global__ void iterate_gpu(float *u_vals, int n, int m, int p);

#endif
