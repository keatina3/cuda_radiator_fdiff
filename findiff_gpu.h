#ifndef _FINDIFF_GPU_H
#ifdef __cplusplus

extern "C" {
	void fd_iterate_gpu(float* u_vals, int block_size, int n, int m, int p, float* tau);
}

#endif

#define _FINDIFF_GPU_H


#endif
