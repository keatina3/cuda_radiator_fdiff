#include <stdlib.h>
#include "utils.h"
#include "findiff_gpu.h"

__device__ void calc_iterate(float *u_vals, float *unew, float *uold, int n, int m, int block_size_X){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
	int ind = threadIdx.x + 2;

	if(idy<m){
		if(idx<n){
			unew[ind] = 1.9*uold[ind-2] + 1.5*uold[ind-1] + uold[ind] + 0.5*uold[ind+1] + 0.1*uold[ind+2];
			unew[ind] /= (float)(5.0);
		}
	}
	__syncthreads();
	
	if(threadIdx.x==0){
		unew[ind-2] = u_vals[(idx+idy*n)-2];
		unew[ind-1] = u_vals[(idx+idy*n)-1];	
		uold[ind-2] = unew[ind-2];
		uold[ind-1] = unew[ind-1];
	}
	if(threadIdx.x==(block_size_X-1)){
		unew[ind+1] = u_vals[(idx+idy*n)+1];
		unew[ind+2] = u_vals[(idx+idy*n)+2];
		uold[ind+1] = unew[ind+1];
		uold[ind+2] = unew[ind+2];
	}
}

__global__ void iterate_gpu(float *u_vals, int n, int m, int p){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
	int ind = threadIdx.x + 2;
	int i;

	__shared__ float unew[32+4];
	__shared__ float uold[32+4];
	
	// initialising shared memory //
	if(idy<m){
		if(idx<n){
			if(threadIdx.x==0){
				unew[ind-2] = u_vals[(idx+idy*n)-2];		// check this index and below //
				unew[ind-1] = u_vals[(idx+idy*n)-1];	
				uold[ind-2] = unew[ind-2];					// need to deal with cyclic factor //
				uold[ind-1] = unew[ind-1];
			}

			unew[ind] = u_vals[idx+idy*n];
			uold[ind] = unew[threadIdx.x+2];
			
			if(threadIdx.x==(32-1)){
				unew[ind+1] = u_vals[(idx+idy*n)+1];
				unew[ind+2] = u_vals[(idx+idy*n)+2];
				uold[ind+1] = unew[ind+1];
				uold[ind+2] = unew[ind+2];
			}
		}
	}
	
	for(i=0;i<p;i++){
		if(i%2==0)
			calc_iterate(u_vals, unew, uold, n, m, 32);
		else
			calc_iterate(u_vals, uold, unew, n, m, 32);
	}
}

extern void fdiff_gpu(float* u_vals, int block_size_X, int block_size_Y, int n, int m, int p, Tau* tau){
	float *u_glob;
	size_t u_glob_size;

	cudaMallocPitch( (void**)&u_glob, &u_glob_size, (size_t)(n*sizeof(float)), m);
	cudaMemcpy2D(u_glob, u_glob_size, u_vals, n*sizeof(float), n*sizeof(float), m, cudaMemcpyHostToDevice);
	
	dim3 dimBlock(block_size_X, block_size_Y);
	dim3 dimGrid ((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));

	iterate_gpu<<<dimGrid, dimBlock>>>(u_vals, n, m, p);
	
	cudaFree(u_glob);
}
