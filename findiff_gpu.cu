#include "utils.h"
#include "findiff_gpu.h"

__device__ void calc_iterate(float *unew, float *uold, int n, int m, int idx, int idy, int ind){
	if(1<idy && idy<m){
		if(idx<n){
			unew[ind] = (1.9*uold[ind-2] + 1.5*uold[ind-1] +
                            uold[ind] + 0.5*uold[ind+1] + 0.1*uold[ind+2]);
			unew[ind] /= (float)(5.0);
		}
	}
	__syncthreads();
}

__device__ void glob_shared_cpy(float *u_glob, float *unew, float *uold, int pitch, int n, int m, int idx, int idy, int ind){
	// READING DATA FROM GLOBAL MEMORY TO SHARED //
    if(idy<m && idx<n){
		if(threadIdx.y==0 && 0<blockIdx.y){
			unew[ind-2] = u_glob[(idy-2) + idx*pitch];
			unew[ind-1] = u_glob[(idy-1) + idx*pitch];	
			uold[ind-2] = unew[ind-2];
			uold[ind-1] = unew[ind-1];
		}

		unew[ind] = u_glob[idy+idx*pitch];
		uold[ind] = unew[ind];
			
		if(threadIdx.y==(blockDim.y-1) || idy==(m-1)){
			unew[ind+1] = u_glob[(idy+1)%m + idx*pitch];
			unew[ind+2] = u_glob[(idy+2)%m + idx*pitch];
            uold[ind+1] = unew[ind+1];
            uold[ind+2] = unew[ind+2];
        }
    }
}

__device__ void shared_glob_cpy(float *u_glob, float *unew, int pitch, int n, int m, int idx, int idy, int ind){
    if(1<idy && idy<m)
		if(idx<n)
		    u_glob[idy+idx*pitch] = unew[ind];
}

__global__ void iterate_gpu(float *u_glob, int pitch, int n, int m){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
	int ind = threadIdx.y + 2;
    float *uold, *unew;
    extern __shared__ float s[];
    
    unew = &(s[0]);
    uold = &(s[blockDim.y+4]);
    
	// initialising shared memory //
    glob_shared_cpy(u_glob, unew, uold, pitch, n, m, idx, idy, ind);     
    
    // iterating and updating unew //
	calc_iterate(unew, uold, n, m, idx, idy, ind);
    
    // sending vals back to global mem //
    shared_glob_cpy(u_glob, unew, pitch, n, m, idx, idy, ind);
}

__global__ void iterate_gpu_slow(float* unew_glob, float* uold_glob, int n, int m){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
    if(1<idy && idy<m){
		if(idx<n){
			unew_glob[idy+idx*m] = (1.9*uold_glob[(idy+idx*m)-2] + 1.5*uold_glob[(idy+idx*m)-1] +
                        uold_glob[idy+idx*m] + 0.5*uold_glob[(idy+1)%m+idx*m] 
                            + 0.1*uold_glob[(idy+2)%m+idx*m]);
			unew_glob[idy+idx*m] /= (float)(5.0);
		}
	}
}

__global__ void red_rows(float* u_glob, float* u_glob_out, int pitch, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    int ind = threadIdx.y;
    extern __shared__ float tmp[];
    int i, disp;

    if(idy<m && idx<n)
        tmp[ind] = u_glob[idy+idx*pitch];         

    disp = (1+blockIdx.y)*blockDim.y;
    i = disp > m ? (blockDim.y - disp%m):blockDim.y;
    for( ; i>1; i>>=1){
        if(ind<(i/2)){
            tmp[ind] += tmp[ind+(i/2)];
            if(ind==0 && i%2!=0)
                tmp[ind] += tmp[ind+i-1]; 
        }
        __syncthreads();
    }

    if(ind==0)
        u_glob_out[blockIdx.y + idx*pitch] = tmp[0];
}

extern "C" {
void fdiff_gpu(float *u_vals, float *temps, int n, int m, int p, int block_size_Y, Tau* tau){
    float *u_glob;
    size_t u_glob_size;
    int i, pitch, m_tmp;
	
    //cudaMalloc( (void**)&u_glob, n*m*sizeof(float));
    //cudaMemcpy(u_glob, u_vals, n*m*sizeof(float), cudaMemcpyHostToDevice);
    cudaMallocPitch( (void**)&u_glob, &u_glob_size, (size_t)(m*sizeof(float)), n);
    cudaMemcpy2D(u_glob, u_glob_size, u_vals, m*sizeof(float), m*sizeof(float), n, cudaMemcpyHostToDevice);
    pitch = (int)u_glob_size/sizeof(float);
    // pitch = m;

    dim3 dimBlock(1, block_size_Y);
    dim3 dimGrid((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));
   
    for(i=0;i<p;i++)
	    iterate_gpu<<<dimGrid,dimBlock,2*(block_size_Y+4)*sizeof(float)>>>(u_glob, pitch, n, m);
    
    //cudaMemcpy(u_vals, u_glob, n*m*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy2D(u_vals, m*sizeof(float), u_glob, u_glob_size, m*sizeof(float), n, cudaMemcpyDeviceToHost);
    
    m_tmp = m;
    while(m_tmp > 1){
        red_rows<<<dimGrid,dimBlock,block_size_Y>>>(u_glob, u_glob, pitch, n, m_tmp);
        m_tmp = (m_tmp/dimBlock.y)+(!(m_tmp%dimBlock.y)?0:1);
    }
    cudaMemcpy2D(temps, sizeof(float), u_glob, u_glob_size, sizeof(float), n, cudaMemcpyDeviceToHost);

    cudaFree(u_glob);
}

void fdiff_gpu_slow(float* u_vals, int n, int m, int p, int block_size, Tau* tau){
	float *uold_glob, *unew_glob;
    int i;
    
    cudaMalloc( (void**)&unew_glob, n*m*sizeof(float));
    cudaMalloc( (void**)&uold_glob, n*m*sizeof(float));

    cudaMemcpy(unew_glob, u_vals, n*m*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(uold_glob, u_vals, n*m*sizeof(float), cudaMemcpyHostToDevice);
    
    dim3 dimBlock(block_size, block_size);
    dim3 dimGrid ((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));

    for(i=0;i<p;i++){
        if(i%2==0)
	        iterate_gpu_slow<<<dimGrid,dimBlock>>>(unew_glob, uold_glob, n, m);
        else
	        iterate_gpu_slow<<<dimGrid,dimBlock>>>(uold_glob, unew_glob, n, m);
    }
    
    if(p%2==0)
        cudaMemcpy(u_vals, unew_glob, n*m*sizeof(float), cudaMemcpyDeviceToHost);
    else
        cudaMemcpy(u_vals, uold_glob, n*m*sizeof(float), cudaMemcpyDeviceToHost);
	
    cudaFree(unew_glob); cudaFree(uold_glob);
}
}
