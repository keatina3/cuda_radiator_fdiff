#include "utils.h"
#include "findiff_gpu.h"
#include "findiff.h"

// applies single propagation to shared memory arrays //
__device__ void calc_iterate(float *unew, float *uold, int n, int m, int idx, int idy, int ind){
	if(1<idy && idy<m){
		if(idx<n){
			unew[ind] = (1.9*uold[ind-2] + 1.5*uold[ind-1] +
                            uold[ind] + 0.5*uold[ind+1] + 0.1*uold[ind+2]);
			unew[ind] /= (float)(5.0);
		}
	}
}

// reading data from global to shared memory //
__device__ void glob_shared_cpy(float *u_glob, float *unew, float *uold, int pitch, int n, int m, int idx, int idy, int ind){
    if(idy<m && idx<n){
		// copying left halo values //
		if(threadIdx.y==0 && 0<blockIdx.y){
			unew[ind-2] = u_glob[(idy-2) + idx*pitch];
			unew[ind-1] = u_glob[(idy-1) + idx*pitch];	
			uold[ind-2] = unew[ind-2];
			uold[ind-1] = unew[ind-1];
		}
		// copying normal grid values //
		unew[ind] = u_glob[idy+idx*pitch];
		uold[ind] = unew[ind];
		//copying right halo values //
		if(threadIdx.y==(blockDim.y-1) || idy==(m-1)){
			unew[ind+1] = u_glob[(idy+1)%m + idx*pitch];
			unew[ind+2] = u_glob[(idy+2)%m + idx*pitch];
            uold[ind+1] = unew[ind+1];
           	uold[ind+2] = unew[ind+2];
        }
    }
}

// sending data from shared back to global memory //
__device__ void shared_glob_cpy(float *u_glob, float *unew, int pitch, int n, int m, int idx, int idy, int ind){
    if(1<idy && idy<m)
		if(idx<n)
		    u_glob[idy+idx*pitch] = unew[ind];
}

// kernel for shared-mem & cudaMallocPitch approach //
__global__ void iterate_gpu(float *u_glob, int pitch, int n, int m){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
	int ind = threadIdx.y + 2;
    float *uold, *unew;
   	extern __shared__ float s[];
    
	// can only dynamically allocate one shared mem array //
	// dividing array into two smaller arrays //
    unew = &(s[0]);
    uold = &(s[blockDim.y+4]);
    
	// initialising shared memory //
    glob_shared_cpy(u_glob, unew, uold, pitch, n, m, idx, idy, ind);     
    
   	// iterating and updating unew //
	calc_iterate(unew, uold, n, m, idx, idy, ind);
    
    // sending vals back to global mem //
    shared_glob_cpy(u_glob, unew, pitch, n, m, idx, idy, ind);
}

// kernel for global-mem/non-optimised approach //
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

// optimised reduce kernel //
__global__ void red_rows(float* u_glob, float* u_glob_out, int pitch, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    int ind = threadIdx.y;
    extern __shared__ float tmp[];
    int i, disp;
	
	// reading from global mem to shared //
    if(idy<m && idx<n)
        tmp[ind] = u_glob[idy+idx*pitch];
        
    // for loop below performs binary reduction on each block //
    
    // disp variable to check if threadId @ end of block > m //
    // if disp > m, then length of last array to be reduced = threadIdx @ end //
    disp = (1+blockIdx.y)*blockDim.y;
    i = (disp > m) ? (blockDim.y - (disp-m)):blockDim.y;
    // perform binary reduction //
    for( ; i>1; i>>=1){
        if(ind<(i/2)){
            tmp[ind] += tmp[ind+(i/2)];	// sum value ind and value (1/2)*size.tmp away //
            if(ind==0 && i%2!=0)		// if odd, add last value also to arr[0] //
                tmp[ind] += tmp[ind+i-1]; 
        }
         __syncthreads();				// sync each block
    }
    if(ind==0)							// write the #blocks values back to global mem //
        u_glob_out[blockIdx.y + idx*pitch] = tmp[0];
}

// unoptimised global memory, binary reduction //
__global__ void red_rows_glob(float* u_glob, float* u_glob_out, int pitch, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;

    if(idy < (m/2) && idx < n){
        u_glob_out[idy + idx*pitch] += u_glob[idy + (m/2) + idx*pitch];
        if(m%2!=0 && idy==0)
            u_glob_out[idy + idx*pitch] += u_glob[idy + (m-1) + idx*pitch];
    }
}

// C function for optimised approach //
extern "C" {
void fdiff_gpu(float *u_vals, float *temps, int n, int m, int p, int block_size_Y, Tau* tau, int mallocPitch, int red){
    float *u_glob;
    size_t u_glob_size;
    int i, pitch, m_tmp;
	cudaEvent_t start, finish;

    cudaEventCreate(&start);
    cudaEventCreate(&finish);
    
    if(!mallocPitch){		// if mallocPitch not flagged, use cudaMalloc
        cudaEventRecord(start, 0);
        cudaMalloc( (void**)&u_glob, n*m*sizeof(float));
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->alloc_GPU, start, finish);
        
        cudaEventRecord(start,0);
        cudaMemcpy(u_glob, u_vals, n*m*sizeof(float), cudaMemcpyHostToDevice);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->transf_GPU, start, finish);
        
        pitch = m;
    } else {
    	// allocating global mem //
        cudaEventRecord(start, 0);
        cudaMallocPitch( (void**)&u_glob, &u_glob_size, (size_t)(m*sizeof(float)), n);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->alloc_GPU, start, finish);
        
        // transfer RAM->DRAM //
        cudaEventRecord(start,0);
        cudaMemcpy2D(u_glob, u_glob_size, u_vals, m*sizeof(float), m*sizeof(float), n, cudaMemcpyHostToDevice);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->transf_GPU, start, finish);
        
        // pitch = offset of values in array, since not contiguous //
        pitch = (int)u_glob_size/sizeof(float);
    }

    dim3 dimBlock(1, block_size_Y);
    dim3 dimGrid((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));
   
   	// calling on kernel for p propagations //
    cudaEventRecord(start, 0);
    for(i=0;i<p;i++)
	    iterate_gpu<<<dimGrid,dimBlock,2*(block_size_Y+4)*sizeof(float)>>>(u_glob, pitch, n, m);
    
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->calc_GPU, start, finish);
	
	// copying back to RAM //
    if(!mallocPitch){
        cudaEventRecord(start, 0);
        cudaMemcpy(u_vals, u_glob, n*m*sizeof(float), cudaMemcpyDeviceToHost);
    } else {
        cudaEventRecord(start, 0);
        cudaMemcpy2D(u_vals, m*sizeof(float), u_glob, u_glob_size, m*sizeof(float), n, cudaMemcpyDeviceToHost);
    }
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->transf_RAM, start, finish);
    
    // if reduce function flagged //
    if(red){  
        m_tmp = m;
        // reduce kernel //
        cudaEventRecord(start, 0);
        while(m_tmp > 1){
            red_rows<<<dimGrid,dimBlock,dimBlock.y*sizeof(float)>>>(u_glob, u_glob, pitch, n, m_tmp);
            m_tmp = (m_tmp/dimBlock.y)+(!(m_tmp%dimBlock.y)?0:1);	// reduced by order blockDim //
        }
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->calc_avgGPU, start, finish);
        
        // copying back to RAM //
        if(!mallocPitch){
            for(i=0;i<n;i++)
                cudaMemcpy(&temps[i], &u_glob[i*m], sizeof(float), cudaMemcpyDeviceToHost);
        } else {
            cudaMemcpy2D(temps, sizeof(float), &u_glob[0], u_glob_size, sizeof(float), n, cudaMemcpyDeviceToHost);
        }
    }

    cudaFree(u_glob);
}

// C code for glob mem approach //
void fdiff_gpu_glob(float* u_vals, float* temps, int n, int m, int p, int block_size, Tau* tau, int red){
	float *uold_glob, *unew_glob, *tmp;
    int i, m_tmp;
	cudaEvent_t start, finish;

    cudaEventCreate(&start);
    cudaEventCreate(&finish);
    
    cudaEventRecord(start, 0);
    cudaMalloc( (void**)&unew_glob, n*m*sizeof(float));
    cudaMalloc( (void**)&uold_glob, n*m*sizeof(float));
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->alloc_GPU, start, finish);

    cudaEventRecord(start, 0);
    cudaMemcpy(unew_glob, u_vals, n*m*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(uold_glob, u_vals, n*m*sizeof(float), cudaMemcpyHostToDevice);
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->transf_GPU, start, finish);
    
    dim3 dimBlock(1, block_size);
    dim3 dimGrid ((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));

    cudaEventRecord(start, 0);
    for(i=0;i<p;i++){
        if(i%2==0)
	        iterate_gpu_slow<<<dimGrid,dimBlock>>>(unew_glob, uold_glob, n, m);
        else
	        iterate_gpu_slow<<<dimGrid,dimBlock>>>(uold_glob, unew_glob, n, m);
    }
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->calc_GPU, start, finish);
    
    cudaEventRecord(start, 0);
    if(p%2==0)
        tmp = uold_glob;
    else
        tmp = unew_glob;
    cudaMemcpy(u_vals, tmp, n*m*sizeof(float), cudaMemcpyDeviceToHost);
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->transf_RAM, start, finish);
    
    m_tmp = m;
    if(red){    
        cudaEventRecord(start, 0);
        for( ; m_tmp>1; m_tmp>>=1)
            red_rows_glob<<<dimGrid,dimBlock>>>(tmp, tmp, m, n, m_tmp);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->calc_avgGPU, start, finish);
        for(i=0;i<n;i++)
            cudaMemcpy(&temps[i], &tmp[i*m], sizeof(float), cudaMemcpyDeviceToHost);
    }
    
    cudaFree(unew_glob); cudaFree(uold_glob);
}
}
