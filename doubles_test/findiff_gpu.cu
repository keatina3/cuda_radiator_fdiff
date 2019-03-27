#include "utils.h"
#include "findiff_gpu.h"
#include "findiff.h"

__device__ void calc_iterate(double *unew, double *uold, int n, int m, int idx, int idy, int ind){
    if(1<idy && idy<m){
        if(idx<n){
            unew[ind] = (1.9*uold[ind-2] + 1.5*uold[ind-1] +
                    uold[ind] + 0.5*uold[ind+1] + 0.1*uold[ind+2]);
            unew[ind] /= (double)(5.0);
        }
    }
}

__device__ void glob_shared_cpy(double *u_glob, double *unew, double *uold, int pitch, int n, int m, int idx, int idy, int ind){
        
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

__device__ void shared_glob_cpy(double *u_glob, double *unew, int pitch, int n, int m, int idx, int idy, int ind){
    if(1<idy && idy<m)
        if(idx<n)
            u_glob[idy+idx*pitch] = unew[ind];
}

__global__ void iterate_gpu(double *u_glob, int pitch, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    int ind = threadIdx.y + 2;
    double *uold, *unew;
    extern __shared__ double s[];

    unew = &(s[0]);
    uold = &(s[blockDim.y+4]);

    // initialising shared memory //
    glob_shared_cpy(u_glob, unew, uold, pitch, n, m, idx, idy, ind);     

    // iterating and updating unew //
    calc_iterate(unew, uold, n, m, idx, idy, ind);

    // sending vals back to global mem //
    shared_glob_cpy(u_glob, unew, pitch, n, m, idx, idy, ind);
}

__global__ void iterate_gpu_slow(double* unew_glob, double* uold_glob, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;

    if(1<idy && idy<m){
        unew_glob[idy+idx*m] = (1.9*uold_glob[(idy+idx*m)-2] + 1.5*uold_glob[(idy+idx*m)-1] +
                uold_glob[idy+idx*m] + 0.5*uold_glob[(idy+1)%m+idx*m] 
                    + 0.1*uold_glob[(idy+2)%m+idx*m]);
        unew_glob[idy+idx*m] /= (double)(5.0);
        }
    }
}

__global__ void red_rows(double* u_glob, double* u_glob_out, int pitch, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    int ind = threadIdx.y;
    extern __shared__ double tmp[];
    int i, disp;

    if(idy<m && idx<n)
        tmp[ind] = u_glob[idy+idx*pitch];

    disp = (1+blockIdx.y)*blockDim.y;
    i = (disp > m) ? (blockDim.y - (disp-m)):blockDim.y;
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

__global__ void red_rows_glob(double* u_glob, double* u_glob_out, int pitch, int n, int m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;

    if(idy < (m/2) && idx < n){
        u_glob_out[idy + idx*pitch] += u_glob[idy + (m/2) + idx*pitch];
        if(m%2!=0 && idy==0)
            u_glob_out[idy + idx*pitch] += u_glob[idy + (m-1) + idx*pitch];
    }
}

extern "C" {
void fdiff_gpu(double *u_vals, double *temps, int n, int m, int p, int block_size_Y, Tau* tau, int mallocPitch, int red){
    double *u_glob;
    size_t u_glob_size;
    int i, pitch, m_tmp;
    cudaEvent_t start, finish;

    cudaEventCreate(&start);
    cudaEventCreate(&finish);

    if(!mallocPitch){
        cudaEventRecord(start, 0);
        cudaMalloc( (void**)&u_glob, n*m*sizeof(double));
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->alloc_GPU, start, finish);

        cudaEventRecord(start,0);
        cudaMemcpy(u_glob, u_vals, n*m*sizeof(double), cudaMemcpyHostToDevice);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->transf_GPU, start, finish);

        pitch = m;
    } else {
        cudaEventRecord(start, 0);
        cudaMallocPitch( (void**)&u_glob, &u_glob_size, (size_t)(m*sizeof(double)), n);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->alloc_GPU, start, finish);

        cudaEventRecord(start,0);
        cudaMemcpy2D(u_glob, u_glob_size, u_vals, m*sizeof(double), m*sizeof(double), n, cudaMemcpyHostToDevice);
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->transf_GPU, start, finish);

        pitch = (int)u_glob_size/sizeof(double);
    }

    dim3 dimBlock(1, block_size_Y);
    dim3 dimGrid((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));

    cudaEventRecord(start, 0);
    for(i=0;i<p;i++)
        iterate_gpu<<<dimGrid,dimBlock,2*(block_size_Y+4)*sizeof(double)>>>(u_glob, pitch, n, m);

    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->calc_GPU, start, finish);

    if(!mallocPitch){
        cudaEventRecord(start, 0);
        cudaMemcpy(u_vals, u_glob, n*m*sizeof(double), cudaMemcpyDeviceToHost);
    } else {
        cudaEventRecord(start, 0);
        cudaMemcpy2D(u_vals, m*sizeof(double), u_glob, u_glob_size, m*sizeof(double), n, cudaMemcpyDeviceToHost);
    }
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->transf_RAM, start, finish);

    if(red){  
        m_tmp = m;
        cudaEventRecord(start, 0);
        while(m_tmp > 1){
            red_rows<<<dimGrid,dimBlock,dimBlock.y*sizeof(double)>>>(u_glob, u_glob, pitch, n, m_tmp);
            m_tmp = (m_tmp/dimBlock.y)+(!(m_tmp%dimBlock.y)?0:1);
        }
        cudaEventRecord(finish, 0);
        cudaEventSynchronize(finish);
        cudaEventElapsedTime(&tau->calc_avgGPU, start, finish);

        if(!mallocPitch){
            for(i=0;i<n;i++)
                cudaMemcpy(&temps[i], &u_glob[i*m], sizeof(double), cudaMemcpyDeviceToHost);
        } else {
            cudaMemcpy2D(temps, sizeof(double), &u_glob[0], u_glob_size, sizeof(double), n, cudaMemcpyDeviceToHost);
        }
    }

    cudaFree(u_glob);
}

void fdiff_gpu_glob(double* u_vals, double* temps, int n, int m, int p, int block_size, Tau* tau, int red){
    double *uold_glob, *unew_glob, *tmp;
    int i, m_tmp;
    cudaEvent_t start, finish;

    cudaEventCreate(&start);
    cudaEventCreate(&finish);

    cudaEventRecord(start, 0);
    cudaMalloc( (void**)&unew_glob, n*m*sizeof(double));
    cudaMalloc( (void**)&uold_glob, n*m*sizeof(double));
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->alloc_GPU, start, finish);

    cudaEventRecord(start, 0);
    cudaMemcpy(unew_glob, u_vals, n*m*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(uold_glob, u_vals, n*m*sizeof(double), cudaMemcpyHostToDevice);
    cudaEventRecord(finish, 0);
    cudaEventSynchronize(finish);
    cudaEventElapsedTime(&tau->transf_GPU, start, finish);

    dim3 dimBlock(1, block_size);
    dim3 dimGrid((n/dimBlock.x)+(!(n%dimBlock.x)?0:1), (m/dimBlock.y)+(!(m%dimBlock.y)?0:1));

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
    cudaMemcpy(u_vals, tmp, n*m*sizeof(double), cudaMemcpyDeviceToHost);
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
            cudaMemcpy(&temps[i], &tmp[i*m], sizeof(double), cudaMemcpyDeviceToHost);
    }

    cudaFree(unew_glob); cudaFree(uold_glob);
}
}
