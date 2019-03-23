#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/time.h>
#include "utils.h"
#include "findiff.h"
#include "radiator.h"

extern void fdiff_gpu(float* u_vals, float* temps, int n, int m, int p, int block_size_Y, Tau* tau, int mallocPitch, int red);
extern void fdiff_gpu_glob(float* u_vals, float* temps, int n, int m, int p, int block_size, Tau* tau, int red);

int main(int argc, char **argv){
	int m = 32, n = 32, p = 10, block_size = 16;
	int t = 0, w = 0, a = 0, option = 0, serial = 0, mallocPitch = 1, glob = 0;
	float **uold, **unew;
	float *u, *uold_vals, *unew_vals, *temps, *uGPU, *tempsGPU, *uglob_GPU, *tempsGPU2;
	struct timeval start, end;
	Tau tau1, tau2;

	while((option=getopt(argc,argv,"n:m:p:b:tawsdg"))!=-1){
		switch(option){
			case 'n': n = atoi(optarg);
				break;
			case 'm': m = atoi(optarg);
				break;
			case 'p': p = atoi(optarg);
				break;
			case 'b': block_size = atoi(optarg);
				break;
			case 't': t = 1;
				break;
			case 'a': a = 1;
				break;
			case 'w': w = 1;	// to write results to file //
				break;
			case 's': serial = 1;
				break;
			case 'd': mallocPitch = 0;
				break;
			case 'g': glob = 1;
				break;
			default:
				printf("Incorrect options entered!\n");
				return 1;
		}
	}	
	if(argc != optind){
		printf("Too many arguments provided, exiting!\n");
		return 1;
	}
	
	uold = (float**)malloc(n*sizeof(float*));
	unew = (float**)malloc(n*sizeof(float*));
	uold_vals = (float*)calloc(n*m,sizeof(float));
	unew_vals = (float*)calloc(n*m,sizeof(float));
	uGPU = (float*)calloc(n*m,sizeof(float));
    uglob_GPU = (float*)calloc(n*m,sizeof(float));
    temps = (float*)calloc(n, sizeof(float));	
    tempsGPU = (float*)calloc(n, sizeof(float));	
    tempsGPU2 = (float*)calloc(n, sizeof(float));	
	
    init_mat(uold, uold_vals, n, m);	
	init_mat(unew, unew_vals, n, m);	
	
	apply_bounds(unew, n, 2, &fleft, &fright, &fzero, &fzero);
	apply_bounds(uold, n, 2, &fleft, &fright, &fzero, &fzero);
    
    memcpy(uGPU, unew_vals, n*m*sizeof(float));
    memcpy(uglob_GPU, unew_vals, n*m*sizeof(float));
    
    if(serial){
	    iterate(uold, unew, p, n, m);
        if(p%2==0)
            u = uold_vals;
        else
            u = unew_vals;
	    printf("\n\n\n");
        print_grid(u, n, m);
	    if(a)
            red_rows(u, temps, n, m);
    }
    if(glob){
        fdiff_gpu_glob(uglob_GPU, tempsGPU2, n, m, p, block_size, &tau2, a);
	    printf("\n\n\n");
        print_grid(uglob_GPU, n, m);
    }

    fdiff_gpu(uGPU, tempsGPU, n, m, p, block_size, &tau1, mallocPitch, a);
	printf("\n\n\n");
    print_grid(uGPU, n, m);
    
    float SSE = sse(temps, tempsGPU, n);
    int i;
    for(i=0;i<n;i++){
        printf("temps[%d] = %f, tempsGPU[%d] = %f, tempsGPUglob[%d]=%f\n",i,temps[i],i,tempsGPU[i],i,tempsGPU[i]);
    } 
    printf("SSE of temps = %f\n",SSE);

	free(uold); free(unew);
    free(temps); free(tempsGPU); free(tempsGPU2);
    free(uGPU); free(uglob_GPU);
	free(uold_vals); free(unew_vals);

	return 0;
}
