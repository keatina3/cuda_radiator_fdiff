#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
//#include <string.h>
//#include <sys/time.h>
#include "utils.h"
#include "findiff.h"
#include "radiator.h"

extern void fd_iterate_gpu(float* u_vals, int block_size_X, int block_size_Y, int n, int m, int p, Tau* tau);

int main(int argc, char **argv){
	int m = 32, n = 32, p = 10;
	int t = 0, w = 0, a = 0, option = 0;
	float **uold, **unew;
	float *uold_vals, *unew_vals, *uoldGPU, *unewGPU;
	struct timeval start, end;
	Tau tau;

	while((option=getopt(argc,argv,"n:m:p:taw"))!=-1){
		switch(option){
			case 'n': n = atoi(optarg);
				break;
			case 'm': m = atoi(optarg);
				break;
			case 'p': p = atoi(optarg);
				break;
			case 't': t = 1;
				break;
			case 'a': a = 1;
				break;
			case 'w': w = 1;	// to write results to file //
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
	
	init_mat(uold, uold_vals, n, m);	
	init_mat(unew, uold_vals, n, m);	
	
	apply_bounds(unew, n, 2, &fleft, &fright, &fzero, &fzero);
	apply_bounds(uold, n, 2, &fleft, &fright, &fzero, &fzero);

	iterate(uold, unew, p, n, m);

	get_grid_avg(unew, n, m);

	print_grid(unew, n, m);

	free(uold); free(unew);
	free(uold_vals); free(unew_vals);

	return 0;
}
