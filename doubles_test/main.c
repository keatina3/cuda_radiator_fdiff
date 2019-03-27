#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/time.h>
#include "utils.h"
#include "findiff.h"
#include "radiator.h"

extern void fdiff_gpu(double* u_vals, double* temps, int n, int m, int p, int block_size_Y, Tau* tau, int mallocPitch, int red);
extern void fdiff_gpu_glob(double* u_vals, double* temps, int n, int m, int p, int block_size, Tau* tau, int red);

int is_empty(FILE* file);
void write_times(char* fname, int n, int m, int p, int block_size, Tau tau, int GPU);

int main(int argc, char **argv){
	int m = 32, n = 32, p = 10, block_size = 16;
	int t = 0, w = 0, a = 0, option = 0, serial = 0, mallocPitch = 1, glob = 0;
	double **uold, **unew;
	double *u, *uold_vals, *unew_vals, *temps, *uGPU, *tempsGPU, *uglob_GPU, *tempsGPU2;
	struct timeval start, end;
	Tau tau_shared, tau_glob;

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
	
    ///////////////////////// SETUP //////////////////////////////////////
	
    uold = (double**)malloc(n*sizeof(double*));
	unew = (double**)malloc(n*sizeof(double*));
	uold_vals = (double*)calloc(n*m,sizeof(double));
	unew_vals = (double*)calloc(n*m,sizeof(double));
	uGPU = (double*)calloc(n*m,sizeof(double));
    uglob_GPU = (double*)calloc(n*m,sizeof(double));
    temps = (double*)calloc(n, sizeof(double));	
    tempsGPU = (double*)calloc(n, sizeof(double));	
    tempsGPU2 = (double*)calloc(n, sizeof(double));	
	
    init_mat(uold, uold_vals, n, m);	
	init_mat(unew, unew_vals, n, m);	
	
	apply_bounds(unew, n, 2, &fleft, &fright, &fzero, &fzero);
	apply_bounds(uold, n, 2, &fleft, &fright, &fzero, &fzero);
    
    memcpy(uGPU, unew_vals, n*m*sizeof(double));
    memcpy(uglob_GPU, unew_vals, n*m*sizeof(double));
    //////////////////////////////////////////////////////////////////////
    
    ///////////////////////// SERIAL /////////////////////////////////////
    
    if(serial){
        gettimeofday(&start, NULL);
	    iterate(uold, unew, p, n, m);
        gettimeofday(&end, NULL);
        tau_shared.calc_CPU = (end.tv_sec-start.tv_sec)*(1E6) +(end.tv_usec - start.tv_usec);
        tau_glob.calc_CPU = tau_shared.calc_CPU;
        
        if(p%2==0)
            u = uold_vals;
        else
            u = unew_vals;
        gettimeofday(&start, NULL);
	    if(a)
            red_rows(u, temps, n, m);
        gettimeofday(&end, NULL);
        tau_shared.calc_avgCPU = (end.tv_sec-start.tv_sec)*(1E6) + (end.tv_usec - start.tv_usec);
        tau_glob.calc_avgCPU = tau_shared.calc_avgCPU;
        tau_shared.tot_CPU = tau_shared.calc_avgCPU + tau_shared.calc_CPU;
        tau_glob.tot_CPU = tau_shared.tot_CPU;

	    if(m <= 10 && n <= 10)
            printf(MAG"\n\n=========================== PRINTING SERIAL GRID "
                            "========================================\n\n" RESET);
        else
            printf(MAG"\n\n===================== PRINTING (%dx%d) OF SERIAL GRID "
                            "===================================\n\n"RESET,10<n?10:n,10<m?10:m);
        print_grid(u,10<n?10:n,10<m?10:m,m);
        printf(MAG"=============================================================="
                            "===========================\n\n"MAG);
        if(w)
            write_times("CPU_times.csv", n, m, p, block_size, tau_shared, 0);
    }
    //////////////////////////////////////////////////////////////////////

    ///////////////////////// GLOBAL /////////////////////////////////////
    
    if(glob){
        gettimeofday(&start, NULL);
        fdiff_gpu_glob(uglob_GPU, tempsGPU2, n, m, p, block_size, &tau_glob, a);
        gettimeofday(&end, NULL);
        tau_glob.tot_GPU = (end.tv_sec-start.tv_sec)*(1E6) + (end.tv_usec - start.tv_usec);

	    if(m <= 10 && n <= 10)
            printf(MAG"\n\n=========================== PRINTING GPU (GLOBAL METHOD) GRID "
                            "===========================\n\n"RESET);
        else
            printf(MAG"\n\n===================== PRINTING (%dx%d) OF GPU (GLOBAL METHOD) GRID "
                            "======================\n\n"RESET,10<n?10:n,10<m?10:m);
        print_grid(uglob_GPU,10<n?10:n,10<m?10:m,m);
        printf(MAG"=============================================================="
                            "===========================\n\n"RESET);

        if(w)
            write_times("GPU_glob_times.csv", n, m, p, block_size, tau_glob, 1);
    }
    //////////////////////////////////////////////////////////////////////
    
    ///////////////////////// SHARED /////////////////////////////////////
    
    gettimeofday(&start, NULL);
    fdiff_gpu(uGPU, tempsGPU, n, m, p, block_size, &tau_shared, mallocPitch, a);
    gettimeofday(&end, NULL);
    tau_shared.tot_GPU = (end.tv_sec-start.tv_sec)*(1E6) + (end.tv_usec - start.tv_usec);

	if(m <= 10 && n <= 10)
        printf(MAG"\n\n=========================== PRINTING GPU (SHAMAG METHOD) GRID "
                            "===========================\n\n"RESET);
    else
        printf(MAG"\n\n===================== PRINTING (%dx%d) OF GPU (SHARED METHOD) GRID " 
                            "======================\n\n"RESET,10<n?10:n,10<m?10:m);
    print_grid(uGPU,10<n?10:n,10<m?10:m,m);
    printf(MAG"=============================================================="
                            "===========================\n\n"RESET);
    
    if(w){
        if(mallocPitch)
            write_times("GPU_sharedpitch_times.csv", n, m, p, block_size, tau_shared, 1);
        else
            write_times("GPU_sharednopitch_times.csv", n, m, p, block_size, tau_shared, 1);
    }
    //////////////////////////////////////////////////////////////////////
    
    ///////////////////////// RESULTS ////////////////////////////////////
    
    if(t){
        printf(MAG"\n========================= RESULTS ==========================\n\n"RESET);
        
        printf(GREEN"Using shared memory...\n"RESET);
        printf(RED"\t\tGPU time\tCPU time\tSpeedup\n"RESET);
        printf(BLUE"Fin diff calc:"RESET"\t%0.6f\t%0.6f\t%0.6f\n", tau_shared.calc_GPU, 
                        tau_shared.calc_CPU, tau_shared.calc_CPU/tau_shared.calc_GPU);
        if(a){printf(BLUE"Get avgs calc:"RESET"\t%0.6f\t%0.6f\t%0.6f\n", tau_shared.calc_avgGPU, 
                        tau_shared.calc_avgCPU, tau_shared.calc_avgCPU/tau_shared.calc_avgGPU);}
        printf(BLUE"Total time:"RESET"\t%0.6f\t%0.6f\t%0.6f\n", tau_shared.tot_GPU, 
                        tau_shared.tot_CPU, tau_shared.tot_CPU/tau_shared.tot_GPU);
        printf(BLUE"Alloc on GPU:"RESET"\t%0.6f\n",tau_shared.alloc_GPU);
        printf(BLUE"Transf RAM/GPU:"RESET"\t%0.6f\t%0.6f\n",tau_shared.transf_RAM, 
                        tau_shared.transf_GPU);

        if(glob){
            printf(GREEN"\nUsing global memory...\n"RESET);
            printf(RED"\t\tGPU(shared)\tGPU(glob)\tSpeedup\n"RESET);
            printf(BLUE"Fin diff calc:"RESET"\t%0.6f\t%0.6f\t%0.6f\n", tau_shared.calc_GPU, 
                        tau_glob.calc_GPU, tau_glob.calc_GPU/tau_shared.calc_GPU);
            if(a){printf(BLUE"Get avgs calc:"RESET"\t%0.6f\t%0.6f\t%0.6f\n", 
                                tau_shared.calc_avgGPU, tau_glob.calc_avgGPU, 
                                tau_glob.calc_avgGPU/tau_shared.calc_avgGPU);}
            printf(BLUE"Total time:"RESET"\t%0.6f\t%0.6f\t%0.6f\n", tau_shared.tot_GPU, 
                        tau_glob.tot_GPU, tau_glob.tot_GPU/tau_shared.tot_GPU);
            printf(BLUE"Alloc on GPU:"RESET"\t%0.6f\t%0.6f\t%0.6f\n",tau_shared.alloc_GPU, 
                        tau_glob.alloc_GPU, tau_glob.alloc_GPU/tau_shared.alloc_GPU);
            printf(BLUE"Transfer GPU:"RESET"\t%0.6f\t%0.6f\t%0.6f\n",tau_shared.transf_GPU, 
                        tau_glob.transf_GPU, tau_glob.transf_GPU/tau_shared.transf_GPU);
            printf(BLUE"Transfer RAM:"RESET"\t%0.6f\t%0.6f\t%0.6f\n",tau_shared.transf_RAM, 
                        tau_glob.transf_RAM, tau_glob.transf_RAM/tau_shared.transf_RAM);
        }
        if(serial){
            double SSEt = sse(temps, tempsGPU, n);
            double SSE2 = sse(uglob_GPU, u, m*n);
            double SSE = sse(uGPU, u, m*n);
            printf(GREEN"\nSSE vals...\n"RESET);
            printf(BLUE"Overall SSE:"RESET"\t%0.8f\n",SSE);
            printf(BLUE"SSE v. Glob:"RESET"\t%0.8f\n",SSE2);
            printf(BLUE"SSE of temp:"RESET"\t%0.8f\n",SSEt);
        }
        
        printf(MAG"============================================================\n\n"RESET);
    }
    //////////////////////////////////////////////////////////////////////
    
	free(uold); free(unew);
    free(temps); free(tempsGPU); free(tempsGPU2);
    free(uGPU); free(uglob_GPU);
	free(uold_vals); free(unew_vals);

	return 0;
}

int is_empty(FILE* file){
    size_t size;

    fseek(file, 0, SEEK_END);
    size=ftell(file);

    return size ? 0 : 1;
}

void write_times(char* fname, int n, int m, int p, int block_size, Tau tau, int GPU){
    FILE* fptr;

    fptr = fopen(fname, "a+");
    if(!fptr)
        printf("Couldn't open file %s\n", fname);
    
    if(is_empty(fptr))
        fprintf(fptr, "Block-size, NxM, p, FD Calc, Avg Calc, Total, Alloc, Transfer GPU, Transfer RAM\n");
    if(GPU){
        fprintf(fptr, "%d, %dx%d, %d, %f, %f, %f, %f, %f, %f\n", block_size, n, m, p, tau.calc_GPU, 
                tau.calc_avgGPU, tau.tot_GPU, tau.alloc_GPU, tau.transf_GPU, tau.transf_RAM);
    } else {
        fprintf(fptr, "%d, %dx%d, %d, %f, %f, %f, %f, %f, %f\n", block_size, n, m, p, tau.calc_CPU, 
                tau.calc_avgCPU, tau.tot_CPU, 0.0, tau.transf_GPU, tau.transf_RAM);
    }

    fclose(fptr);
}
