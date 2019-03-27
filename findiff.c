#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "findiff.h"

void apply_bounds(float **u, int nx, int ny,
	float (*lbound)(int, int, int, int),
	float (*rbound)(int, int, int, int),
	float (*ubound)(int, int, int, int),
	float (*bbound)(int, int, int, int))
	{
	int i;

	// left & right boundary //
	for(i=0;i<nx;i++){
		u[i][0] = lbound(i, 0, nx, ny);
		u[i][ny-1] = rbound(i, (ny-1), nx, ny);
	}
	
	// upper & bottom boundary //
	for(i=1;i<(ny-1);i++){
		u[0][i] = bbound(0,i,nx,ny);
		u[nx-1][i] = ubound((nx-1),i,nx,ny);
	}
}

// serial propagation code //
void iterate(float **uold, float **unew, int p, int nx, int ny){
	int i,j,k;
	float **tmp;

	for(k=0;k<p;k++){
		for(i=0;i<nx;i++){
			for(j=2;j<ny;j++){
				unew[i][j] = (1.9*uold[i][(j-2)%ny] + 1.5*uold[i][(j-1)%ny] +
                    uold[i][j] + 0.5*uold[i][(j+1)%ny] + 0.1*uold[i][(j+2)%ny]);
				unew[i][j] /= (float)(5.0);
			}
		}
		tmp = uold;
		uold = unew;
		unew = tmp;
	}
    unew=uold; 
}

// reduce grid to one tmp //
float get_grid_avg(float *u, int nx, int ny){
	float *row_sum, temp;
	
	row_sum = (float*)calloc(nx,sizeof(float));
	red_rows(u, row_sum, nx, ny);
	temp = vec_reduce(row_sum, nx);
	free(row_sum);
	
	return temp;
}

void print_grid(float *u, int nx, int ny, int offset){
	int i,j;
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++)
			printf("%f ", u[j+i*offset]);
		printf("\n");
	}
}
