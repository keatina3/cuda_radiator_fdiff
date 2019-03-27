#ifndef _FINDIFF_H_
#define _FINDIFF_H_

void apply_bounds(double **u, int nx, int ny,
	double (*lbound)(int, int, int, int),
	double (*rbound)(int, int, int, int),
	double (*ubound)(int, int, int, int),
	double (*bbound)(int, int, int, int));
void iterate(double **uold, double **unew, int p, int nx, int ny);
double get_grid_avg(double *u, int nx, int ny);
void print_grid(double *u, int nx, int ny, int offset);

#endif
