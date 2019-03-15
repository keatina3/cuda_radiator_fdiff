#ifndef _FINDIFF_H_
#define _FINDIFF_H_

void apply_bounds(float **u, int nx, int ny,
	float (*lbound)(int, int, int, int),
	float (*rbound)(int, int, int, int),
	float (*ubound)(int, int, int, int),
	float (*bbound)(int, int, int, int));
void iterate(float **uold, float **unew, int p, int nx, int ny);
float get_grid_avg(float **u, int nx, int ny);
void print_grid(float **u, int nx, int ny);

#endif
