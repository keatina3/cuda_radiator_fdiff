#include "radiator.h"

double fleft(int xind, int yind, int nx, int ny){
    return 1.0*(double)(1+xind)/(float)(nx);
}

double fright(int xind, int yind, int nx, int ny){
    return 0.8*(double)(1+xind)/(float)(nx);
}

double fzero(int xind, int yind, int nx, int ny){
    return 0.0;
}
