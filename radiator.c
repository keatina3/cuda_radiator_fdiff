#include "radiator.h"

float fleft(int xind, int yind, int nx, int ny){
    return 1.0*(float)(1+xind)/(float)(nx);
}

float fright(int xind, int yind, int nx, int ny){
    return 0.8*(float)(1+xind)/(float)(nx);
}

float fzero(int xind, int yind, int nx, int ny){
    return 0.0;
}
