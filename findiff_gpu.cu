#include "utils.h"
#include "findiff_gpu.h"
 
extern "C" {
void fd_iterate_gpu(float* u_vals, int block_size_X, int block_size_Y, int n, int m, int p, Tau* tau);

}
}
