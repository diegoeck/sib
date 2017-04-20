#include <stdint.h>
typedef uint16_t char16_t;
#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include "lapack.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* 3x3 matrix A
     * 76 25 11
     * 27 89 51
     * 18 60 32
     */
    double A[4] = {1, 2, 3, 4};
    double b[2] = {5, 6};

    long N = 2;
    long M = 1;

    long ipiv[2];
    long info;
    
    dgesv(&N, &M, A, &N, ipiv, b, &N, &info);
 
    mexPrintf("%1.10f\n",b[0]);
    mexPrintf("%1.10f\n",b[1]);

    mexPrintf("%d\n",ipiv[0]);
    mexPrintf("%d\n",ipiv[1]);

    
}
