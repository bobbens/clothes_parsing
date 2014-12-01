#define INF 1E20
#include <math.h>
#include <sys/types.h>
#include "mex.h"

/*
 * Generalized distance transforms based on Felzenswalb and Huttenlocher.
 */

static inline int square(int x) { return x*x; }

void dt1d(double *src, double *dst, int *ptr, int step, int n, double a, double b) {
  int   *v = new int[n];
  float *z = new float[n+1];
  int k = 0;
  v[0] = 0;
  z[0] = -INF;
  z[1] = +INF;
  for (int q = 1; q <= n-1; q++) {
    float s = ((src[q*step] - src[v[k]*step]) - b*(q - v[k]) + a*(square(q) - square(v[k]))) / (2*a*(q-v[k]));
    while (s <= z[k]) {
      // Update pointer
      k--;
      s  = ((src[q*step] - src[v[k]*step]) - b*(q - v[k]) + a*(square(q) - square(v[k]))) / (2*a*(q-v[k]));
    }
    k++;
    v[k]   = q;
    z[k]   = s;
    z[k+1] = +INF;
  }

   k = 0;
   for (int q = 0; q <= n-1; q++) {
     while (z[k+1] < q)
       k++;
     dst[q*step] = a*square(q-v[k]) + b*(q-v[k]) + src[v[k]*step];
     ptr[q*step] = v[k];
  }

  delete [] v;
  delete [] z;
}



// matlab entry point
// [M, Ix, Iy] = dt(vals, ax, bx, ay, by)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if (nrhs != 5)
    mexErrMsgTxt("Wrong number of inputs"); 
  if (nlhs != 3)
    mexErrMsgTxt("Wrong number of outputs");
  if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Invalid input");

  const int *dims = mxGetDimensions(prhs[0]);
  double *vals = (double *)mxGetPr(prhs[0]);
  double ax = mxGetScalar(prhs[1]);
  double bx = mxGetScalar(prhs[2]);
  double ay = mxGetScalar(prhs[3]);
  double by = mxGetScalar(prhs[4]);
  
  mxArray *mxM = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  mxArray *mxIx = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  mxArray *mxIy = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
  double *M = (double *)mxGetPr(mxM);
  int32_t *Ix = (int32_t *)mxGetPr(mxIx);
  int32_t *Iy = (int32_t *)mxGetPr(mxIy);

  double *tmpM = (double *)mxCalloc(dims[0]*dims[1], sizeof(double));
  int32_t *tmpIx = (int32_t *)mxCalloc(dims[0]*dims[1], sizeof(int32_t));
  int32_t *tmpIy = (int32_t *)mxCalloc(dims[0]*dims[1], sizeof(int32_t));

  for (int x = 0; x < dims[1]; x++)
    dt1d(vals+x*dims[0], tmpM+x*dims[0], tmpIy+x*dims[0], 1, dims[0], -ay, -by);
  
  for (int y = 0; y < dims[0]; y++)
    dt1d(tmpM+y, M+y, tmpIx+y, dims[0], dims[1], -ax, -bx);

  // get argmins and adjust for matlab indexing from 1
  for (int x = 0; x < dims[1]; x++) {
    for (int y = 0; y < dims[0]; y++) {
      int p = x*dims[0]+y;
      Ix[p] = tmpIx[p]+1;
      Iy[p] = tmpIy[tmpIx[p]*dims[0]+y]+1;
    }
  }

  mxFree(tmpM);
  mxFree(tmpIx);
  mxFree(tmpIy);
  plhs[0] = mxM;
  plhs[1] = mxIx;
  plhs[2] = mxIy;
}

/*
%%% DEBUGGING CODE %%%
N = 500;
A = rand(N,N);
w = rand(4,1);

tic; [res,ix,iy] = dt(A,w(1),w(2),w(3),w(4)); toc;
tic; [res2,ix2,iy2] = dt2(A,w(1),w(2),w(3),w(4)); toc;

norm(res(:) - res2(:))
mean(ix(:) ~= ix2(:))
mean(iy(:) ~= iy2(:))

*/
