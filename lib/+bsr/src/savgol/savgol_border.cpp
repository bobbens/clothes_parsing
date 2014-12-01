#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <float.h>
using namespace std;

#include "mex.h"



/*****************************************/
void  savgol(double *a_in, double *z, double ra, double rb, const double & theta, const int & h, const int & w, double *a_out) {
    
    ra = max(1.5, ra);
    rb = max(1.5, rb);
    const double ira2 = 1 / pow(ra, 2);
    const double irb2 = 1 / pow(rb, 2);
    const int wr = (int) floor(max(ra, rb));
    const double sint = sin(theta);
    const double cost = cos(theta);
    double d0, d1, d2, d3, d4, v0, v1, v2;
    int xi, yi, x, y, cpi;
    double di, ei, zi, di2, detA, invA, param;
    const double eps = exp(-300);
    
    for (int cp = 0; cp<(w*h); cp++){
        y = cp%h; x=cp/h;
        if ((x>=wr) && (x<(w-wr)) && (y>=wr) && (y<(h-wr))){
            a_out[cp] = a_in[cp];
        }
        else{
            
            d0=0; d1=0; d2=0; d3=0; d4=0;
            v0=0; v1=0; v2=0;
            for ( int u = -wr; u <= wr; u++ ){
                xi = x + u;
                if ((xi<0) || (xi>=w))
                    continue;
                for ( int v = -wr; v <= wr; v++ ){
                    yi = y + v;
                    if ((yi<0) || (yi>=h))
                        continue;
                    di = -u*sint + v*cost;
                    ei = u*cost + v*sint;
                    if ( (di*di*ira2 + ei*ei*irb2) > 1)
                        continue;
                    cpi = yi+xi*h;
                    zi = z[cpi];
                    di2 = di*di;
                    d0 = d0 + 1;
                    d1 = d1 + di;
                    d2 = d2 + di2;
                    d3 = d3 + di*di2;
                    d4 = d4 + di2*di2;
                    v0 = v0 + zi;
                    v1 = v1 + zi*di;
                    v2 = v2 + zi*di2;
                }
            }
            
            detA = -d2*d2*d2 + 2*d1*d2*d3 - d0*d3*d3 - d1*d1*d4 + d0*d2*d4;
            if (detA>eps)
                a_out[cp] = ((-d3*d3+d2*d4)*v0 + (d2*d3-d1*d4)*v1 + (-d2*d2+d1*d3)*v2)/ detA;
            
        }
    }
    
}
/*****************************************/
/*           MEX INTERFACE               */
/*****************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 5) mexErrMsgTxt("INPUT: (a_in, z, ra, rb, theta)");
    if (nlhs != 1) mexErrMsgTxt("OUTPUT: a_out");
    
    double* a_in = mxGetPr(prhs[0]);
    
    int h = mxGetM(prhs[0]);
    int w = mxGetN(prhs[0]);
    double* z = mxGetPr(prhs[1]);
    
    double ra = mxGetScalar(prhs[2]);
    double rb = mxGetScalar(prhs[3]);
    double theta = mxGetScalar(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
    double* a_out = mxGetPr(plhs[0]);
    
    savgol(a_in, z, ra, rb, theta, h, w, a_out);
}
