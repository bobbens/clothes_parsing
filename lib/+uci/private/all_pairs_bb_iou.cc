#include "mex.h"
#include <float.h>
#include <memory.h>

#if !defined(max)
#define	max(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(min)
#define	min(A, B)	((A) < (B) ? (A) : (B))
#endif


double intervals_overlap(double a1, double a2, double b1, double b2)
{
  double o1 = a2-b1;
  double o2 = b2-a1;
  double o3 = a2-a1;
  double o4 = b2-b1;
  double om12 = min(o1,o2);
  double om34 = min(o3,o4);
  double om = min(om12,om34);
  return max(om,0);
}

double bb_area(double minx, double miny, double maxx, double maxy)
{
  return (maxx-minx)*(maxy-miny);
}

// bounding-boxes expected in format: min_x, min_y, max_x, max_y
// computes intersection/union for all pairs
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int				nx, ny, dim;
	const double	*X, *Y, *px;
	double			*d2;
	double			d, dmin;
	int				i, j, k;
	double isect, ovx, ovy, ovar, sumar;

	if (nrhs != 2)
		mexErrMsgTxt("two input arguments expected.");
	if (nlhs != 1)
		mexErrMsgTxt("one output arguments expected.");

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
		mxGetNumberOfDimensions(prhs[0]) != 2)
		mexErrMsgTxt("input 1 (X) must be a real double matrix");

	dim = mxGetM(prhs[0]);
	nx = mxGetN(prhs[0]);

	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
		mxGetNumberOfDimensions(prhs[1]) != 2 ||
		mxGetM(prhs[1]) != dim)
		mexErrMsgTxt("input 2 (Y) must be a real double matrix compatible with input 1 (X)");

	ny = mxGetN(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(nx, ny, mxREAL);
	d2 = mxGetPr(plhs[0]);

	X = mxGetPr(prhs[0]);
	Y = mxGetPr(prhs[1]);


	for (i = 1; i <= ny; i++, Y += dim)
	{
		dmin = DBL_MAX;
		for (j = 1, px = X; j <= nx; j++, px += dim)
		{
			ovx = intervals_overlap(Y[0], Y[2], px[0], px[2]); // overlap in x
			ovy = intervals_overlap(Y[1], Y[3], px[1], px[3]); // overlap in y
			ovar = ovx*ovy;
			sumar = bb_area(Y[0], Y[1], Y[2], Y[3]) + bb_area(px[0], px[1], px[2], px[3]); // possible speedup: compute bb areas once and reuse
			*(d2++) = ovar/(sumar-ovar);
		}
	}
}
