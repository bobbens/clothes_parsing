/* 
    Constrained segmentation by front propagation on Ultrametric Contour Map
    Source Code

    By Pablo Arbelaez
    arbelaez@eecs.berkeley.edu
    March 2008
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <queue>
#include <vector>
#ifndef __APPLE__
#include <values.h>
#else
#include <float.h>
#endif
using namespace std;

#include "mex.h"

#ifndef Active_h
#define Active_h

/*****************************************/
class Active
{
    public:
        double nrg;
        double lbl;
        int px;
        
        
        Active() { nrg = 0.0; lbl = -1; px = 0;}
        
        Active( const double& e, const double& l , const int& p) { nrg = e; lbl = l; px = p;}
        
        bool operator < ( const Active& x ) const { return ( (nrg > x.nrg) || ( (nrg == x.nrg) && ( lbl > x.lbl ) ) ); } 
        
};
#endif

/*****************************************/
void  UVT( double *ucm, double* markers, const int& tx, const int& ty,
double* labels, double* boundaries)
{
    //const double MAXDOUBLE = 1.7976931348623158e+308;
    // initialization
    priority_queue<Active, vector<Active>, less<Active> > band;
    bool* used = new bool[tx*ty];
    double* dist = new double[tx*ty]; 
		
	for (int p = 0; p < tx*ty; p++)
    {
        if( markers[p] > 0 )
        {
            labels[p] = markers[p];
            dist[p] = 0.0;
			boundaries[p]=0.0;
            band.push( Active(dist[p], labels[p], p) );
        }
        else
        {
            labels[p] = -1;
            dist[p] = DBL_MAX;
        }
        used[p] = false;
    }
    
    
    // propagation
    int vx[4] = { 1,  0, -1,  0};
    int vy[4] = { 0,  1,  0, -1};
    int cp, nxp, nyp, cnp;
    double u;
    
    while ( !band.empty() )
    {
        cp = band.top().px; band.pop();
        if (used[cp] == false)
        {
            for(int v = 0; v < 4; v++ )
            {
                nxp = (cp%tx) + vx[v]; nyp = (cp/tx) + vy[v]; cnp = nxp + nyp*tx;
                if ( (nyp >= 0) && (nyp < ty) && (nxp < tx) && (nxp >= 0) )
                {
 						u = max(dist[cp], ucm[cnp]);

					if (((u < dist[cnp])&&(labels[cnp]==-1)) || ( (u == dist[cnp]) && (labels[cnp] < labels[cp] ) ) )
                    {
                        labels[cnp] = labels[cp];
                        
						dist[cnp] = u;

                        band.push( Active(dist[cnp],labels[cnp], cnp) );
                    }
                }
            }
            used[cp] = true;
        }
    }
    
    delete[] used; delete[] dist;
	
	for (int cp = 0; cp < tx*ty; cp++)
		for(int v = 0; v < 4; v++ )
        {
			nxp = (cp%tx) + vx[v]; nyp = (cp/tx) + vy[v]; cnp = nxp + nyp*tx;
  			if ( (nyp >= 0) && (nyp < ty) && (nxp < tx) && (nxp >= 0) && (labels[cnp] < labels[cp]))
				boundaries[cp]=1;
		}


}
/*****************************************/
/*           MEX INTERFACE               */
/*****************************************/

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgTxt("INPUT: (ucm, seeds)");
    if (nlhs != 2) mexErrMsgTxt("OUTPUT: [ boundaries, labels]");
    
    double* ucm= mxGetPr(prhs[0]);
    
    int rows = mxGetM(prhs[0]);
    int cols = mxGetN(prhs[0]);
    
    double* markers = mxGetPr(prhs[1]);
    
    if( (mxGetM(prhs[1]) != rows) || (mxGetN(prhs[1]) != cols) )
        mexErrMsgTxt("ERROR: Input image and sources should have the same size.");
    
    plhs[0] = mxCreateDoubleMatrix(rows,cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows,cols, mxREAL);
    
    double* boundaries = mxGetPr(plhs[0]);
	double* labels = mxGetPr(plhs[1]);
    
    UVT(ucm, markers, rows, cols, labels, boundaries);
}
