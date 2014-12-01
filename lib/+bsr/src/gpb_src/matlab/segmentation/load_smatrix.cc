#include <string.h>
#include <math.h>
#include <mex.h>

/*
  code for reading in a sparse matrix

  n = fread(fid,1,'int');
  nnz = fread(fid,1,'int');
  x = zeros(nnz,1);
  y = zeros(nnz,1);
  z = zeros(nnz,1);
  ct = 1;
  for i = 1:n
    waitbar(i/n,h1);
    nz(i) = fread(fid,1,'int');
    vals = fread(fid,nz(i),'double');
    cols = fread(fid,nz(i),'int');
    x(ct:ct+nz(i)-1) = cols+1;
    y(ct:ct+nz(i)-1) = i*ones(nz(i),1);
    z(ct:ct+nz(i)-1) = vals;
    ct = ct + nz(i);
  end;
*/

void
mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // check number of arguments
    if (nlhs < 3) {
        mexErrMsgTxt("Too few output arguments.");
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (nrhs < 1) {
        mexErrMsgTxt("Too few input arguments.");
    }
    if (nrhs > 1) {
        mexErrMsgTxt("Too many input arguments.");
    }

    char *filename = static_cast<char*>(mxCalloc(mxGetN(prhs[0])+1, sizeof(char))); //mxCalloc is similar to malloc in C
    mxGetString(prhs[0],filename,mxGetN(prhs[0])+1);
    FILE* fp = fopen(filename,"r");
    if (fp != NULL)
    {
      int n = 0;
      fread(&n,sizeof(int),1,fp);
      int nnz = 0;
      fread(&nnz,sizeof(int),1,fp);

      plhs[0] = mxCreateDoubleMatrix(nnz,1,mxREAL);
      plhs[1] = mxCreateDoubleMatrix(nnz,1,mxREAL);
      plhs[2] = mxCreateDoubleMatrix(nnz,1,mxREAL);
      double* x = mxGetPr(plhs[0]); //col index
      int* xint = new int[nnz];
      double* y = mxGetPr(plhs[1]); //row index
      double* z = mxGetPr(plhs[2]); //values

      int ct = 0;
      int nz = 0;
      for (int row = 0; row < n; row++) 
      {
        fread(&nz,sizeof(int),1,fp); //number of entries in this row
        fread(&(z[ct]),sizeof(double),nz,fp); //value
        fread(&(xint[ct]),sizeof(int),nz,fp);    //col index
        for (int col = 0; col < nz; col++)
        {
          x[ct+col] = static_cast<double>(xint[ct+col] + 1); //add one for matlab indexing
          y[ct+col] = static_cast<double>(row + 1);
        } 
        ct = ct + nz;
      }
      fclose(fp);
      delete[] xint;
    }
    else
    {
      mexErrMsgTxt("Unable to open file for input.");
    }
}




