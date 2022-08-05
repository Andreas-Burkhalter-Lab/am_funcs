#include "mex.h"
typedef short int16;
typedef unsigned long uint32;

int readwave(int16 *wave, char *filename, int numch, long firstscan, long lastscan);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int mrows,ncols;
  int buflen,status;
  char *filename;
  int numch;
  long firstscan,lastscan;
  int dims[2];
  double *sP;
  int16 *wave;

  if (nrhs != 3)
    mexErrMsgTxt("Must have 3 input arguments");

  // Argument syntax: filename, number of channels, 2-vector of scan range
  // filename
  if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be the filename");
  if (mxGetM(prhs[0]) != 1)
    mexErrMsgTxt("Filename must be a row vector");
  buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  filename = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[0], filename, buflen);
  if(status != 0) 
    mexWarnMsgTxt("Not enough space. String is truncated.");

  // number of channels
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      !(mrows==1 && ncols==1) )
    mexErrMsgTxt("number of channels must be a noncomplex scalar double.");
  numch = mxGetScalar(prhs[1]);

  // scan numbers
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(mrows*ncols == 2) )
    mexErrMsgTxt("scan range must be a noncomplex 2-vector.");
  sP = mxGetPr(prhs[2]);
  firstscan = sP[0];
  lastscan = sP[1];

  // Allocate needed storage
  dims[0] = numch;
  dims[1] = lastscan-firstscan+1;
  plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS,mxREAL);
  wave = (int16*) mxGetPr(plhs[0]);

  // Do the call
  readwave(wave,filename,numch,firstscan,lastscan);
}

int readwave(int16 *wave, char *filename, int numch, long firstscan, long lastscan)
{
  FILE* fp;
  uint32 datastart;

  fp = fopen(filename,"rb");
  fread(&datastart,sizeof(uint32),1,fp);
  fseek(fp,datastart+firstscan*numch*sizeof(int16),SEEK_SET);
  fread(wave,sizeof(int16),(lastscan-firstscan+1)*numch,fp);
  fclose(fp);
}
