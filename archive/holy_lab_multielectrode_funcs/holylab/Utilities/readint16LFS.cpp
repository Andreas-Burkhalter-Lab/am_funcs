//read int16 and return matrix to matlab

#include "../commong/lfs_common.h"

#include "mex.h"

#if defined(UINT16SAMPLE)
typedef unsigned short int16;
mxClassID sampleClass=mxUINT16_CLASS;
#else
typedef short int16;
mxClassID sampleClass=mxINT16_CLASS;
#endif

typedef unsigned long uint32;

//--------------------------------------------------------------------------------------------
void readwave(int16 *wave, int fd, int numch, long firstscan, long lastscan, long long datastart)
{
   off_t pos=datastart+off_t(firstscan)*numch*sizeof(int16);
   off_t posNow=lseek(fd, pos, SEEK_SET);
   if(posNow==off_t(-1)){
     mexPrintf("fd %d, numch %d, firstscan %d, lastscan %d, datastart %lld, pos %lld\n",fd,numch,firstscan,lastscan,datastart,(long long)(pos));
     mexErrMsgTxt("readint16LFS: error when lseek() inside readwave()."); }

   ssize_t status=read(fd, wave, sizeof(int16)*(lastscan-firstscan+1)*numch );
   if(status==-1){mexErrMsgTxt("error when read() inside readwave()."); }

   //todo: endian conversion

}

//--------------------------------------------------------------------------------------------
//arg: fd, number of channels, 2-vector of scan range, datastart

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

  if (nrhs != 4)
    mexErrMsgTxt("Must have 4 input arguments");

  // Argument syntax: fd, number of channels, 2-vector of scan range, datastart:

  // 0: file descriptor:
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !(mrows==1 && ncols==1) )
     mexErrMsgTxt("fd must be a noncomplex scalar double.");
  int fd = mxGetScalar(prhs[0]);
  
  // 1: number of channels
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !(mrows==1 && ncols==1) )
    mexErrMsgTxt("number of channels must be a noncomplex scalar double.");
  numch = mxGetScalar(prhs[1]);

  // 2: scan numbers
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(mrows*ncols == 2) )
    mexErrMsgTxt("scan range must be a noncomplex 2-vector.");
  sP = mxGetPr(prhs[2]);
  firstscan = sP[0];
  lastscan = sP[1];

  // 3: datastart
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || !(mrows==1 && ncols==1) )
     mexErrMsgTxt("datastart must be a noncomplex scalar double.");
  long long datastart = mxGetScalar(prhs[3]);
  

  // Allocate needed storage
  dims[0] = numch;
  dims[1] = lastscan-firstscan+1;
  plhs[0] = mxCreateNumericArray(2,dims, sampleClass, mxREAL);
  wave = (int16*) mxGetPr(plhs[0]);

  // Do the call
  readwave(wave,fd,numch,firstscan,lastscan, datastart);
}

