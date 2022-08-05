//read int8's and return matrix to matlab

#include "../commong/lfs_common.h"

#include "mex.h"
typedef char int8;

//--------------------------------------------------------------------------------------------
void readint8s(int8 *buf, int fd, int unitsize, long firstunit, long lastunit, long datastart)
{
   off_t pos=datastart+off_t(firstunit)*unitsize*sizeof(int8);
   off_t posNow=lseek(fd, pos, SEEK_SET);
   if(posNow==off_t(-1)){mexErrMsgTxt("error when lseek() inside readchars()."); }

   ssize_t status=read(fd, buf, sizeof(int8)*(lastunit-firstunit+1)*unitsize );
   if(status==-1){mexErrMsgTxt("error when read() inside readchars()."); }

   //todo: endian conversion

}

//--------------------------------------------------------------------------------------------
//arg: fd, unitsize, 2-vector of unit range, datastart

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int mrows,ncols;
  int buflen,status;
  char *filename;
  int unitsize;
  long firstunit,lastunit;
  int dims[2];
  double *sP;
  int8 *buf;

  if (nrhs != 4)
    mexErrMsgTxt("Must have 4 input arguments");

  // Argument syntax: fd, unitsize, 2-vector of unit range, datastart:

  // 0: file descriptor:
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !(mrows==1 && ncols==1) )
     mexErrMsgTxt("fd must be a noncomplex scalar double.");
  int fd = mxGetScalar(prhs[0]);
  
  // 1: unitsize
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !(mrows==1 && ncols==1) )
    mexErrMsgTxt("number of channels must be a noncomplex scalar double.");
  unitsize = mxGetScalar(prhs[1]);

  // 2:unit range
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(mrows*ncols == 2) )
    mexErrMsgTxt("scan range must be a noncomplex 2-vector.");
  sP = mxGetPr(prhs[2]);
  firstunit = sP[0];
  lastunit = sP[1];

  // 3: datastart
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || !(mrows==1 && ncols==1) )
     mexErrMsgTxt("datastart must be a noncomplex scalar double.");
  long datastart = mxGetScalar(prhs[3]);
  

  // Allocate needed storage
  dims[0] = unitsize;
  dims[1] = lastunit-firstunit+1;
  plhs[0] = mxCreateNumericArray(2,dims,mxINT8_CLASS,mxREAL);
  buf = (int8*) mxGetPr(plhs[0]);

  // Do the call
  readint8s(buf,fd,unitsize,firstunit,lastunit, datastart);
}

