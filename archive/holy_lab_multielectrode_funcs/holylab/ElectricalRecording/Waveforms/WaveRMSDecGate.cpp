// MATLAB gateway for WaveRMSDec
// From a raw waveform file, compute the mean(abs(waveform))
// in time bins
#include "mex.h"
#include "iotypes.h"
#include "FileHeaders.cp"
#include "Utils.h"
#include <stdio.h>

extern void WaveRMSDec(FilePtr &fpin,double *dout,const int32 numCh,const int32 nscans,const vector<int> &chIndex,int decimate);

const int printing = 0;

// Gateway routine: do the argument parsing, open the file, and
// allocate storage.
// Calling syntax: WaveRMSDec(filename,decfactor,trange,chanlist)
// where trange is optional (quoted in units of seconds) (default entire file)
// and chanlist is also optional (defaults to all channels)
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  char *filename = "";
  int decimate;
  long itrange[2];
  vector<int> chanIndex;

  FilePtr fpin;
  AIHeader aih;
  long nscans,ndec;
  double *dout,*dblP;

  int mrows,ncols,i;
  int cbufflen;
  int status;

  if (nrhs > 4 )
    mexErrMsgTxt( "Too many input arguments." );
  if (nrhs < 2 )
    mexErrMsgTxt( "Too few input arguments." );
  if (nlhs > 1 )
    mexErrMsgTxt( "Too many output arguments." );

  // filename: row vector of chars, convert to array of chars
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if (!mxIsChar(prhs[0]) || mrows != 1)
    mexErrMsgTxt("filename input must be a row vector of characters");
  cbufflen = mrows*ncols*sizeof(mxChar)+1;
  filename = (char*) mxCalloc(cbufflen,sizeof(char));
  status = mxGetString(prhs[0],filename,cbufflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space, filename string is truncated");

  // decimate: noncomplex double scalar, turn into int
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mrows*ncols != 1)
    mexErrMsgTxt("decimate input must be non-complex scalar double");
  if (mxGetScalar(prhs[1]) <= 0)
    mexErrMsgTxt("decimate must be positive");
  decimate = (int) mxGetScalar(prhs[1]);

  // Open the file now,
  // so can figure out default itrange & chanIndex
  fpin.open(filename,"rb");
  fpin >> aih;

  // trange: a 2-vector of doubles (convert to longs, in units of scans)
  itrange[0] = 0;
  itrange[1] = aih.nscans-1;
  if (nrhs >= 3) {
    mrows = mxGetM(prhs[2]);
    ncols = mxGetN(prhs[2]);
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mrows*ncols != 2)
      mexErrMsgTxt("trange input must be non-complex double 2-vector");
    dblP = mxGetPr(prhs[2]);
    // Now convert to integer scans
    itrange[0] = round(dblP[0]*aih.scanrate);
    itrange[1] = round(dblP[1]*aih.scanrate)-1;
    if (itrange[0] < 0)
      itrange[0] = 0;
    if (itrange[1] >= aih.nscans)
      itrange[1] = aih.nscans-1;
    if (itrange[1] < itrange[0]) {
    	// Return empty matrix
		plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
  }

  // Channel list: vector of doubles
  if (nrhs >= 4) {
    mrows = mxGetM(prhs[3]);
    ncols = mxGetN(prhs[3]);
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || (mrows > 1 && ncols > 1) )
      mexErrMsgTxt("chanlist input must be non-complex double vector");
    dblP = mxGetPr(prhs[3]);
    vector<int16> chanlist(mrows*ncols);
    for (i = 0; i < mrows*ncols; i++)
      chanlist[i] = dblP[i];
    // Now match up the selected channels to the recorded channels,
    // and arrange as an index
    MatchChannels(aih.channel,chanlist,chanIndex);
  }
  else {
    chanIndex.resize(aih.numCh);
    for (i = 0; i < aih.numCh; i++)
      chanIndex[i] = i;
  }


  // Allocate storage for output
  nscans = itrange[1]-itrange[0]+1;
  ndec = floor(float(nscans)/decimate);		// ideally this should be ceil, but needs fix in WaveRMSDec
  if (printing) {
    mexPrintf("file %s, decimate = %d, itrange [%d %d], outsize %d %d\n",filename,decimate,itrange[0],itrange[1],chanIndex.size(),ndec);
  }
  plhs[0] = mxCreateDoubleMatrix(chanIndex.size(),ndec,mxREAL);
  dout = mxGetPr(plhs[0]);

  // Advance file pointer to first desired scan
  fpin.skip(aih.numCh*itrange[0]);


  // Now do the real work
  WaveRMSDec(fpin,dout,aih.numCh,nscans,chanIndex,decimate);
}
