// filterint16.c: MEX source for filtering electrophys. data
// Tim Holy, 06-21-01

// Can get speed improvements by commenting out the
// "#define float double" statement below, but you should be really
// careful that your filters are sufficiently stable. Use MATLAB's
// ZPLANE to check this, and/or compare against the result of using
// MATLAB's FILTER function.
#include "mex.h"
#include <math.h>
typedef short int16;
typedef unsigned long uint32;

const int printing = 0;

typedef int16 outtype;
#define float double

int filterwave(float *filtb, float *filta, int16 *wavein, int *chanIndex,
	       float *z, outtype *waveout, int totnchannels, int nchannels,
	       int nscans, int nfilt, int filtblen);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int totnchannels,nscans;
  int nchannels;
  int nfilt,filtblen,filtalen;
  int nrows,ncols;
  int status;
  int argnum;
  double *dfiltb,*dfilta;
  float *ffiltb,*ffilta;        // convert filters to single-precision (speed)
  double *dz;
  float *fz;
  double *dchanIndex;
  int *chanIndex;
  int16 *wavein;
  outtype *waveout;
  int dims[2];
  int i,j;
  int supplyInitialValues;
  int saveFinalValues;

  // Input argument syntax: filterb, filtera, wavein, chanIndex. Optionally: zi
  if (nrhs < 4 || nrhs > 5)
    mexErrMsgTxt("Must have 4 or 5 input arguments");
  supplyInitialValues = 0;
  if (nrhs == 5)
    supplyInitialValues = 1;
  // Output argument syntax: waveout. Optionally: zf.
  if (nlhs < 1 || nlhs > 2)
    mexErrMsgTxt("Must have 1 or 2 output arguments");
  saveFinalValues = 0;
  if (nlhs == 2)
    saveFinalValues = 1;

  // Parse the input arguments

  // filterb
  argnum = 0;
  filtblen = mxGetM(prhs[argnum]);
  nfilt = mxGetN(prhs[argnum]);
  if( !mxIsDouble(prhs[argnum]) || mxIsComplex(prhs[argnum]))
    mexErrMsgTxt("Filter b be a noncomplex double.");
  // a vector (either row or column) is assumed to be a single filter
  if (mxIsEmpty(prhs[argnum]))
    mexErrMsgTxt("Filter b cannot be empty");
  if (filtblen == 1)
    {
      filtblen = nfilt;
      nfilt = 1;
    }
  dfiltb = mxGetPr(prhs[argnum]);

  //filtera
  argnum = 1;
  if (!mxIsEmpty(prhs[argnum])) {
    filtalen = mxGetM(prhs[argnum]);
    ncols = mxGetN(prhs[argnum]);
    if( !mxIsDouble(prhs[argnum]) || mxIsComplex(prhs[argnum]))
      mexErrMsgTxt("Filter a be a noncomplex double.");
    if (filtalen == 1)
      {
	filtalen = ncols;
	ncols = 1;
      }
    if ( ((ncols != nfilt) || (filtalen != filtblen)) && filtalen > 1 )
      mexErrMsgTxt("Filters a and b are not of the same size, and a\nis not a scalar or empty.");
    dfilta = mxGetPr(prhs[argnum]);
  }
  else {
    dfilta = 0;
    filtalen = 0;
    if (printing)
      mexPrintf("filtera is empty\n");
  }

  //wave
  argnum = 2;
  totnchannels = mxGetM(prhs[argnum]);
  nscans = mxGetN(prhs[argnum]);
  if( !mxIsInt16(prhs[argnum]) || mxIsComplex(prhs[argnum]))
    mexErrMsgTxt("wave be a noncomplex int16.");
  wavein = (int16*) mxGetData(prhs[argnum]);

  //chanIndex
  argnum = 3;
  nrows = mxGetM(prhs[argnum]);
  ncols = mxGetN(prhs[argnum]);
  if( !mxIsDouble(prhs[argnum]) || mxIsComplex(prhs[argnum]))
    mexErrMsgTxt("chanIndex be a noncomplex double.");
  if ((nrows != 1) && (ncols != 1))
    mexErrMsgTxt("chanIndex must be a vector.");
  nchannels = nrows*ncols;
  if ((nchannels != nfilt) && (nfilt != 1))
    {
      mexPrintf("nchan %d, nfilt %d\n",nchannels,nfilt);
      mexErrMsgTxt("Number of filters does not agree with number of keeper channels.");
    }
  dchanIndex = mxGetPr(prhs[argnum]);
  
  //zi
  if (supplyInitialValues) {
    argnum = 4;
    nrows = mxGetM(prhs[argnum]);
    ncols = mxGetN(prhs[argnum]);
    if( !mxIsDouble(prhs[argnum]) || mxIsComplex(prhs[argnum]))
      mexErrMsgTxt("zi be a noncomplex double.");
    if ((ncols != nchannels) || (nrows != filtblen-1))
      mexErrMsgTxt("zi must be of size nchannels-by-filtlength-1.");
    dz = mxGetPr(prhs[argnum]);
  }

  // We've done the argument checking. Now convert
  // them into the desired data types.

  // For the filters, also divide everything by a(1)
  ffiltb = (float*) mxCalloc(nfilt*filtblen,sizeof(float));
  for (i = 0; i < nfilt*filtblen; i++)
    ffiltb[i] = dfiltb[i];


  if (filtalen > 1)
    ffilta = (float*) mxCalloc(nfilt*filtalen,sizeof(float));
  else
    ffilta = 0;         // if a is a scalar, we'll handle it by scaling
  // Scale filters by a(1)
  if (filtalen > 0)
    for (i = 0; i < nfilt; i++)
      for (j = 0; j < filtblen; j++)
	ffiltb[i*filtblen+j] /= dfilta[i*filtalen];
  if (ffilta)
    for (i = 0; i < nfilt; i++)
      for (j = 0; j < filtalen; j++)
	ffilta[i*filtalen+j] = dfilta[i*filtalen+j]/dfilta[i*filtalen];

  // Convert the channel index
  chanIndex = (int*) mxCalloc(nchannels,sizeof(int));
  for (i = 0; i < nchannels; i++)
    chanIndex[i] = dchanIndex[i] - 1;  // To account for unit-offset

  // Set up the history values
  // (initialized to 0 if not supplied)
  fz = (float*) mxCalloc(nchannels*(filtblen-1),sizeof(float));
  if (supplyInitialValues)
    for (i = 0; i < nchannels*(filtblen-1); i++)
      fz[i] = dz[i];

  if (printing)
    {
      mexPrintf("totnchannels %d, nscans %d, nchannels %d\n",totnchannels,
		nscans, nchannels);
      mexPrintf("nfilt %d, filtblen %d, filtalen %d\n", nfilt, filtblen,
		filtalen);
      mexPrintf("filterb: ");
      for (i = 0; i < nfilt; i++)
	{
	  for (j = 0; j < filtblen; j++)
	    mexPrintf("%g ",ffiltb[i*filtblen + j]);
	  mexPrintf("\n");
	}
      if (ffilta)
	{
	  mexPrintf("filtera: ");
	  for (i = 0; i < nfilt; i++)
	    {
	      for (j = 0; j < filtalen; j++)
		mexPrintf("%g ",ffilta[i*filtalen + j]);
	      mexPrintf("\n");
	    }
	}
    }

  // Allocate needed storage for waveout
  dims[0] = nchannels;
  dims[1] = nscans;
  plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS,mxREAL);
  waveout = (outtype*) mxGetData(plhs[0]);

  // Do the call
  filterwave(ffiltb,ffilta,wavein,chanIndex,fz,waveout,totnchannels,
	     nchannels,nscans,nfilt,filtblen);

  // If preserving the final z values, create the matrix
  // and copy over
  if (saveFinalValues) {
    plhs[1] = mxCreateDoubleMatrix(filtblen-1,nchannels,mxREAL);
    dz = mxGetPr(plhs[1]);
    for (i = 0; i < nchannels*(filtblen-1); i++)
      dz[i] = fz[i];
  }

  // Free temporary space
  mxFree(ffiltb);
  if (ffilta)
    mxFree(ffilta);
  mxFree(chanIndex);
  mxFree(fz);
}

int filterwave(float *filtb, float *filta, int16 *wavein, int *chanIndex,
	       float *z, outtype *waveout, int totnchannels, int nchannels,
	       int nscans, int nfilt, int filtblen)
{
  // This filtering is implemented by a direct form II transposed structure
  // See help for MATLAB's filter function
  int16 *currin,*endin;
  outtype *currout;
  float *currfiltb,*currfilta, *currz;
  float tempout;
  int i,j;

  if (printing)
    mexPrintf("filtb %0x, filta %0x, wavein %0x, chanIndex %0x, z %0x, waveout %0x\ntotnchannels %d, nchannels %d, nscans %d, nfilt %d, filtblen %d",
	      filtb,filta,wavein,chanIndex,z,waveout,
	      totnchannels,nchannels,nscans,nfilt,filtblen);

  endin = wavein + nscans*totnchannels;

  if (filta)
    {
      for (i = 0; i < nchannels; i++)
	{
	  if (nfilt > 1)
	    {
	      currfiltb = filtb + i*filtblen;
	      currfilta = filta + i*filtblen;
	    }
	  else
	    {
	      currfiltb = filtb;
	      currfilta = filta;
	    }
	  currz = z + i*(filtblen-1);
	  // Special case: if input filter(s) are actually scalars
	  if (filtblen < 2)
	    for (currin = wavein + chanIndex[i], currout = waveout + i;
		 currin < endin; currin += totnchannels, currout += nchannels)
	      *currout = (int16) floor(*currin * *currfiltb + 0.5);
	  else {
	    for (currin = wavein + chanIndex[i], currout = waveout + i;
		 currin < endin; currin += totnchannels, currout += nchannels)
	      {
		tempout = *currfiltb * *currin + *currz;
		/*
		if (printing > 1)
		  mexPrintf("*currfiltb %g, *currin %d, *currz %g, tempout %g\n",
			    *currfiltb, *currin, *currz, tempout);
			    */
		for (j = 0; j < filtblen-2; j++)
		  currz[j] = currfiltb[j+1] * *currin + currz[j+1]
		    - currfilta[j+1] * tempout;
		currz[filtblen-2] = currfiltb[filtblen-1] * *currin
		  - currfilta[filtblen-1] * tempout;
		*currout = (int16) floor(tempout+0.5);
	        //*currout = (int16) tempout;
	      }   // Loop over current channel
	  }    // else (filtblen >= 2)
	}   // Loop over channels
    }
  else   // filtera is all zeros, so save some operations
    {
      mexPrintf("Using filta = 0\n");
      for (i = 0; i < nchannels; i++)
	{
	  if (nfilt > 1)
	      currfiltb = filtb + i*filtblen;
	  else
	      currfiltb = filtb;
	  currz = z + i*(filtblen-1);
	  // Special case: if input filter(s) are actually scalars
	  if (filtblen < 2)
	    for (currin = wavein + chanIndex[i], currout = waveout + i;
		 currin < endin; currin += totnchannels, currout += nchannels)
	      *currout = (int16) floor(*currin * *currfiltb + 0.5);
	  else {
	    for (currin = wavein + chanIndex[i],currout = waveout + i;
		 currin < endin; currin += totnchannels, currout += nchannels)
	      {
		tempout = *currfiltb * *currin + *currz;
		/*
		if (printing > 1)
		  mexPrintf("*currfiltb %g, *currin %d, *currz %g, tempout %g\n",
			    *currfiltb, *currin, *currz, tempout);
			    */
		for (j = 0; j < filtblen-2; j++)
		  currz[j] = currfiltb[j+1] * *currin + currz[j+1];
		currz[filtblen-2] = currfiltb[filtblen-1] * *currin;
		*currout = (int16) floor(tempout+0.5);
		//*currout = (int16) tempout;
	      }   // Loop over current channel
	  }    // else (filtblen >= 2)
	}   // Loop over channels
    }
}
