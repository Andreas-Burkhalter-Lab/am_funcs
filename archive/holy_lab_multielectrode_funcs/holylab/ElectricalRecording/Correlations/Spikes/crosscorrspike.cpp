// Compute the crosscorrelation function for a set of spikes, represented
// by vectors of their arrival times
#include "mex.h"
#include <math.h>
#include <vector>

using namespace std;

template <class dataType>
void CrossCorrEnum(dataType *t1,dataType *t2,long n1,long n2,dataType tmax,vector<dataType> &tcc,vector<long> *tccindx)
{
  long i,j;
  j = 0;
  dataType dt;
  for (i = 0; i < n1; i++) {
    if (j >= n2) j = n2-1;
    // Back up if necessary
    while (j > 0 && t2[j]+tmax > t1[i]) j--;
    // Advance if necessary
    while (j < n2 && t2[j]+tmax < t1[i]) j++;
    //mexPrintf("i = %d: start j = %d,",i,j);
    // Grab the appropriate range
    while (j < n2 && fabs(dt = t2[j]-t1[i]) < tmax) {
      tcc.push_back(dt);
      tccindx[0].push_back(i+1);
      tccindx[1].push_back(j+1);
      j++;
    }
    //mexPrintf("finish j = %d\n",j);
  }
}

template <class dataType>
vector<long> CrossCorrBin(dataType *t1,dataType *t2,long n1,long n2,dataType tmax,int nbins)
{
  dataType binnumfac = double(nbins)/(2*tmax);
  vector<long> nInBins(nbins+1,0);		// One extra for dt = tmax (roundoff errors)
	
  long i,j;
  j = 0;
  dataType dt;
  if (n1*n2 == 0)                               // TEH 2004-07-23 crash fix
    return nInBins;                             // If one is empty...
  for (i = 0; i < n1; i++) {
    if (j >= n2) j = n2-1;
    // Back up if necessary
    while (j > 0 && t2[j]+tmax > t1[i]) j--;
    // Advance if necessary
    while (j < n2 && t2[j]+tmax < t1[i]) j++;
    //mexPrintf("i = %d: start j = %d,",i,j);
    // Grab the appropriate range
    while (j < n2 && fabs(dt = t2[j]-t1[i]) <= tmax) {
      nInBins[int(binnumfac*(dt+tmax))]++;
      j++;
    }
    //mexPrintf("finish j = %d\n",j);
  }
  //nInBins.erase(--nInBins.end());		// Get rid of the extra one
  //nInBins.erase(nInBins.end()-1);		// Get rid of the extra one
  return nInBins;
}

// The gateway routine
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  void *t1,*t2;
  mxClassID dataType;

  // Argument parsing
  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("crosscorrspike requires 3 or 4 inputs");


  // t1
  int bad1 =  (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1));
  if (bad1 && !mxIsEmpty(prhs[0]))
    mexErrMsgTxt("The first input to crosscorrspike must be a real numeric vector");
  t1 = mxGetData(prhs[0]);
  long n1 = mxGetNumberOfElements(prhs[0]);
  dataType = mxGetClassID(prhs[0]);

  // t2
  bad1 =  (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetN(prhs[1]) != 1 && mxGetM(prhs[1]) != 1));
  if (bad1 && !mxIsEmpty(prhs[1]))
    mexErrMsgTxt("The second input to crosscorrspike must be a real numeric vector");
  if (mxGetClassID(prhs[1]) != dataType)
    mexErrMsgTxt("The class types of t1 and t2 must agree");
  t2 = mxGetData(prhs[1]);
  long n2 = mxGetNumberOfElements(prhs[1]);

  // tmax
  if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
    mexErrMsgTxt("The third input to Crosscorrspike must be a real numeric scalar");
  double tmax = mxGetScalar(prhs[2]);

  // nbins
  int binning;
  int nbins;
  binning = 0;
  if (nrhs == 4) {
    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
      mexErrMsgTxt("The fourth input (if present) to crosscorrspike must be a real scalar");
    else {
      binning = 1;
      nbins = (int) mxGetScalar(prhs[3]);
    }
  }

  /*
  // Echo the inputs
  mexPrintf("First vector: ");
  for (int i = 0; i < n1; i++)
    mexPrintf("%g ",t1[i]);
  mexPrintf("\nSecond vector: ");
  for (int i = 0; i < n2; i++)
    mexPrintf("%g ",t2[i]);
  mexPrintf("\ntmax = %g, binning = %d",tmax,binning);
  if (binning)
    mexPrintf(", nbins = %d",nbins);
  mexPrintf("\n");
  */

  if (binning && nlhs != 1)
    mexErrMsgTxt("When binning, one output is required");
  else if(nlhs < 1 || nlhs > 2)
    mexErrMsgTxt("When not binning, one or two outputs are required");
  if (binning) {
    vector<long> nInBin;
    if (dataType == mxDOUBLE_CLASS)
      nInBin = CrossCorrBin<double>((double *)t1,(double *)t2,n1,n2,(double) tmax,nbins);
    else if (dataType == mxSINGLE_CLASS)
      nInBin = CrossCorrBin<float>((float *)t1,(float *)t2,n1,n2,(float) tmax,nbins);
    else
      mexErrMsgTxt("Data type not yet implemented");
    plhs[0] = mxCreateDoubleMatrix(1,nbins,mxREAL);
    double *outp = mxGetPr(plhs[0]);
    for (long i = 0; i < nbins; i++)
      outp[i] = nInBin[i];
    return;
  }
  else {
    vector<long> tccindx[2];
    if (dataType == mxDOUBLE_CLASS) {
      vector<double> tcc;
      CrossCorrEnum<double>((double *) t1,(double *) t2,n1,n2,(double) tmax,tcc,tccindx);
      plhs[0] = mxCreateDoubleMatrix(1,tcc.size(),mxREAL);
      double *outp = mxGetPr(plhs[0]);
      for (long i = 0; i < tcc.size(); i++)
	outp[i] = tcc[i];
    } else if (dataType == mxSINGLE_CLASS) {
      vector<float> tcc;
      CrossCorrEnum<float>((float *) t1,(float *) t2,n1,n2,(float) tmax,tcc,tccindx);
      plhs[0] = mxCreateNumericMatrix(1,tcc.size(),mxSINGLE_CLASS,mxREAL);
      float *outp = (float *) mxGetData(plhs[0]);
      for (long i = 0; i < tcc.size(); i++)
	outp[i] = tcc[i];
    }
    if (nlhs == 2) {
      long n_pairs = tccindx[0].size();
      plhs[1] = mxCreateDoubleMatrix(2,n_pairs,mxREAL);
      double *outp = mxGetPr(plhs[1]);
      for (long i = 0; i < n_pairs; i++) {
	outp[2*i] = tccindx[0][i];
	outp[2*i+1] = tccindx[1][i];
      }	
    }
  }
}
