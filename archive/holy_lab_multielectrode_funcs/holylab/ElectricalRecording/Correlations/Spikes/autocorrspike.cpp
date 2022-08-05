// Compute the autocorrelation function for a set of spikes, represented
// by a vector of their arrival times
#include "mex.h"
#include <vector>

using namespace std;

template<class dataType>
void AutoCorrEnum(dataType *t,long n,dataType tmax,vector<dataType> &tac,vector<long> *tacindx)
{
  //tac.erase(tac.begin(),tac.end());
  //tacindx.erase(tacindx.begin(),tacindx.end());
  long i,j;
  i = 0;
  dataType dt;
  while (i < n-1) {
    j = 1;
    while (i+j < n && (dt = t[i+j]-t[i]) < tmax) {
      tac.push_back(dt);
      tacindx[0].push_back(i+1);
      tacindx[1].push_back(i+j+1);
      j += 1;
    }
    i += 1;
  }
}

template<class dataType>
vector<long> AutoCorrBin(dataType *t,long n,dataType tmax,int nbins)
{
  dataType binnumfac = (double) nbins/tmax;
  vector<long> nInBins(nbins+1,0);		// One extra for dt = tmax
	
  long i,j;
  i = 0;
  dataType dt;
  while (i < n-1) {
    j = 1;
    while (i+j < n && (dt = t[i+j]-t[i]) <= tmax) {
      nInBins[int(binnumfac*dt)]++;
      j += 1;
    }
    i += 1;
  }
  //nInBins.erase(--nInBins.end());		// Get rid of the extra one
  nInBins.erase(nInBins.end()-1);		// Get rid of the extra one
  return nInBins;
}


//
// The gateway routine
//
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
  bool using_double,using_single;
  const mxArray *curarg;

  // Argument parsing
  if (nrhs < 2 || nrhs > 3)
    mexErrMsgTxt("autocorrspike requires 2 or 3 inputs");

  // t
  curarg = prhs[0];
  int bad1 =  (!mxIsNumeric(curarg) || mxIsComplex(curarg)
	       || (mxGetN(curarg) != 1 && mxGetM(curarg) != 1));
  if (bad1 && !mxIsEmpty(curarg))
    mexErrMsgTxt("The first input to autocorrspike must be a real numeric vector");
  void *t = mxGetData(curarg);
  long n = mxGetN(curarg)*mxGetM(curarg);
  using_double = mxIsDouble(curarg);
  using_single = mxIsSingle(curarg);

  // tmax
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("The second input to autocorrspike must be a real scalar");
  double tmax = mxGetScalar(curarg);

  // nbins (if present)
  int binning = 0;
  int nbins;
  if (nrhs == 3) {
    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
      mexErrMsgTxt("The third input (if present) to autocorrspike must be a real scalar");
    else {
      binning = 1;
      nbins = (int) mxGetScalar(prhs[2]);
    }
  }

  if (binning && nlhs != 1)
    mexErrMsgTxt("When binning, one output is required");
  else if(nlhs < 1 || nlhs > 2)
    mexErrMsgTxt("When not binning, one or two outputs are required");
  vector<long> nInBin;
  if (binning) {
    if (using_double)
      nInBin = AutoCorrBin<double>((double*) t,n,(double) tmax,nbins);
    else if (using_single)
      nInBin = AutoCorrBin<float>((float *) t,n,(float) tmax,nbins);
    else
      mexErrMsgTxt("Data type not yet implemented");
    plhs[0] = mxCreateDoubleMatrix(1,nbins,mxREAL);
    double *outp = mxGetPr(plhs[0]);
    for (long i = 0; i < nbins; i++)
      outp[i] = nInBin[i];
    return;
  }
  else {
    vector<long> tacindx[2];
    if (using_double) {
      vector<double> tac;
      AutoCorrEnum<double>((double *) t,n,(double) tmax,tac,tacindx);
      plhs[0] = mxCreateDoubleMatrix(1,tac.size(),mxREAL);
      double *outp = mxGetPr(plhs[0]);
      for (long i = 0; i < tac.size(); i++)
	outp[i] = tac[i];
    } else if (using_single) {
      vector<float> tac;
      AutoCorrEnum<float>((float *) t,n,(float) tmax,tac,tacindx);
      plhs[0] = mxCreateNumericMatrix(1,tac.size(),mxSINGLE_CLASS,mxREAL);
      float *outp = (float *) mxGetData(plhs[0]);
      for (long i = 0; i < tac.size(); i++)
	outp[i] = tac[i];
    }
    if (nlhs == 2) {
      int n_pairs = tacindx[0].size();
      plhs[1] = mxCreateDoubleMatrix(2,n_pairs,mxREAL);
      double *outp = mxGetPr(plhs[1]);
      for (long i = 0; i < n_pairs; i++) {
	outp[2*i] = tacindx[0][i];
	outp[2*i+1] = tacindx[1][i];
      }	
    }
  }
}
