#include "mex.h"
#include <pthread.h>
#include <unistd.h>  // for # of CPUs

int cfw_thread_launcher(const float *v,int d, const double *t,long n_spikes,const int *sniprange,const double *W,const double *A,int n_components,int n_threads,double *C);


/*
 * This is the Matlab wrapper
 */
// Syntax:
//   [C,badFlag] = covar_from_waveform(v,t,sniprange)
//   [C,badFlag] = covar_from_waveform(v,t,sniprange,W,A)

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *v;
  int d;
  const double *t;
  long n_snips;
  const double *sniprangeP;
  int sniprange[2];
  int sniplen;
  long ndimC;
  const double *W, *A;
  int n_components;
  double *C;
  const mxArray *curarg;
  int n_cpus,n_threads,ret;

  if (nrhs != 3 && nrhs != 5)
    mexErrMsgTxt("covar_from_waveform: requires three or five inputs");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("covar_from_waveform: v must be a real single-precision matrix");
  v = (const float *) mxGetData(curarg);
  d = mxGetM(curarg);
  v -= d;  // to handle the unit-offset of MATLAB

  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("covar_from_waveform: t must be a real double-precision vector");
  t = mxGetPr(curarg);
  n_snips = mxGetNumberOfElements(curarg);

  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg) || mxGetNumberOfElements(curarg) != 2)
    mexErrMsgTxt("covar_from_waveform: sniprange must be a real double-precision 2-element vector");
  sniprangeP = mxGetPr(curarg);
  sniprange[0] = (int) sniprangeP[0];
  sniprange[1] = (int) sniprangeP[1];
  sniplen = sniprange[1] - sniprange[0] + 1;

  ndimC = (long) d * (long) sniplen;

  //mexPrintf("n_snips %d, d %d, sniprange [%d %d], sniplen %d, ndimC %d\n",n_snips,d,sniprange[0],sniprange[1],sniplen,ndimC);

  W = NULL;
  A = NULL;
  n_components = 0;
  if (nrhs > 3) {
    curarg = prhs[3];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
      mexErrMsgTxt("covar_from_waveform: W must be a real double-precision matrix");
    W = mxGetPr(curarg);
    n_components = mxGetN(curarg);
    if (mxGetM(curarg) != ndimC)
      mexErrMsgTxt("covar_from_waveform: the size of W does not match sniprange and the dimensionality of the waveform");
    
    curarg = prhs[4];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
      mexErrMsgTxt("covar_from_waveform: A must be a real double-precision matrix");
    A = mxGetPr(curarg);
    if (mxGetM(curarg) != n_components)
      mexErrMsgTxt("covar_from_waveform: the number of components in A does not match W");
    if (mxGetN(curarg) != n_snips)
      mexErrMsgTxt("covar_from_waveform: the number of events in A does not match t");
  }
  
  // Set up the outputs
  //plhs[0] = mxCreateDoubleMatrix(1,ndimC*(ndimC+1)/2,mxREAL);
  plhs[0] = mxCreateDoubleMatrix(ndimC,ndimC,mxREAL);
  C = mxGetPr(plhs[0]);
  
  n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
  n_threads = n_cpus;

  // Do the actual work
  ret = cfw_thread_launcher(v,d,t,n_snips,sniprange,W,A,n_components,n_threads,C);

  if (ret == 0)
    mexErrMsgTxt("covar_from_waveform: error launching threads");

  return;
}


void cfw_work(const float *v,int d,const double *t,long n_spikes,const int *sniprange,const double *W,const double *A,int n_components,double *C)
{
  int sniplen;
  long ndimC,szC,spikeCounter;
  double *csnip, *csnipP, *csnipP2, *csnipEnd;
  const double *WP, *WEnd, *AP;
  const float *vL, *vP;
  double *CP, *CEnd;

  sniplen = sniprange[1] - sniprange[0] + 1;
  ndimC = (long) d * (long) sniplen;
  szC = ndimC*(ndimC+1)/2;
  csnip = new double[ndimC];
  csnipEnd = csnip + ndimC;
  if (W != NULL)
    WEnd = W + ndimC*n_components;
  CEnd = C+szC;

  for (CP = C; CP < CEnd; CP++)
    *CP = 0;

  for (spikeCounter = 0; spikeCounter < n_spikes; spikeCounter++) {     // Loop over the events in t
    vL = v + d*((long) t[spikeCounter] + sniprange[0]);
    // Copy the spike data into the current snippet, converting to
    // double (so it only has to be done once)
    for (csnipP = csnip, vP = vL; csnipP < csnipEnd; csnipP++,vP++)
      *csnipP = (double) *vP;
    if (W != NULL) {
      // Add the component-encoded part of the snippet
      for (WP = W,AP = A+n_components*spikeCounter; WP < WEnd; AP++)
	for (csnipP = csnip; csnipP < csnipEnd; csnipP++,WP++)
	  *csnipP += *WP * (*AP);
    }
    //for (csnipP = csnip; csnipP < csnipEnd; csnipP++)
    //  mexPrintf("%g ",*csnipP);
    //mexPrintf("\n");
    // Add the contribution to the covariance matrix
    for (csnipP = csnip,CP = C; csnipP < csnipEnd; csnipP++)
      for (csnipP2 = csnipP; csnipP2 < csnipEnd; csnipP2++, CP++)
	*CP += (*csnipP)*(*csnipP2);
  }

  delete[] csnip;
}


// A structure for passing data to pthreads
struct thread_data_type {
  const float *v;
  int d;
  const double *t;
  long n_spikes;
  const int *sniprange;
  const double *W;
  const double *A;
  int n_components;
  double *C;
};

// A wrapper for multithreading
void *cfw_thread(void *p)
{
  thread_data_type *td;

  td = (thread_data_type*) p;
  cfw_work(td->v,td->d,td->t,td->n_spikes,td->sniprange,td->W,td->A,td->n_components,td->C);
}

void vec2mtrx(double *Cvec,double *Cmtrx,long ndimC,long n_spikes)
{
  long k,km,szCvec;
  int i,j;
  double *CvecEnd;
  double Ctmp;

  szCvec = ndimC*(ndimC+1)/2;
  CvecEnd = Cvec+szCvec;
  for (i = 0, j = 0; Cvec < CvecEnd; Cvec++, i++) {
    if (i >= ndimC) {
      j++;
      i = j;
    }
    Ctmp = *Cvec/n_spikes;
    km = j*ndimC+i;
    Cmtrx[km] = Ctmp;
    km = i*ndimC+j;
    Cmtrx[km] = Ctmp;
  }
}

// The launcher for multithreading
int cfw_thread_launcher(const float *v,int d, const double *t,long n_spikes,const int *sniprange,const double *W,const double *A,int n_components,int n_threads,double *C)
{
  int i,sniplen,n_launched;
  long j,szC,ndimC;
  long *splits,extras,n_pts_per_thread;
  double *Ctmp;
  thread_data_type *td;
  pthread_t *thread;

  sniplen = sniprange[1] - sniprange[0] + 1;
  ndimC = (long) d * (long) sniplen;
  szC = ndimC*(ndimC+1)/2;  // The 1-d version of C

  if (n_threads == 1) {
    Ctmp = new double[szC];
    cfw_work(v,d,t,n_spikes,sniprange,W,A,n_components,Ctmp);
    vec2mtrx(Ctmp,C,ndimC,n_spikes);
    delete[] Ctmp;
    return 1;
  }

  // Decide how the data will be split among threads
  n_pts_per_thread = n_spikes/n_threads;
  extras = n_spikes - n_threads*n_pts_per_thread;
  splits = new long[n_threads+1];
  splits[0] = 0;
  for (i = 0; i < n_threads; i++)
    splits[i+1] = splits[i] + n_pts_per_thread + (extras-- > 0);

  // Pack arguments into thread data type, splitting out the various
  // points in t to be handled by different CPUs.
  td = new thread_data_type[n_threads];
  Ctmp = new double[szC*n_threads];
  for (i = 0; i < n_threads; i++) {
    td[i].v = v;
    td[i].d = d;
    td[i].t = t + splits[i];
    td[i].n_spikes = splits[i+1]-splits[i];
    td[i].sniprange = sniprange;
    td[i].W = W;
    if (W != NULL)
      td[i].A = A + n_components*splits[i];
    else
      td[i].A = NULL;
    td[i].n_components = n_components;
    td[i].C = Ctmp + i*szC;
  }

  delete[] splits;

  // Launch threads
  thread = new pthread_t[n_threads];
  for (i = 0; i < n_threads && td[i].n_spikes > 0; i++)
    if (pthread_create(thread+i,NULL,cfw_thread,(void*) (td+i)) != 0)
      return 0;  // error creating thread

  // Wait for threads to return
  // Only do this for the ones that got launched; i should be
  // initialized to the first thread that _didn't_ run (or to
  // n_threads)
  n_launched = i;
  for (i = n_launched-1; i >= 0; i--)
    pthread_join(thread[i],NULL);

  // Consolidate results from the different threads
  for (j = 0; j < szC; j++)
    for (i = 1; i < n_launched; i++)
      td[0].C[j] += td[i].C[j];
  // Reform into symmetric matrix
  vec2mtrx(td[0].C,C,ndimC,n_spikes);
    

  delete[] thread;
  delete[] td;
  delete[] Ctmp;

  return 1;  // success!
}
