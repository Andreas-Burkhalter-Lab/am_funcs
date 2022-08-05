#include "mex.h"
#include <pthread.h>
#include <unistd.h>  // for # of CPUs

int war_thread_launcher(const float *v,int d,const double *w,int sniplen,int n_components,const double *t,long n_spikes,int n_threads,double *b);


/*
 * This is the Matlab wrapper
 */
// Syntax:
//   b = wdecomp_amplitude_rhs(v,w,t)
//   Someday: b = wdecomp_amplitude_rhs(v,w,t,tfine)

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const float *v;
  int d;
  const double *w,*t;
  double *b;
  long n_spikes;
  int sniplen,n_components;
  const mxArray *curarg;
  const mwSize *w_size;
  int n_dims;
  int n_cpus,n_threads,ret;

  if (nrhs != 3)
    mexErrMsgTxt("wdecomp_amplitude_rhs: requires three inputs");

  // Parse the inputs
  // Voltage waveform
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("wdecomp_amplitude_rhs: v must be a real single-precision matrix");
  v = (const float *) mxGetData(curarg);
  d = mxGetM(curarg);
  v -= d;  // to handle the unit-offset of MATLAB

  // Components
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("wdecomp_amplitude_rhs: w must be a real double-precision three-dimensional array");
  w = mxGetPr(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
  w_size = mxGetDimensions(curarg);
  if (w_size[0] != d)
    mexErrMsgTxt("wdecomp_amplitude_rhs: w must have the same number of channels as the data");
  sniplen = w_size[1];
  n_components = 1;
  if (n_dims > 2)
    n_components = w_size[2];
  if (n_dims > 3)
    mexErrMsgTxt("wdecomp_amplitude_rhs: w must be a three-dimensional array");

  // Spike times
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg))
    mexErrMsgTxt("wdecomp_amplitude_rhs: t must be a real double-precision vector");
  t = mxGetPr(curarg);
  n_spikes = mxGetNumberOfElements(curarg);

  
  // Set up the outputs
  plhs[0] = mxCreateDoubleMatrix(n_components,n_spikes,mxREAL);
  b = mxGetPr(plhs[0]);
  
  n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
  n_threads = n_cpus;

  // Do the actual work
  ret = war_thread_launcher(v,d,w,sniplen,n_components,t,n_spikes,n_threads,b);

  if (ret == 0)
    mexErrMsgTxt("wdecomp_amplitude_rhs: error launching threads");

  return;
}


// Computes the "dot product" between the voltage waveform and the components
void war_work(const float *v,int d,const double *w,int sniplen,int n_components,const double *t,long n_spikes,double *b)
{
  double sum;
  long snipIndex,componentLen;
  int componentIndex;
  const float *vP, *vStart;
  const double *wP, *wEnd;

  componentLen = ((long) d)*((long) sniplen);
  for (snipIndex = 0; snipIndex < n_spikes; snipIndex++) {
    vStart = v + ((long) d)*((long) t[snipIndex]);
    wP = w;
    for (componentIndex = 0; componentIndex < n_components; componentIndex++,b++) {
      wEnd = wP + componentLen;
      sum = 0;
      for (vP = vStart; wP < wEnd; vP++, wP++)
	sum += *vP * *wP;
      *b = sum;
    }
  }
}

// A structure for passing data to pthreads
struct thread_data_type {
  const float *v;
  int d;
  const double *w;
  int sniplen;
  int n_components;
  const double *t;
  long n_spikes;
  double *b;
};

// A wrapper for multithreading
void *war_thread(void *p)
{
  thread_data_type *td;

  td = (thread_data_type*) p;
  war_work(td->v,td->d,td->w,td->sniplen,td->n_components,td->t,td->n_spikes,td->b);
}

// The launcher for multithreading
int war_thread_launcher(const float *v,int d,const double *w,int sniplen,int n_components,const double *t,long n_spikes,int n_threads,double *b)
{
  int i,n_launched;
  long *splits,extras,n_pts_per_thread;
  thread_data_type *td;
  pthread_t *thread;

  if (n_threads == 1) {
    war_work(v,d,w,sniplen,n_components,t,n_spikes,b);
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
  for (i = 0; i < n_threads; i++) {
    td[i].v = v;
    td[i].d = d;
    td[i].w = w;
    td[i].sniplen = sniplen;
    td[i].n_components = n_components;
    td[i].t = t + splits[i];
    td[i].n_spikes = splits[i+1]-splits[i];
    td[i].b = b + ((long) n_components)*splits[i];
  }

  delete[] splits;

  // Launch threads
  thread = new pthread_t[n_threads];
  for (i = 0; i < n_threads && td[i].n_spikes > 0; i++)
    if (pthread_create(thread+i,NULL,war_thread,(void*) (td+i)) != 0)
      return 0;  // error creating thread

  // Wait for threads to return
  // Only do this for the ones that got launched; i should be
  // initialized to the first thread that _didn't_ run (or to
  // n_threads)
  n_launched = i;
  for (i = n_launched-1; i >= 0; i--)
    pthread_join(thread[i],NULL);

  delete[] thread;
  delete[] td;

  return 1;  // success!
}
