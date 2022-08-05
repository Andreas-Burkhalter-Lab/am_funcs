#include "mex.h"
#include <pthread.h>
#include <unistd.h>  // for # of CPUs

template <class Tdata>
int mindist_thread_launcher(const Tdata *x,const Tdata *y,int N,int d,int q,double *md,double *ident,int n_threads);

template <class Tdata>
void mindist_generic(const Tdata *x,const Tdata *y,int N,int d,int q,double *md,double *ident);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const double *x, *y;
  double *md, *ident;
  int N,d,q;
  const mxArray *curarg;
  int ncpus,nthreads,ret;

  if (nrhs < 2)
    mexErrMsgTxt("mindist: requires two inputs");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("mindist: x must be a real matrix");
  x = mxGetPr(curarg);
  d = mxGetM(curarg);
  N = mxGetN(curarg);

  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("mindist: y must be a real matrix");
  y = mxGetPr(curarg);
  if (d != mxGetM(curarg))
    mexErrMsgTxt("mindist: the dimensionality of x and y disagree");
  q = mxGetN(curarg);

  //mexPrintf("d %d, N %d, q %d, nlhs %d, nrhs %d\n",d,N,q,nlhs,nrhs);

  // Set up the outputs
  plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
  md = mxGetPr(plhs[0]);
  ident = NULL;
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
    ident = mxGetPr(plhs[1]);
  }

  ncpus = sysconf(_SC_NPROCESSORS_ONLN);
  /*
  if (ncpus > 4)
    ncpus = 4;  // performance is not very linear, so don't go crazy
  */
  nthreads = (int) ((double(d)*double(q)*double(N))/100000);
  if (nthreads > ncpus)
    nthreads = ncpus;
  if (nthreads > N)
    nthreads = N;  // since we thread on N, no point if N == 1
  if (nthreads < 1)
    nthreads = 1;
  //mexPrintf("nthreads = %d (ncpus = %d)\n",nthreads,ncpus);

  // Do the actual work
  if (mxIsDouble(prhs[0]) && mxIsDouble(prhs[1])) {
    ret = mindist_thread_launcher<double>(x,y,N,d,q,md,ident,nthreads);
    //mindist_generic<double>(x,y,N,d,q,md,ident);
  }
  else if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1]))
    ret = mindist_thread_launcher<float>((float *)x,(float *)y,N,d,q,md,ident,nthreads);
  else
    mexErrMsgTxt("mindist: x and y must both be either single or double");

  if (ret == 0)
    mexErrMsgTxt("mindist: error launching threads");

  return;
}


/*
 * Given x and y, compute the distance (md) between x and its closest
 * point in y, and the index (ident) of that closest point.  x is
 * d-by-N, y is d-by-q, and both are organized by column in memory.
 */
template <class Tdata>
void mindist_generic(const Tdata *x,const Tdata *y,int N,int d,int q,double *md,double *ident)
{
  double distsq,mdtemp,dtemp;
  Tdata sqrtmp;
  const Tdata *yp;
  int i,j,k,identtemp;

  for (i = 0; i < N; i++) {     // Loop over the points of x
    mdtemp = DBL_MAX;
    identtemp = -1;
    yp = y;                     // We're using pointer arithmetic to
				// avoid a multiplication
    for (j = 0; j < q; j++) {   // Loop over the points of y
      dtemp = 0;
      for (k = 0; k < d; k++,yp++) {  // Loop over coordinates
	sqrtmp = x[k]-*yp;
	dtemp += sqrtmp*sqrtmp;
      }
      if (dtemp < mdtemp) {     // Keep track of closest-yet
	mdtemp = dtemp;
	identtemp = j+1;        // j+1 because Matlab is unit-offset
      }
    }
    md[i] = mdtemp;
    if (ident != NULL)
      ident[i] = (double) identtemp;   // Because Matlab uses doubles...
    x += d;                     // Go on to the next point of x
  }
  return;
}

// A structure for passing data to pthreads
template <class Tdata>
struct thread_data_type {
  const Tdata *x;
  const Tdata *y;
  int d;
  int N;
  int q;
  double *md;
  double *ident;
};

// A wrapper for multithreading
template <class Tdata>
void *mindist_thread(void *p)
{
  thread_data_type<Tdata> *td;

  td = (thread_data_type<Tdata>*) p;
  mindist_generic<Tdata>(td->x,td->y,td->N,td->d,td->q,td->md,td->ident);
}

// The launcher for multithreading
template <class Tdata>
int mindist_thread_launcher(const Tdata *x,const Tdata *y,int N,int d,int q,double *md,double *ident,int n_threads)
{
  int i,j;
  int *splits,extras,n_pts_per_thread;
  thread_data_type<Tdata> *td;
  pthread_t *thread;

  if (n_threads == 1) {
    mindist_generic(x,y,N,d,q,md,ident);
    return 1;
  }

  // Decide how the data will be split among threads
  n_pts_per_thread = N/n_threads;
  extras = N - n_threads*n_pts_per_thread;
  splits = new int[n_threads+1];
  splits[0] = 0;
  for (i = 0; i < n_threads; i++)
    splits[i+1] = splits[i] + n_pts_per_thread + (extras-- > 0);

  // Pack arguments into thread data type, splitting out the various
  // points in x to be handled by different CPUs.
  td = new thread_data_type<Tdata>[n_threads];
  for (i = 0; i < n_threads; i++) {
    td[i].x = x + d*splits[i];
    td[i].y = y;
    td[i].N = splits[i+1]-splits[i];
    td[i].d = d;
    td[i].q = q;
    td[i].md = md+splits[i];
    if (ident != NULL)
      td[i].ident = ident+splits[i];
    else
      td[i].ident = NULL;
  }

  delete[] splits;

  // Launch threads
  thread = new pthread_t[n_threads];
  for (i = 0; i < n_threads && td[i].q > 0; i++)
    if (pthread_create(thread+i,NULL,mindist_thread<Tdata>,(void*) (td+i)) != 0)
      return 0;  // error creating thread

  // Wait for threads to return
  // Only do this for the ones that got launched; i should be
  // initialized to the first thread that _didn't_ run (or to
  // n_threads, if q > n_threads)
  for (i--; i >= 0; i--)
    pthread_join(thread[i],NULL);

  delete[] thread;
  delete[] td;

  return 1;  // success!
}

