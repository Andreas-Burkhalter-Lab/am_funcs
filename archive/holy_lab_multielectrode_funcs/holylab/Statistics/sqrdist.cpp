#include "mex.h"
#include <pthread.h>
#include <unistd.h>  // for # of CPUs

template <class Tdata>
int sqrdist_thread_launcher(const Tdata *x,const Tdata *y,int N,int d,int q,double *distmtrx,int n_threads);

/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  double *x, *y, *distmtrx;
  int N,d,q;
  const mxArray *curarg;
  int ncpus,nthreads,ret;

  if (nrhs < 2)
    mexErrMsgTxt("sqrdist: requires two inputs");

  // Parse the inputs
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("sqrdist: x must be a real matrix");
  x = mxGetPr(curarg);
  d = mxGetM(curarg);
  N = mxGetN(curarg);

  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("sqrdist: y must be a real matrix");
  y = mxGetPr(curarg);
  if (d != mxGetM(curarg))
    mexErrMsgTxt("sqrdist: the dimensionality of x and y disagree");
  q = mxGetN(curarg);

  ncpus = sysconf(_SC_NPROCESSORS_ONLN);
  if (ncpus > 4)
    ncpus = 4;  // performance is not very linear, so don't go crazy
  nthreads = (int) (double(d)*double(q)*double(N)/100000);
  if (nthreads > ncpus)
    nthreads = ncpus;
  if (nthreads > q)
    nthreads = q;  // since we thread on q, no point if q == 1
  if (nthreads < 1)
    nthreads = 1;
  //nthreads = 1;
  //mexPrintf("nthreads = %d\n",nthreads);

  //mexPrintf("d %d, N %d, q %d, nlhs %d, nrhs %d\n",d,N,q,nlhs,nrhs);

  // Set up the output
  plhs[0] = mxCreateDoubleMatrix(N,q,mxREAL);
  distmtrx = mxGetPr(plhs[0]);
  // Do the actual work
  if (mxIsDouble(prhs[0]) && mxIsDouble(prhs[1]))
    sqrdist_thread_launcher<double>(x,y,N,d,q,distmtrx,nthreads);
  else if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1]))
    sqrdist_thread_launcher<float>((float *)x,(float *)y,N,d,q,distmtrx,nthreads);
  else
    mexErrMsgTxt("sqrdist: both x and y must either be single or double");
  return;
}


/*
 * Given x and y, compute the distance between them.  x is
 * d-by-N, y is d-by-q, and both are organized by column in memory.
 */
// with multithreading one might be able to get a performance boost by
// switching the order of the loops (fewer cache misses on x?), but
// this hasn't been tried yet.
template <class dataT>
void sqrdist_generic(const dataT *x,const dataT *y,int N,int d,int q,double *distmtrx)
{
  double distsq,distmtrxtemp,dtemp,sqrtmp;
  const dataT *xp;
  int i,j,k,identtemp;

  for (j = 0; j < q; j++) {   // Loop over points of y
    xp = x;                   // Restart at first x point
    for (i = 0; i < N; i++,distmtrx++) { // Loop over the points of x
      dtemp = 0;
      for (k = 0; k < d; k++,xp++) {  // Loop over coordinates
	sqrtmp = y[k]-*xp;
	dtemp += sqrtmp*sqrtmp;
      }
      *distmtrx = dtemp;
    }
    y += d;                     // Go on to the next point of y
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
  double *distmtrx;
};

// A wrapper for multithreading
template <class Tdata>
void *sqrdist_thread(void *p)
{
  thread_data_type<Tdata> *td;

  td = (thread_data_type<Tdata>*) p;
  sqrdist_generic<Tdata>(td->x,td->y,td->N,td->d,td->q,td->distmtrx);
}

// The launcher for multithreading
template <class Tdata>
int sqrdist_thread_launcher(const Tdata *x,const Tdata *y,int N,int d,int q,double *distmtrx,int n_threads)
{
  int i,j;
  int *splits,extras,n_pts_per_thread;
  thread_data_type<Tdata> *td;
  pthread_t *thread;

  if (n_threads == 1) {
    sqrdist_generic(x,y,N,d,q,distmtrx);
    return 1;
  }

  // Decide how the data will be split among threads
  n_pts_per_thread = q/n_threads;
  extras = q - n_threads*n_pts_per_thread;
  splits = new int[n_threads+1];
  splits[0] = 0;
  for (i = 0; i < n_threads; i++)
    splits[i+1] = splits[i] + n_pts_per_thread + (extras-- > 0);

  // Pack arguments into thread data type, splitting out the various
  // points in y to be handled by different CPUs.
  td = new thread_data_type<Tdata>[n_threads];
  for (i = 0; i < n_threads; i++) {
    td[i].x = x;
    td[i].y = y + d*splits[i];
    td[i].N = N;
    td[i].d = d;
    td[i].q = splits[i+1]-splits[i];
    td[i].distmtrx = distmtrx+N*splits[i];
  }

  delete[] splits;

  // Launch threads
  thread = new pthread_t[n_threads];
  for (i = 0; i < n_threads && td[i].q > 0; i++)
    if (pthread_create(thread+i,NULL,sqrdist_thread<Tdata>,(void*) (td+i)) != 0)
      return 0;  // error creating thread

  // Wait for threads to return
  // Only do this for the ones that got launched; i should be
  // initialized to the first thread that _didn't_ run (or to
  // n_threads, if q > n_threads)
  for (i--; i >= 0; i--)
    pthread_join(thread[i],NULL);

  delete[] td;
  delete[] thread;

  return 1;  // success!
}
