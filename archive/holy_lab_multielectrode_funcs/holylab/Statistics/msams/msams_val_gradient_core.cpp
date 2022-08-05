#include <pthread.h>
#include <math.h>

// Copyright 2008 by Timothy E. Holy

// Types
// Data type for holding information to be passed to & from threads.
template <class Tdata>
struct thread_data_type {
  const Tdata *x;
  const Tdata *x0;
  const Tdata *kvec;
  int d;
  int N;
  Tdata *w;
  Tdata sumw;
  Tdata *step;
  Tdata *var;
  Tdata *gradsterm;
  Tdata *gradvterm;
};

// Function prototypes
template <class Tdata>
void msams_vg1(thread_data_type<Tdata>* p);
template <class Tdata>
void msams_vg2(thread_data_type<Tdata>* p);

// Wrappers for functions called by pthread.
// These just cast the void* to an appropriate pointer type
template <class Tdata>
void *msams_thread1(void *p)
{
  msams_vg1((thread_data_type<Tdata>*) p);
}

template <class Tdata>
void *msams_thread2(void *p)
{
  msams_vg2((thread_data_type<Tdata>*) p);
}

// This does housekeeping and launches the threads, and consolidates results.
template <class Tdata>
int msams_val_gradient_core(const Tdata *x,const Tdata *x0,const Tdata *kvec,int d,int N,Tdata *step,Tdata *var,Tdata *grad,Tdata *w,int n_threads)
{
  int i,j;
  int *splits,extras,n_pts_per_thread;
  int free_w;
  thread_data_type<Tdata> *td;
  pthread_t *thread;
  Tdata numerator,denominator,E,k2,gradsterm,gradvterm;
  Tdata sumw;

  // Decide how the data will be split among threads
  n_pts_per_thread = N/n_threads;
  extras = N - n_threads*n_pts_per_thread;
  splits = new int[n_threads+1];
  splits[0] = 0;
  for (i = 0; i < n_threads; i++)
    splits[i+1] = splits[i] + n_pts_per_thread + (extras-- > 0);

  // Allocate storage for weights, if need be
  free_w = 0;
  if (w == NULL) {
    w = new Tdata[N];
    free_w = 1;
  }

  // Pack arguments into thread data type
  td = new thread_data_type<Tdata>[n_threads];
  for (i = 0; i < n_threads; i++) {
    td[i].d = d;
    td[i].x0 = x0;
    td[i].kvec = kvec;
    /*
    // We make a copy of x0 and kvec for each thread to reduce the
    // chance that memory will have to move between CPUs for these
    // often-accessed values.
    // Same goes for var & step & grad: give each its own output, then
    // consolidate.
    td.x0 = new Tdata[d];
    td.kvec = new Tdata[d];
    memcpy(td.x0,x0,d*sizeof(Tdata));
    memcpy(td.kvec,kvec,d*sizeof(Tdata));
    */
    // The next inputs represent the portion of the data for which the
    // given thread is responsible
    td[i].x = x + d*splits[i];
    td[i].w = w + splits[i];
    td[i].N = splits[i+1]-splits[i];
    // Allocate storage & initialize partial outputs from this thread
    td[i].step = new Tdata[d];
    td[i].var = new Tdata[d];
    if (grad != NULL) {
      td[i].gradsterm = new Tdata[d];
      td[i].gradvterm = new Tdata[d];
    }
  }

  // Launch first round of threads
  thread = new pthread_t[n_threads];
  for (i = 0; i < n_threads; i++)
    if (pthread_create(thread+i,NULL,msams_thread1<Tdata>,(void*) (td+i)) != 0)
      return 0;  // error creating thread
  // Wait for threads to return
  for (i = 0; i < n_threads; i++)
    pthread_join(thread[i],NULL);
  // Consolidate the results from the first round
  for (j = 0; j < d; j++) {
    step[j] = 0;
    var[j] = 0;
    for (i = 0; i < n_threads; i++) {
      step[j] += td[i].step[j];
      var[j] += td[i].var[j];
    }
  }
  sumw = 0;
  for (i = 0; i < n_threads; i++)
    sumw += td[i].sumw;
  // Free the temporary storage that won't be needed again
  for (i = 0; i < n_threads; i++) {
    delete[] td[i].step;
    delete[] td[i].var;
  }

  // Calculation of the gradient: launch a new round of threads
  if (grad != NULL) {
    // Calculate the function value (needed in the gradient computation)
    numerator = 0;
    denominator = 0;
    for (j = 0; j < d; j++) {
      k2 = kvec[j]*kvec[j];
      numerator += k2*step[j]*step[j];
      denominator += k2*var[j];
    }
    E = numerator/denominator;
    for (i = 0; i < n_threads; i++) {
      td[i].step = step;
      td[i].var = var;
      if (pthread_create(thread+i,NULL,msams_thread2<Tdata>,(void*) (td+i)) != 0)
	return 0;  // error creating thread
    }
    // Wait for threads to return
    for (i = 0; i < n_threads; i++)
      pthread_join(thread[i],NULL);
    // Consolidate the results
    for (j = 0; j < d; j++) {
      gradsterm = 0;
      gradvterm = 0;
      for (i = 0; i < n_threads; i++) {
	gradsterm += td[i].gradsterm[j];
	gradvterm += td[i].gradvterm[j];
      }
      grad[j] = 2 * kvec[j] * ((step[j]*step[j] - 2*gradsterm) - E * (var[j] - 2*gradvterm)) / denominator;
    }
    // Delete the temporary storage
    for (i = 0; i < n_threads; i++) {
      delete[] td[i].gradsterm;
      delete[] td[i].gradvterm;
    }
  }

  // Properly scale each by sumw so user doesn't get confused in using
  // the step by itself.
  Tdata sumw2 = sumw*sumw;
  for (j = 0; j < d; j++) {
    step[j] /= sumw;
    var[j] /= sumw2;
  }

  // Clean up
  delete[] splits;
  if (free_w)
    delete[] w;
  delete[] td;
  delete[] thread;

  return 1;   // success!
}


// Calculate w, components of step & var
template <class Tdata>
void msams_vg1(thread_data_type<Tdata>* p)
{
  int d, N;
  int j, ell;
  const Tdata *xp;
  Tdata *dx,*k2;
  Tdata kdx2,wtmp,wdx,sumw;

  d = p->d;
  N = p->N;

  dx = new Tdata[d];
  k2 = new Tdata[d];

  for (j = 0; j < d; j++) {
    p->step[j] = 0;
    p->var[j] = 0;
    k2[j] = p->kvec[j]*p->kvec[j];
  }
  sumw = 0;

  xp = p->x;
  for (ell = 0; ell < N; ell++) {
    // Calculation of w
    kdx2 = 0;
    for (j = 0; j < d; j++, xp++) {
      dx[j] = *xp - (p->x0[j]);
      kdx2 += dx[j]*dx[j]*k2[j];
    }
    wtmp = exp(-kdx2);
    p->w[ell] = wtmp;
    sumw += wtmp;
    // Calculation of step & var component contributions
    for (j = 0; j < d; j++) {
      wdx = wtmp*dx[j];
      p->step[j] += wdx;
      p->var[j] += wdx*wdx;
    }
  }
  p->sumw = sumw;

  delete[] dx;
  delete[] k2;
}


// Calculate gradient terms (the portions that loop over all points)
template <class Tdata>
void msams_vg2(thread_data_type<Tdata>* p)
{
  int d, N;
  int i, j, ell;
  const Tdata *xp;
  Tdata *dx,*k2;
  Tdata ki2dxj2dxiw;

  d = p->d;
  N = p->N;

  dx = new Tdata[d];
  k2 = new Tdata[d];

  for (j = 0; j < d; j++) {
    p->gradsterm[j] = 0;
    p->gradvterm[j] = 0;
    k2[j] = p->kvec[j] * p->kvec[j];
  }

  xp = p->x;
  for (ell = 0; ell < N; ell++) {
    for (j = 0; j < d; j++, xp++)
      dx[j] = *xp - p->x0[j];
    for (j = 0; j < d; j++)
      for (i = 0; i < d; i++) {
	ki2dxj2dxiw = k2[i] * dx[j] * dx[j] * dx[i] * p->w[ell];
	p->gradsterm[j] += ki2dxj2dxiw * p->step[i];
	p->gradvterm[j] += ki2dxj2dxiw * p->w[ell] * dx[i];
      }
  }

  delete[] dx;
  delete[] k2;
}


