// Define MAIN while compiling to create a stand-alone command-line
// program. See Makefile.
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#include <callgrind.h>   // for profiling
#else
#include "mex.h"
#include "newdelete.h"
#endif

// Use the following for performance optimization
// (This must be before the include imiterators.cxx statement)
//#define MAXDIMS 4

#include <string.h>  // for memcpy and imiterators
#include <pthread.h>
#include <unistd.h>  // for # of CPUs
#include "imiterators.cxx"
#include "image_utilities.h"


/* Syntax:
 *   Ap = array_prolong(A,outSize)
 *   Ap = array_prolong(A,outSize,n_threads)
 */

// forward declarations
bool prolong_work_double(const double *A,int n_dims,const int *Asz,double *Aout,const int *Aoutsz,int n_threads);
bool prolong_work_float(const float *A,int n_dims,const int *Asz,float *Aout,const int *Aoutsz,int n_threads);


/*
 * This is the Matlab wrapper
 */

void mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  const void *A;
  void *Aout;
  int n_dims,dimIndex;
  const int *Asz;     // the size of A
  int *Aoutsz;        // the output size
  double *outSize;    // pointer to the desired output size (before sanity-checking)
  const mxArray *curarg;
  const mxLogical* dimFlag;
  bool isdouble;
  int n_threads,n_cpus;

  if (nrhs < 1 || nrhs > 3)
    mexErrMsgTxt("array_prolong: requires one to three inputs");
  if (nlhs != 1)
    mexErrMsgTxt("array_prolong: requires one output");

  // Parse the input
  // A
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("array_prolong: A must be a real array");
  if (mxIsSingle(curarg))
    isdouble = false;
  else if (mxIsDouble(curarg))
    isdouble = true;
  else
    mexErrMsgTxt("array_prolong: A must be single or double");
  if (mxIsEmpty(curarg))
    mexErrMsgTxt("array_prolong: does not work on empty arrays");

  A = mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
#ifdef MAXDIMS
  if (n_dims > MAXDIMS)
    mexErrMsgTxt("array_prolong: too many dimensions in input. Recompile?");
#endif
  Asz = mxGetDimensions(curarg);

  // Aoutsz
  curarg = prhs[1];
  if (!mxIsDouble(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("array_prolong: outSize must be a real double vector");
  if (mxGetNumberOfElements(curarg) != n_dims)
    mexErrMsgTxt("array_prolong: number of dimensions in A must match length of outSize");
  outSize = mxGetPr(curarg);
  // Convert to int and do sanity checking
  Aoutsz = (int *) mxMalloc(n_dims*sizeof(int));
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    Aoutsz[dimIndex] = (int) outSize[dimIndex];
    if (!(Aoutsz[dimIndex] == Asz[dimIndex] || 
	 Aoutsz[dimIndex] == 2*Asz[dimIndex]-2 ||
	  Aoutsz[dimIndex] == 2*Asz[dimIndex]-1))
      mexErrMsgTxt("array_prolong: requested output size is not feasible, given the input size");
  }

  // n_threads
  n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
  n_threads = n_cpus;
  //n_threads = 1;  // Currently not much benefit to multiple threads
  if (nrhs > 2) {
    curarg = prhs[2];
    if (mxGetNumberOfElements(curarg) != 1)
      mexErrMsgTxt("array_prolong: n_threads must be a scalar");
    n_threads = (int) mxGetScalar(curarg);
    if (n_threads < 1)
      n_threads = 1;
    if (n_threads > n_cpus)
      n_threads = n_cpus;
  }

  // Set up the output and do the work
  if (isdouble) {
    plhs[0] = mxCreateNumericArray(n_dims,Aoutsz,
				   mxDOUBLE_CLASS,mxREAL);
    Aout = mxGetData(plhs[0]);
    if (!prolong_work_double((double *) A,n_dims,Asz,(double *) Aout,Aoutsz,n_threads))
      mexErrMsgTxt("Error creating threads!");
  } else {
    plhs[0] = mxCreateNumericArray(n_dims,Aoutsz,
				   mxSINGLE_CLASS,mxREAL);
    Aout = mxGetData(plhs[0]);
    if (!prolong_work_float((float *) A,n_dims,Asz,(float *) Aout,Aoutsz,n_threads))
      mexErrMsgTxt("Error creating threads!");
  }

  mxFree(Aoutsz);

  return;
}



// A structure for passing data to pthreads
template <class dataType>
struct thread_data_type {
  // Frequently-written data
  int *dx;
  int *x;
  // Output data
  dataType *out;  
  // Thread-specific information (not written by thread)
  long outStart,outEnd;  // range to process
  // Never-written data (within a thread)
  int n_dims;          // number of dimensions
  const dataType *in;
  const int *szIn,*szOut;
  const long *pixelSkipIn,*pixelSkipOut;
  const int *expandsz;
  const long *expandszCum;
  dataType **coef;
  long **memoffset;
  const int *n;
};

// Wrapper for multithreading
template <class dataType>
void *prolong_thread_wrapper(void *p)
{
  thread_data_type<dataType> *td;

  td = (thread_data_type<dataType>*) p;
  prolong_thread(td);
}


void print_pixI_coords(const pixIterator &pI)
{
  for (int dimIndex = 0; dimIndex < pI.nDims(); dimIndex++)
    printf("%d ",pI.coord(dimIndex));
  printf(" (current %ld)\n",long(pI));
}

int dx2int(const int *dx,const long *pixelSkip,int n_dims)
{
  int ret = 0;
  for (int dimIndex = 1; dimIndex < n_dims; dimIndex++)
    ret += dx[dimIndex] * pixelSkip[dimIndex];
  return ret;
}

// This is the main calculation
template <class dataType>
void prolong_thread(thread_data_type<dataType> *td)
{
  const dataType *inP,*colStartP,*outEnd;
  dataType *outP;
  int dimIndex;
  long ptIndexIn,ptIndexOut,tmpI,tmpII;
  int baseCategory,thisCategory;
  const dataType *coefP,*coefEnd;
  long *offsetP;
  int n_nbrs;
  dataType val;

  const int n_dims = td->n_dims;
  int *dx = td->dx;
  int *x = td->x;

  // Initialize: convert the starting index to appropriate coordinates
  ptIndexOut = td->outStart;
  ptIndexIn = 0;
  for (dimIndex = n_dims-1; dimIndex >= 0; dimIndex--) {
    tmpI = x[dimIndex] = ptIndexOut/td->pixelSkipOut[dimIndex];
    ptIndexOut -= tmpI*td->pixelSkipOut[dimIndex];
    if (td->expandsz[dimIndex] > 1) {
      tmpII = tmpI >> 1;  // this dimension is expanding; divide size by 2
      dx[dimIndex] = tmpI - (tmpII << 1);  // remainder from divide by 2
    } else {
      tmpII = tmpI;  // this dimension is not expanding
      dx[dimIndex] = 0;
    }
    /*
    if (dimIndex == 0)
      colStartP = td->in + ptIndexIn;  // don't include first coord in colStart
    */
    ptIndexIn += tmpII * td->pixelSkipIn[dimIndex];
  }
  ptIndexOut = td->outStart;

  /*
  mexPrintf("\nptIndexOut = %d, ptIndexIn = %d\n",ptIndexOut,ptIndexIn);
  mexPrintf("x: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    mexPrintf("%d ",x[dimIndex]);
  mexPrintf("\ndx: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    mexPrintf("%d ",dx[dimIndex]);
  mexPrintf("\n\n");
  */

  // Loop over all output pixels
  inP = td->in + ptIndexIn;
  outP = td->out + ptIndexOut;
  outEnd = td->out + td->outEnd;
  const int *dxEnd = dx + n_dims;
  baseCategory = dx2int(dx,td->expandszCum,n_dims);
  while (outP < outEnd) {
    // Compute the output pixel value
    thisCategory = baseCategory + dx[0];
    coefP = td->coef[thisCategory];
    offsetP = td->memoffset[thisCategory];
    n_nbrs = td->n[thisCategory];
    //mexPrintf("ptIndexIn = %d, ptIndexOut = %d, thisCategory = %d, n_nbrs = %d\n",inP-td->in,outP-td->out,thisCategory,n_nbrs);
    coefEnd = coefP + n_nbrs;
    val = 0;
    //mexPrintf("Neighbor values: ");
    for (; coefP < coefEnd; coefP++, offsetP++) {
      val += *coefP * inP[*offsetP];
      //mexPrintf("%g ",inP[*offsetP]);
    }
    *outP = val;
    //mexPrintf("\n");
    // Increment: we tie the increment of the input to the increment
    // of the output, to keep everything in register
    outP++;
    if (outP < outEnd) {
      x[0]++;
      if (x[0] < td->szOut[0]) {
	// We're progressing down a single column, but we haven't
	// reached the end yet
	dx[0]++;
	if (dx[0] >= td->expandsz[0]) {
	  // We have finished the 1 or 2 output categories associated
	  // with this input base point and output column.  Therefore,
	  // progress to the next input point.
	  inP++;
	  dx[0] = 0;
	}
      } else {
	// We finished a "column" of the output, so go to the
	// beginning of the next output column. This is an
	// increment-with-carry.
	// Note in what follows we know there are at least two
	// dimensions, because we got here knowing that outP <
	// outEnd. So we can assume dx[1] exists, and we don't have to
	// check whether dimIndex >= n_dims.
	dimIndex = 0;
	bool test = true;
	while (test) {
	  dx[dimIndex] = x[dimIndex] = 0;
	  dimIndex++;
	  x[dimIndex]++;
	  dx[dimIndex]++;
	  test = (x[dimIndex] >= td->szOut[dimIndex]);
	}
	if (dx[dimIndex] >= td->expandsz[dimIndex])
	  dx[dimIndex] = 0;
	// Recalculate inP
	ptIndexIn = 0;
	for (dimIndex = n_dims-1; dimIndex >= 0; dimIndex--) {
	  tmpI = x[dimIndex];
	  if (td->expandsz[dimIndex] > 1)
	    tmpII = tmpI >> 1;  // expanding, divide size by 2
	  else
	    tmpII = tmpI;  // this dimension is not expanding
	  ptIndexIn += tmpII * td->pixelSkipIn[dimIndex];
	}
	inP = td->in+ptIndexIn;
	/*
	mexPrintf("Recalculate: ptIndexIn = %d\n",ptIndexIn);
	mexPrintf("x: ");
	for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
	  mexPrintf("%d ",x[dimIndex]);
	mexPrintf("\ndx: ");
	for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
	  mexPrintf("%d ",dx[dimIndex]);
	mexPrintf("\n\n");
	*/
	baseCategory = dx2int(dx,td->expandszCum,n_dims);
      }
    }
  }
}

template <class dataType>
bool prolong_work(const dataType *A,int n_dims_in,const int *Asz_in,dataType *Aout,const int *Aoutsz_in,int n_threads)
{
  // This function initializes everything: it sets up the weights, and
  // then configures the threads

  #ifdef MAIN
  CALLGRIND_START_INSTRUMENTATION;
  #endif

  int n_dims;
  int dimIndex,index;
  long nptsOut;
  
  // Eliminate singleton dimensions
  int *dimList = new int[n_dims_in];
  n_dims = skip_unity_dimensions_index(Aoutsz_in,n_dims_in,dimList,n_dims_in);
  int *Asz = new int[n_dims];
  int *Aoutsz = new int[n_dims];
  int max_nbrs = 1;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    Asz[dimIndex] = Asz_in[dimList[dimIndex]];
    Aoutsz[dimIndex] = Aoutsz_in[dimList[dimIndex]];
    max_nbrs *= 2;
  }
  delete[] dimList;

  // Calculate pixel skips
  long *Askip = new long[n_dims];
  long *Aoutskip = new long[n_dims];
  calc_pixel_skip(Asz,n_dims,Askip);
  nptsOut = calc_pixel_skip(Aoutsz,n_dims,Aoutskip);
  
  // Determine which dimensions are to be interpolated, and set up the
  // one-dimensional coefficient & n_nbrs table.  Table structure:
  // dim0:
  //   output0
  //   output1
  // dim1:
  //   output0
  //   output1
  // etc.
  // The coefficients are twice as big as the n_nbrs table because
  // there are up to 2 neighbors that need coefficients.
  int *expandsz = new int[n_dims];
  dataType *coef1d = new dataType[4*n_dims];
  int *n1d = new int[2*n_dims];
  dataType *coef1dP;
  int *n1dP;
  int n_out_per_cell = 1;  // total number of output "categories" (see below)
  for (dimIndex = 0, coef1dP = coef1d, n1dP = n1d; dimIndex < n_dims; dimIndex++, coef1dP += 4, n1dP += 2) {
    expandsz[dimIndex] = 1+(Asz[dimIndex] < Aoutsz[dimIndex]);
    n_out_per_cell *= expandsz[dimIndex];
    if (expandsz[dimIndex]>1) {
      if (Aoutsz[dimIndex] % 2 == 0) { // even case
	n1dP[0] = n1dP[1] = 2;          // two nbrs contribute in both cases
	coef1dP[0] = coef1dP[3] = 0.75; // output0,nbr0 and output1,nbr1
	coef1dP[1] = coef1dP[2] = 0.25; // output0,nbr1 and output1,nbr0
      } else { // odd case
	n1dP[0] = 1;             // only 1 nbr contributes for on-grid case
	n1dP[1] = 2;             // 2 nbrs contribute at half-grid case
	coef1dP[0] = 1;          // output0,nbr0 (the on-grid case)
	coef1dP[2] = coef1dP[3] = 0.5;  //output1,nbr0 and output1,nbr1
      }
    } else { // no expansion case (just copy)
      n1dP[0] = 1;           // only 1 output, only 1 nbr
      coef1dP[0] = 1;
    }
  }

  // Calculate the multidimensional coefficients & offsets
  //
  // The coefficients are products of the one-dimensional
  // coefficients. For any given "cell" of input points there will
  // multiple output points, each with different sets of coefficients:
  // these are called output "categories."
  //
  //  We use a somewhat complex scheme to avoid any zero
  // coefficients.  More important than saving computation, this also
  // ensures that no out-of-range components are addressed at edges of
  // the image.  One consequence is that each "output category" has its
  // own unique set of neighbors.  Therefore, the list of contributing
  // neighbors will vary with output category.
  //
  // There will be n_out_per_cell = prod(expandsz) different output categories.
  // We will look these up based on the even/odd status of the source
  // point and whether the given dimension is being expanded.
  dataType **coef = new dataType*[n_out_per_cell];
  dataType *coefList = new dataType[n_out_per_cell*max_nbrs];
  long **memoffset = new long*[n_out_per_cell];
  long *memoffsetList = new long[n_out_per_cell*max_nbrs];
  int *n = new int[n_out_per_cell];  // # of contributing neighbors
  pixIterator dx(expandsz,n_dims,false);  // iterator over output categories
  pixIterator nbrI;   // iterator over inputs contributing to given output
  int *thisn = new int[n_dims];
  dataType *coefP,coeftmp;
  long *memoffsetP,memoffsettmp;
  for (index = 0,coefP = coefList,memoffsetP = memoffsetList; !dx.at_end(); dx++,index++) {  // loop over the different "categories" of output pixels
    coef[index] = coefP;
    memoffset[index] = memoffsetP;
    //mexPrintf("coefP %d, memoffsetP %d\n",coefP-coefList,memoffsetP-memoffsetList);
    // Calculate the number of nbrs that contribute to the current output category
    for (dimIndex = 0, n1dP = n1d; dimIndex < n_dims; dimIndex++,n1dP += 2)
      thisn[dimIndex] = n1dP[dx.coord(dimIndex)];
    nbrI.initialize(thisn,n_dims,false);
    n[index] = nbrI.numel();
    //mexPrintf("nnbrs = %d\n",n[index]);
    // Loop over all of the nbrs and assign coefs & memoffsets
    for (; !nbrI.at_end(); nbrI++, coefP++, memoffsetP++) {
      coeftmp = 1;
      memoffsettmp = 0;
      for (dimIndex = 0, coef1dP = coef1d; dimIndex < n_dims; dimIndex++,coef1dP+=4) {
	coeftmp *= coef1dP[2*dx.coord(dimIndex)+nbrI.coord(dimIndex)];
	memoffsettmp += nbrI.coord(dimIndex) * Askip[dimIndex];
      }
      *coefP = coeftmp;
      *memoffsetP = memoffsettmp;
    }
  }
  delete[] thisn;
  delete[] coef1d;
  delete[] n1d;

  /*
  // Test it
  int tmpI;
  for (tmpI = 0, dx.restart(); !dx.at_end(); tmpI++,dx++) {
    mexPrintf("Output category: ");
    print_pixI_coords(dx);
    int ntmp = n[tmpI];
    mexPrintf("Neighbor memoffset (%d total): ",ntmp);
    memoffsetP = memoffset[tmpI];
    for (index = 0; index < ntmp; index++)
      mexPrintf("%d ",memoffsetP[index]);
    mexPrintf("\nCoef: ");
    coefP = coef[tmpI];
    for (index = 0; index < ntmp; index++)
      mexPrintf("%g ",coefP[index]);
    mexPrintf("\n\n");
  }
  */

  // Prepare for multithreading: pack data into the thread structure
  thread_data_type<dataType> tdbase;
  tdbase.n_dims = n_dims;
  tdbase.in = A;
  tdbase.out = Aout;
  tdbase.szIn = Asz;
  tdbase.szOut = Aoutsz;
  tdbase.pixelSkipIn = Askip;
  tdbase.pixelSkipOut = Aoutskip;
  tdbase.expandsz = expandsz;
  long *expandszCum = new long[n_dims];
  calc_pixel_skip(expandsz,n_dims,expandszCum);
  tdbase.expandszCum = expandszCum;
  tdbase.coef = coef;
  tdbase.memoffset = memoffset;
  tdbase.n = n;

  // Create the thread-specific info
  tdbase.dx = new int[n_dims];
  tdbase.x = new int[n_dims];
  tdbase.outStart = 0;
  double n_per_thread = (double) nptsOut / n_threads;
  tdbase.outEnd = long(n_per_thread);

  bool isgood = true;
  if (n_threads <= 1) {
    prolong_thread(&tdbase);
  }
  else {
    //mexPrintf("%d threads being used; n_per_thread = %g.\n",n_threads,n_per_thread);
    thread_data_type<dataType> *td;
    td = new thread_data_type<dataType>[n_threads];
    int threadIndex;
    //mexPrintf("Thread 0: dx = %p, x = %p, outStart %d, outEnd %d\n",tdbase.dx,tdbase.x,tdbase.outStart,tdbase.outEnd);
    for (threadIndex = 1; threadIndex < n_threads; threadIndex++) {
      memcpy(td+threadIndex,&tdbase,sizeof(tdbase));
      td[threadIndex].dx = new int[n_dims];
      td[threadIndex].x = new int[n_dims];
      td[threadIndex].outStart = long(threadIndex*n_per_thread);
      td[threadIndex].outEnd = long((threadIndex+1)*n_per_thread);
      //mexPrintf("Thread %d: dx = %p, x = %p, outStart = %d, outEnd = %d\n",threadIndex,td[threadIndex].dx,td[threadIndex].x,td[threadIndex].outStart,td[threadIndex].outEnd);
    }
    td[n_threads-1].outEnd = nptsOut;  // to avoid roundoff error issues
    //mexPrintf("nptsOut = %d (outEnd of last thread changed to this)\n",nptsOut);

    // Launch the threads. We don't create a thread 0, we run
    // that one in this thread
    pthread_t *thread;
    thread = new pthread_t[n_threads];
    for (threadIndex = 1; threadIndex < n_threads && isgood; threadIndex++)
      if (pthread_create(thread+threadIndex,NULL,prolong_thread_wrapper<dataType>,(void*) (td+threadIndex)) != 0)
	isgood = false;
    if (isgood)
      prolong_thread<dataType>(&tdbase);  // run thread 0

    if (isgood) {
      // Wait for threads to return
      for (threadIndex = 1; threadIndex < n_threads; threadIndex++)
	pthread_join(thread[threadIndex],NULL);
    }

    // Clean up
    for (threadIndex = 1; threadIndex < n_threads; threadIndex++) {
      delete[] td[threadIndex].dx;
      delete[] td[threadIndex].x;
    }
    delete[] thread;
    delete[] td;
  }
  delete[] tdbase.dx;
  delete[] tdbase.x;

  delete[] Asz;
  delete[] Aoutsz;
  delete[] Askip;
  delete[] Aoutskip;
  delete[] expandsz;
  delete[] coef;
  delete[] coefList;
  delete[] memoffset;
  delete[] memoffsetList;
  delete[] n;
  delete[] expandszCum;

  #ifdef MAIN
  CALLGRIND_STOP_INSTRUMENTATION;
  #endif

  return isgood;
}

bool prolong_work_double(const double *A,int n_dims,const int *Asz,double *Aout,const int *Aoutsz,int n_threads)
{
  return prolong_work<double>(A,n_dims,Asz,Aout,Aoutsz,n_threads);
}

bool prolong_work_float(const float *A,int n_dims,const int *Asz,float *Aout,const int *Aoutsz,int n_threads)
{
  return prolong_work<float>(A,n_dims,Asz,Aout,Aoutsz,n_threads);
}


#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 3;
  const int n_outputs = 1;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "A",
    "Aoutsz",
    "n_threads"
  };
  const char *output_names[] = {
    "Aout"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables A, Aoutsz, and n_threads\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  printf("Output file just before calling mexfcn: %s\n",fileout);
  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray**) input);

  // Save the outputs
  printf("Output file just before save: %s\n",fileout);
  mat_save_variables(fileout,output_names,n_outputs,output);

  // Clear the allocated variables
  /*
  for (int i = 0; i < n_inputs; i++)
    mxFree(input[i]);
  for (int i = 0; i < n_outputs; i++)
    mxFree(output[i]);
  */

  return EXIT_SUCCESS;
}
#endif
