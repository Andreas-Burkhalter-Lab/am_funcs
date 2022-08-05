// Define MAIN while compiling to create a stand-alone command-line
// program. See Makefile.
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
#include <callgrind.h>   // for profiling
#else
#include "mex.h"
#endif

// Use the following for performance optimization
// (This must be before the include imiterators.cxx statement)
#define MAXDIMS 4

#include <string.h>  // for memcpy and imiterators
#include <pthread.h>
#include <unistd.h>  // for # of CPUs
#include "imiterators.cxx"
#include "image_utilities.h"


/* Syntax:
 *   Ar = array_restrict(A)
 *   Ar = array_restrict(A,dimFlag)
 *   Ar = array_restrict(A,dimFlag,n_threads)
 */

// forward declarations
int restrict_work_double(const double *A,int n_dims,const int *Asz,const bool *restrict_dim,double *Aout,int n_threads);
int restrict_work_float(const float *A,int n_dims,const int *Asz,const bool *restrict_dim,float *Aout,int n_threads);


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
  const int *Asz;  // the size of A
  int *Aoutsz;     // the size of Aout
  bool *restrict_dim;  // whether to restrict a given dimension
  const mxArray *curarg;
  const mxLogical* dimFlag;
  bool isdouble;
  int n_threads,n_cpus;

  if (nrhs < 1 || nrhs > 3)
    mexErrMsgTxt("array_restrict: requires one to three inputs");
  if (nlhs != 1)
    mexErrMsgTxt("array_restrict: requires one output");

  // Parse the input
  // A
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg))
    mexErrMsgTxt("array_restrict: A must be a real array");
  if (mxIsSingle(curarg))
    isdouble = false;
  else if (mxIsDouble(curarg))
    isdouble = true;
  else
    mexErrMsgTxt("array_restrict: A must be single or double");
  if (mxIsEmpty(curarg))
    mexErrMsgTxt("array_restrict: does not work on empty arrays");

  A = mxGetData(curarg);
  n_dims = mxGetNumberOfDimensions(curarg);
#ifdef MAXDIMS
  if (n_dims > MAXDIMS)
    mexErrMsgTxt("array_restrict: too many dimensions in input");
#endif
  Asz = mxGetDimensions(curarg);
  Aoutsz = (int *) mxMalloc(n_dims*sizeof(int));
  restrict_dim = (bool *) mxMalloc(n_dims*sizeof(bool));

  // dimFlag
  if (nrhs > 1) {
    curarg = prhs[1];
    if (!mxIsLogical(curarg))
      mexErrMsgTxt("array_restrict: dimFlag must be a real logical vector");
    if (mxGetNumberOfElements(curarg) != n_dims)
      mexErrMsgTxt("array_restrict: number of dimensions in A must match length of dimFlag");
    dimFlag = mxGetLogicals(curarg);
    // Copy it so can be modified (e.g., so we protect non-reduceable dimensions)
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
      restrict_dim[dimIndex] = (dimFlag[dimIndex] > 0);
  } else
    for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
      restrict_dim[dimIndex] = true;
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (Asz[dimIndex] < 2)
      restrict_dim[dimIndex] = false;  // protect non-reduceable dimensions


  // n_threads
  n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
  //n_threads = n_cpus;
  n_threads = 1;  // Currently not much benefit to multiple threads
  if (nrhs > 2) {
    curarg = prhs[2];
    if (mxGetNumberOfElements(curarg) != 1)
      mexErrMsgTxt("array_restrict: n_threads must be a scalar");
    n_threads = (int) mxGetScalar(curarg);
    if (n_threads < 1)
      n_threads = 1;
    if (n_threads > n_cpus)
      n_threads = n_cpus;
  }

  // Calculate the size of the output
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (restrict_dim[dimIndex]) {
      if (Asz[dimIndex] % 2)
	Aoutsz[dimIndex] = (Asz[dimIndex]+1)/2;  // odd case
      else {                           // even case
	if (Asz[dimIndex] == 2)
	  Aoutsz[dimIndex] = 1;
	else
	  Aoutsz[dimIndex] = Asz[dimIndex]/2 + 1;
      }
    }
    else
      Aoutsz[dimIndex] = Asz[dimIndex];

  // Set up the output and do the work
  if (isdouble) {
    plhs[0] = mxCreateNumericArray(n_dims,Aoutsz,
				   mxDOUBLE_CLASS,mxREAL);
    Aout = mxGetData(plhs[0]);
    if (!restrict_work_double((double *) A,n_dims,Asz,restrict_dim,(double *) Aout,n_threads))
      mexErrMsgTxt("Error creating threads!");
  } else {
    plhs[0] = mxCreateNumericArray(n_dims,Aoutsz,
				   mxSINGLE_CLASS,mxREAL);
    Aout = mxGetData(plhs[0]);
    if (!restrict_work_float((float *) A,n_dims,Asz,restrict_dim,(float *) Aout,n_threads))
      mexErrMsgTxt("Error creating threads!");
  }

  return;
}

// A structure for passing data to pthreads
template <class Tdata>
struct thread_data_type {
  // Frequently-written data
  pixIterator outI;
  // Rarely- or never-written data
  int n_dims;          // number of dimensions
  int kern_numel;      // total number of "pixels" in the kernel
  Tdata *kern_vector;  // the amplitude assigned to each neighbor pixel
  int *source_offsets; // the memory displacement to find neighbors in source
  int *coord_offsets;  // the coordinates of each element of the kernel
  const int *Asz;
  const long *Askip;   // the pointer increment for each coordinate of A
  const Tdata *Ain;
  Tdata *Aout;         // the array is frequently-written, but the ptr isn't
  const bool *restrict_dim;

  thread_data_type<Tdata>() {n_dims = 0; kern_numel = 0; kern_vector = NULL; source_offsets = NULL; coord_offsets = NULL;}
  //thread_data_type<Tdata>(thread_data_type<Tdata> &td);
  ~thread_data_type<Tdata>();

  //void copy(thread_data_type<Tdata> &td);
};

/*
template <class Tdata>
thread_data_type<Tdata>::thread_data_type(thread_data_type<Tdata> &td)
{
  copy(td);
}

template <class Tdata>
void thread_data_type<Tdata>::copy(thread_data_type<Tdata> &td)
{
  n_dims = td.n_dims;
  kern_numel = td.kern_numel;
  Ain = td.Ain;
  Aout = td.Aout;
  kern_vector = new Tdata[kern_numel];
  source_offsets = new int[kern_numel];
  coord_offsets = new int[kern_numel*n_dims];
  restrict_dim = new bool[n_dims];
  
  memcpy(kern_vector,td.kern_vector,kern_numel*sizeof(Tdata));
  memcpy(source_offsets,td.source_offsets,kern_numel*sizeof(int));
  memcpy(coord_offsets,td.coord_offsets,kern_numel*n_dims*sizeof(int));
  memcpy(restrict_dim,td.restrict_dim,n_dims*sizeof(bool));

  outI.copy(td.outI);
}
*/

template <class Tdata>
thread_data_type<Tdata>::~thread_data_type()
{
  ;
}
    

// Wrapper for multithreading
template <class Tdata>
void *restrict_thread_wrapper(void *p)
{
  thread_data_type<Tdata> *td;

  td = (thread_data_type<Tdata>*) p;
  restrict_thread(td);
}


void print_pixI_coords(const pixIterator &pI)
{
  for (int dimIndex = 0; dimIndex < pI.nDims(); dimIndex++)
    printf("%d ",pI.coord(dimIndex));
  printf(" (current %ld)\n",long(pI));
}

template <class Tdata>
struct rtvars {
  int dimIndex;
  Tdata tmp;
  bool is_inside;
  int c;
  const Tdata *Aincur;
  Tdata *kern_vector;
  int *source_offset;
  int *coord_offset;
#ifdef MAX_DIMS
  int coord_in[MAX_DIMS];
#else
  int *coord_in;
#endif

  rtvars(int n_dims) {
#ifndef MAX_DIMS
    coord_in = new int[n_dims];
#endif
  }
  ~rtvars() {
#ifndef MAX_DIMS
    delete[] coord_in;
#endif
  }
} __attribute__((aligned(64)));

// This is the main calculation
template <class Tdata>
void restrict_thread(thread_data_type<Tdata> *td)
{
  // Frequently-written variables
  rtvars<Tdata> rtv(td->n_dims);
  // Rarely-written/read-only variables
  const Tdata *Ainend;
  Tdata *kern_end;
  int n_dims;
  int inc;

  kern_end = td->kern_vector+td->kern_numel;
  n_dims = td->n_dims;
  Ainend = td->Ain + td->Asz[n_dims-1]*td->Askip[n_dims-1];
  /*
  printf("n_dims: %d\n",n_dims);
  printf("numel(A) = %d\n",Ainend-td->Ain);
  printf("restrict_dim: ");
  for (rtv.dimIndex = 0; rtv.dimIndex < n_dims; rtv.dimIndex++)
    printf("%d ",td->restrict_dim[rtv.dimIndex]);
  printf("\n");
  */
  inc = td->restrict_dim[0]+1;
  for (; !td->outI.at_end(); td->outI++, rtv.Aincur+=inc, rtv.coord_in[0]+=inc) {
    if (td->outI.coord(0) == 0) {
      // We're at the beginning of a column; initialize the input coordinates
      rtv.Aincur = td->Ain;
      for (rtv.dimIndex = 0; rtv.dimIndex < n_dims; rtv.dimIndex++) {
	rtv.coord_in[rtv.dimIndex] = td->outI.coord(rtv.dimIndex) * (1 + td->restrict_dim[rtv.dimIndex]);
	rtv.Aincur += rtv.coord_in[rtv.dimIndex] * td->Askip[rtv.dimIndex];
      }
    }
    /*
    printf("outI coords: ");
    print_pixI_coords(td->outI);
    printf("inI coords: ");
    for (rtv.dimIndex = 0; rtv.dimIndex < n_dims; rtv.dimIndex++)
      printf("%d ",rtv.coord_in[rtv.dimIndex]);
    printf("(current %d)\n",rtv.Aincur-td->Ain);
    */
    if (!td->outI.on_edge(td->restrict_dim)) {
      // Not on an edge, can do the calculation without checks
      //printf("Not on edge\n");
      rtv.tmp = 0;
      for (rtv.source_offset = td->source_offsets, rtv.kern_vector = td->kern_vector;
	   rtv.kern_vector != kern_end; rtv.source_offset++, rtv.kern_vector++) {
	/*
	if (rtv.Aincur+*source_offset < td->Ain || rtv.Aincur+*rtv.source_offset >= Ainend) {
	  printf("Bad memory access in input: ");
	  for (rtv.dimIndex = 0; rtv.dimIndex < n_dims; rtv.dimIndex++)
	    printf("%d ",rtv.coord_in[rtv.dimIndex]);
	  printf("Current offset: %d, current kernentry: %g\n",*rtv.source_offset,*rtv.kern_vector);
	}
	*/
	rtv.tmp += *(rtv.Aincur+*rtv.source_offset) * *rtv.kern_vector;
	//printf("%d/%g/%g ",*rtv.source_offset,*(rtv.Aincur+*rtv.source_offset),*rtv.kern_vector);
      }
      td->Aout[td->outI] = rtv.tmp;
    } else {
      // We're on an edge, must check each element of the kernel to
      // see if it's within the coordinate bounds
      rtv.tmp = 0;
      for (rtv.source_offset = td->source_offsets, rtv.kern_vector = td->kern_vector, rtv.coord_offset = td->coord_offsets;
	   rtv.kern_vector != kern_end; rtv.source_offset++, rtv.kern_vector++, rtv.coord_offset+=n_dims) {
	rtv.is_inside = true;
	for (rtv.dimIndex = 0; rtv.dimIndex < n_dims; rtv.dimIndex++) {
	  rtv.c = rtv.coord_in[rtv.dimIndex] + rtv.coord_offset[rtv.dimIndex];
	  if (rtv.c < 0 || rtv.c >= td->Asz[rtv.dimIndex]) {
	    rtv.is_inside = false;
	    break;
	  }
	}
	if (rtv.is_inside) {
	  rtv.tmp += *(rtv.Aincur+*rtv.source_offset) * *rtv.kern_vector;
	  //printf("Inside\n");
	}
	//else
	//  printf("Not inside\n");
      }
      td->Aout[td->outI] = rtv.tmp;
    }
  }
}

template <class Tdata>
int restrict_work(const Tdata *A,int n_dims,const int *Asz,const bool* restrict_dim,Tdata *Aout, int n_threads)
{
  // This function initializes everything: it sets up the kernel, and
  // then configures the threads

  #ifdef MAIN
  CALLGRIND_START_INSTRUMENTATION;
  #endif

  // Eliminate unity dimensions
  int dimIndex;
  /*
  int *Asz;
  bool *restrict_dim;
  int *indxkeep;
  n_dims = skip_unity_dimensions_index(Asz_in,n_dims,indxkeep,n_dims);
  Asz = new int[n_dims];
  restrict_dim = new bool[n_dims];
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    Asz[dimIndex] = Asz_in[indxkeep[dimIndex]];
    restrict_dim[dimIndex] = restrict_dim_in[indxkeep[dimIndex]];
  }
  Asz = Asz_in;
  restrict_dim = restrict_dim_in;
  */
  /*
  printf("n_dims: %d\nAsz: ",n_dims);
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    printf("%d ",Asz[dimIndex]);
  printf("\n");
  */

  // Different kernel options depending on whether a dimension of A is
  // even, odd, or of size < 3
  const Tdata kernodd[3] = {0.25, 0.5, 0.25};
  const Tdata kerneven[4] = {0.125, 0.375, 0.375, 0.125};
  const Tdata kernunity[1] = {1.0};
  const Tdata kerntwo[3] = {0, 0.5, 0.5};

  const Tdata **kern_array = new const Tdata*[n_dims];
  int *kernsz = new int[n_dims];
  int *Aoutsz = new int[n_dims];
  int kern_numel = 0;
  Tdata tmp;

  // Calculate size of Aout (we do this again in case this function
  // gets called from a different wrapper...) and assign the
  // "kernel" along each dimension
  // (depending on restrict_dim and whether we're even/odd)
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    if (restrict_dim[dimIndex]) {
      if (Asz[dimIndex] % 2) {  // odd case
	Aoutsz[dimIndex] = (Asz[dimIndex]+1)/2;
	kernsz[dimIndex] = 3;
	kern_array[dimIndex] = kernodd;
      }
      else {  // even case
	if (Asz[dimIndex] == 2) {
	  Aoutsz[dimIndex] = 1;
	  kernsz[dimIndex] = 3;
	  kern_array[dimIndex] = kerntwo;
	} else {
	  Aoutsz[dimIndex] = Asz[dimIndex]/2 + 1;
	  kernsz[dimIndex] = 4;
	  kern_array[dimIndex] = kerneven;
	}
      }
    }
    else {
      Aoutsz[dimIndex] = Asz[dimIndex];
      kernsz[dimIndex] = 1;
      kern_array[dimIndex] = kernunity;
    }

  kern_numel = kernsz[0];
  for (dimIndex = 1; dimIndex < n_dims; dimIndex++)
    kern_numel = kern_numel * kernsz[dimIndex];
  
  /*
  printf("Kernel dimensions: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    printf("%d ",kernsz[dimIndex]);
  printf("\n");
  */
  /*
  printf("Ain dimensions: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    printf("%d ",Asz[dimIndex]);
  printf("\n");
  printf("Aout dimensions: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    printf("%d ",Aoutsz[dimIndex]);
  printf("\n");
  printf("restrict_dim: ");
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++)
    printf("%d ",restrict_dim[dimIndex]);
  printf("\n");
  */

  pixIterator kernI(kernsz,n_dims,false);

  // Compute the kernel as a product of all the contributions from
  // each dimension. Also assign "coordinate offsets" to each element
  // of the kernel, so that at the edges we can determine which
  // entries to keep and which to ignore.  The coordinate system is
  // defined so that entries with offsets 0 and higher are used at the
  // "left" corner; offsets less than 0 are used at the "right"
  // corner. (Interior points use all entries in the kernel.)
  int *coord_offset = new int[n_dims];
  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    coord_offset[dimIndex] = -(kernsz[dimIndex]/2); // 0 for unity, -1 for odd, -2 for even
  }
  long *Askip = new long[n_dims];
  calc_pixel_skip(Asz,n_dims,Askip);

  thread_data_type<Tdata> tdbase;
  tdbase.outI.initialize(Aoutsz,n_dims,false);
  tdbase.kern_vector = new Tdata[kern_numel];
  tdbase.source_offsets = new int[kern_numel];
  tdbase.coord_offsets = new int[kern_numel*n_dims];
  tdbase.Asz = Asz;
  tdbase.Askip = Askip;
  tdbase.Ain = A;
  tdbase.Aout = Aout;
  tdbase.restrict_dim = restrict_dim;
  tdbase.n_dims = n_dims;
  tdbase.kern_numel = kern_numel;

  int source_offset;
  int indx;
  for (; !kernI.at_end(); kernI++) {
    //printf("kernel coords: ");
    //print_pixI_coords(kernI);
    tmp = kern_array[0][kernI.coord(0)];
    for (dimIndex = 1; dimIndex < n_dims; dimIndex++)
      tmp *= kern_array[dimIndex][kernI.coord(dimIndex)];
    //printf("Kernel value: %g\n",tmp);
    tdbase.kern_vector[kernI] = tmp;
    for (source_offset = 0, dimIndex = 0; dimIndex < n_dims; dimIndex++) {
      indx = kernI.coord(dimIndex)+coord_offset[dimIndex];
      source_offset += indx * Askip[dimIndex];
      tdbase.coord_offsets[kernI*n_dims+dimIndex] = indx;
    }
    tdbase.source_offsets[kernI] = source_offset;
  }

  // Do what cleaning up we can at this point
  delete[] kern_array;
  delete[] kernsz;
  delete[] coord_offset;

  // Partition the work into the different threads
  // Divide the work by "columns". Don't create any more threads than
  // there are columns.
  long out_numcol = 1;  // the # of "columns" in the output
  for (dimIndex = 1; dimIndex < n_dims; dimIndex++)
    out_numcol *= Aoutsz[dimIndex];
  if (out_numcol < n_threads)
    n_threads = out_numcol;

  if (n_threads <= 1) {
    restrict_thread(&tdbase);
  }
  else {
    //printf("%d threads being used.\n",n_threads);
    float thread_ncol = ((float) out_numcol) / n_threads;
    thread_data_type<Tdata> *td;
    td = new thread_data_type<Tdata>[n_threads];
    int threadIndex;
    long colstart,colend;
    for (threadIndex = 0; threadIndex < n_threads; threadIndex++) {
      td[threadIndex] = tdbase;
      colstart = long(threadIndex*thread_ncol);
      colend = long((threadIndex+1)*thread_ncol);
      //printf("colstart %d, colend %d\n",colstart,colend);
      td[threadIndex].outI.setrange(colstart*Aoutsz[0],colend*Aoutsz[0]);
    }
    
    // Launch the threads. We don't create a thread 0, we run
    // that one in this thread
    pthread_t *thread;
    thread = new pthread_t[n_threads];
    for (threadIndex = 1; threadIndex < n_threads; threadIndex++)
      if (pthread_create(thread+threadIndex,NULL,restrict_thread_wrapper<Tdata>,(void*) (td+threadIndex)) != 0) {
	delete[] Aoutsz;
	delete[] Askip;
	delete[] tdbase.kern_vector;
	delete[] tdbase.source_offsets;
	delete[] tdbase.coord_offsets;
	delete[] thread;
	delete[] td;
	return 0; // error creating thread
      }
    restrict_thread_wrapper<Tdata>((void *) td);  // run thread 0

    // Wait for threads to return
    for (threadIndex = 1; threadIndex < n_threads; threadIndex++)
      pthread_join(thread[threadIndex],NULL);

    // Clean up
    delete[] thread;
    delete[] td;
  }

  delete[] Aoutsz;
  delete[] Askip;
  delete[] tdbase.kern_vector;
  delete[] tdbase.source_offsets;
  delete[] tdbase.coord_offsets;

  #ifdef MAIN
  CALLGRIND_STOP_INSTRUMENTATION;
  #endif

  return 1;
}

int restrict_work_double(const double *A,int n_dims,const int *Asz,const bool *restrict_dim,double *Aout,int n_threads)
{
  return restrict_work<double>(A,n_dims,Asz,restrict_dim,Aout,n_threads);
}

int restrict_work_float(const float *A,int n_dims,const int *Asz,const bool *restrict_dim,float *Aout,int n_threads)
{
  return restrict_work<float>(A,n_dims,Asz,restrict_dim,Aout,n_threads);
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
    "dimFlag",
    "n_threads"
  };
  const char *output_names[] = {
    "Aout"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables A, dimFlag, and n_threads\n",argv[0]);
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

  return EXIT_SUCCESS;
}
#endif
