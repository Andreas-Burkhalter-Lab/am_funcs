// Define MAIN while compiling to create a stand-alone command-line
// program. See Makefile.
#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#include <stdio.h>
//#include <callgrind.h>   // for profiling
#else
#include "mex.h"
#endif

#include <stdint.h>
#include <string.h>
#include <arpa/inet.h>

#if defined __LITTLE_ENDIAN
#define swapbytes(x) ntohl(x)  /* use system byteswap (ntohl is a noop on bigendian systems) */
#else
uint32_t swapbytes(uint32_t x) {
  return ((x & 0x000000ffU) << 24) |
         ((x & 0x0000ff00U) <<  8) |
         ((x & 0x00ff0000U) >>  8) |
         ((x & 0xff000000U) >> 24);
}
#endif

typedef union {
   uint32_t u32;
   float flt;
} U32;

static const unsigned char *b64_tbl = (const unsigned char*) "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static const unsigned char b64_pad = '=';


void b64_decode_mio ( char *dest,  char *src )
{
   char *temp;

  temp = dest;

  while (*src)
    {
      int register a;
      int register b;
      /*
      int32_t a;
      int32_t b;
      */
      int t1,t2,t3,t4;

      t1 = src[0];
      t2 = src[1];
      t3 = src[2];
      t4 = src[3];

     
      if( t1 > 96 )		// [a-z]
	a = (t1 - 71);
      else if( t1 > 64 )		// [A-Z]
	a = (t1 - 65);
      else if( t1 > 47 ) {		// [0-9], or '='
         if (t1 == 61 )	{	// if == '='
         	return;
      }
	a = (t1 + 4);
      } else if( t1 == 43 )
	a = 62;
      else				// src[0] == '/'
	a = 63;     


      if( t2 > 96 )		// [a-z]
	b = (t2 - 71);
      else if( t2 > 64 )		// [A-Z]
	b = (t2 - 65);
      else if( t2 > 47 )		// [0-9]
	b = (t2 + 4);
      else if( t2 == 43 )
	b = 62;
      else				// src[1] == '/'
	b = 63;     
    
      *temp++ = ( a << 2) | ( b >> 4);
     
      if( t3 > 96 )		// [a-z]
	a = (t3 - 71);
      else if( t3 > 64 )		// [A-Z]
	a = (t3 - 65);
      else if( t3 > 47 ) {		// [0-9], or '='
         if (t3 == 61) {
            return;
         }
	a = (t3 + 4);
      } else if( t3 == 43 )
	a = 62;
      else				// src[2] == '/'
	a = 63;     


      *temp++ = ( b << 4) | ( a >> 2);

      if (t4 == 61)
	return;

      if( t4 > 96 )		// [a-z]
	b = (t4 - 71);
      else if( t4 > 64 )		// [A-Z]
	b = (t4 - 65);
      else if( t4 > 47 )		// [0-9]
	b = (t4 + 4);
      else if( t4 == 43 )
	b = 62;
      else				// src[3] == '/'
	b = 63;    

      *temp++ = ( a << 6) | ( b );

      src += 4;
    }
}


void convert_work(float *output,char *pdecoded,char *pdata,int n_peaks)
{
  int i;
  U32 tmp;

  // Convert from network format
  b64_decode_mio(pdecoded,pdata);


  // Do the endian conversion and load into the result
  for (i = 0; i < 2*n_peaks; i++) {
    tmp.u32 = swapbytes(((uint32_t *) pdecoded)[i]);
    output[i] = tmp.flt;
  }
}




void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mxArray *curarg;
  char *pdata;
  char *pdecoded;
  int len;
  int status;
  int n_peaks;
  float *output;
  int buffer = 0;

  if (nrhs != 2)
    mexErrMsgTxt("mzxml_decode: requires two inputs");
  if (nlhs != 1)
    mexErrMsgTxt("mzxml_decode: requires one output");

  // Parse the input
  // The character string
  curarg = prhs[0];
  if (!mxIsUint8(curarg))
    mexErrMsgTxt("mzxml_decode: requires uint8 input");
  len = mxGetNumberOfElements(curarg);
  if (4*(len/4) != len)
    mexErrMsgTxt("mzxml_decode: input must have length that is a multiple of 4");
  pdata = (char *) mxMalloc((len+buffer+1)*sizeof(char));
  memcpy(pdata,mxGetData(curarg),len);
  pdata[len] = '\0';
  /*
  if (status) {
    mexPrintf("len = %d, m %d, n %d\n",len,mxGetM(curarg),mxGetN(curarg));
    mexErrMsgTxt("mzxml_decode: unknown problem in copying string");
  }
  */

  // The # of peaks
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || mxGetNumberOfElements(curarg)!=1)
    mexErrMsgTxt("mzxml_decode: n_peaks must be a real scalar");
  n_peaks = (int) mxGetScalar(curarg);

  // Allocate space for temporary storage
  pdecoded = (char *) mxMalloc(8*n_peaks+1);
  pdecoded[8*n_peaks] = '\0';

  // Allocate space for the output
  plhs[0] = mxCreateNumericMatrix(2,n_peaks,mxSINGLE_CLASS,mxREAL);
  output = (float *) mxGetData(plhs[0]);

  convert_work(output,pdecoded,pdata,n_peaks);

  mxFree(pdecoded);
  mxFree(pdata);
}


#ifdef MAIN
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 2;
  const int n_outputs = 1;
  mxArray *input[n_inputs];
  mxArray *output[n_outputs];
  const char *input_names[] = {
    "str",
    "n_peaks"
  };
  const char *output_names[] = {
    "peakdata"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables str and n_peaks\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];
  printf("Input file: %s\nOutput file: %s\n",filein,fileout);

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);

  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs,(const mxArray**) input);

  // Save the outputs
  mat_save_variables(fileout,output_names,n_outputs,output);

  return EXIT_SUCCESS;
}
#endif
