
#include "matlab_arg.h"

namespace matlab_arg {

#include <string>

using namespace std;

#include <assert.h>

//------------------------------------------------------------------

//usage sample:
//  getDoubleArg(prhs[0]);
double getDoubleArg(const mxArray * arg)
{
   /* The input must be a noncomplex scalar double.*/
   int m = mxGetM(arg);
   int n = mxGetN(arg);
   if( !mxIsDouble(arg) || mxIsComplex(arg) || !(m==1 && n==1) ) {
      mexErrMsgTxt("Input must be a noncomplex scalar double.");
   }
  
   return *mxGetPr(arg); //or:  return mxGetScalar(arg);
}

//------------------------------------------------------------------
//usage sample:
//plhs[0]=setDoubleArg(10);
//or:  setDoubleArg(plhs[0], 10);
mxArray * setDoubleArg(double value)
{
   mxArray * result = mxCreateDoubleMatrix(1, 1, mxREAL); 
           //mxCreateDoubleMatrix() creates mxArray, like class string in c++
   *mxGetPr(result) = value;       
           //mxGetPr() returns internal data buf of a mxArray, like string::data() in c++

   return result; 
}
//---------------------------------------------------------
mxArray * & setDoubleArg(mxArray * & rArg, double value)
{
   rArg = mxCreateDoubleMatrix(1, 1, mxREAL);
   *mxGetPr(rArg) = value;

   return rArg; 
}

//------------------------------------------------------------------

string getStringArg(const mxArray * arg)
{
    /* input must be a string */
    if( mxIsChar(arg) != 1)  mexErrMsgTxt("Input must be a string.");

    /* input must be a row vector or column vector */
    int m=mxGetM(arg); int n=mxGetN(arg);
    if(m!=1 && n!=1)  mexErrMsgTxt("Input must be a row or column vector.");
    
    /* get the length of the input string */
    int buflen = m * n + 1;

    /* allocate memory for input string */
    char * input_buf=(char*)mxCalloc(buflen, sizeof(char)); //if fail, MEX file will be terminated
    assert(input_buf);
    
    /* copy the string data from arg into a C string input_buf.
     * If the string array contains several rows, they are copied,
     * one column at a time, into one long string array.
     */
    int status = mxGetString(arg, input_buf, buflen);
    if(status != 0) mexWarnMsgTxt("Not enough space. String is truncated.");
    
    string result(input_buf);

    mxFree(input_buf); //this can be skipped. Matlab will auto reclaim the mem

    return result;
}

//------------------------------------------------------------------
mxArray * setStringArg(string& value)
{
   mxArray * result = mxCreateString(value.c_str());

   return result; 
}
//---------------------------------------------------------
mxArray * & setStringArg(mxArray * & rArg, string& value)
{
   rArg = mxCreateString(value.c_str());

   return rArg; 
}

//------------------------------------------------------------------

TDoubleMatrixParam getDoubleMatrixArg(const mxArray * arg)
{
   /* The input must be a noncomplex double matrix.*/
   int m = mxGetM(arg);
   int n = mxGetN(arg);
   if( !mxIsDouble(arg) || mxIsComplex(arg) ) {
      mexErrMsgTxt("Input must be a noncomplex double matrix.");
   }

   TDoubleMatrixParam result;
   result.p = mxGetPr(arg);
   result.m = m;
   result.n = n;
  
   return result; 
}

//------------------------------------------------------------------

} //namespace, matlab_arg
