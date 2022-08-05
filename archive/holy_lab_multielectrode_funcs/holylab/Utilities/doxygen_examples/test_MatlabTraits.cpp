#include <iostream>
#include "mat.h"
#include "MatlabTraits.h"
#include <error.h>
#include <errno.h>

// Compile this with
//  mex -f $MATLABINSTALLDIR/bin/matopts.sh -I.. test_MatlabTraits.cpp -o test_MatlabTraits

#define mexPrintf printf
#define mexErrMsgIdAndTxt(a, b, ...) error (1, 0, b, ##__VA_ARGS__)

template <typename T>
void myFunction(const mxArray *pMxInput,mxArray *pMxOutput)
 {
   // Check the input
   if (!MatlabTraits<T>::IsA(pMxInput))
     mexErrMsgIdAndTxt("mymexfile:typeError","Input is not of type %s %s.",MatlabTraits<T>::ComplexName(),MatlabTraits<T>::RealName());
   // Prepare the output
   mwSize szout[2] = {2, 2};
   pMxOutput = MatlabTraits<T>::mxAllocator(2,szout);
   if (pMxOutput == NULL)
     mexErrMsgIdAndTxt("mymexfile:allocationError","Failed to allocate output.");
   // Do whatever processing you intend. We have two cases here
     // because the storage format of Matlab complex types is
     // different from std::complex; for complex data you may have to do
     // extra work
   if (!MatlabTraits<T>::IsComplex) {
     // Get your hands on the input & output data
     const T* pInput = (const T*) mxGetData(pMxInput);
     T* pOutput = (T*) mxGetData(pMxOutput);
     // Now do something to process the data in pInput and stuff the
       // results into pOutput
   } else {
     // Complex case; here you may need to do more "munging" of the
       // data, but that will be problem-specific
   }
   mexPrintf("Function finished without errors\n");
 }

using namespace std;
int main()
{
  cout << "float = " << MatlabTraits<float>::RealName() << endl;
  typedef std::complex<int16_t> Tcomplex;
  if (MatlabTraits<Tcomplex>::IsComplex)
    cout << "It's complex" << endl;
  else
    cout << "It's real" << endl;
  cout << "complexity string is " << MatlabTraits<Tcomplex>::ComplexName() << endl;
  cout << "complex data type is " << MatlabTraits<Tcomplex>::RealName() << endl;

  mwSize sz[2] = {5, 3};
  mxArray *pMxInput = mxCreateNumericArray(2, sz, mxSINGLE_CLASS, mxCOMPLEX);
  mxArray *pMxOutput;
  
  cout << "\nThis next should succeed:" << endl;
  myFunction<std::complex<float> >(pMxInput,pMxOutput);
  cout << "\nThis next should fail:" << endl;
  myFunction<std::complex<double> >(pMxInput,pMxOutput);
}
