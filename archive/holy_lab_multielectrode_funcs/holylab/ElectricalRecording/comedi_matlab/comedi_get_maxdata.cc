//compile:
// mex -L /usr/local/matlab13/extern/lib/glnx86 comedi_get_maxdata.cc matlab_arg.cpp -lcomedi -lcomedipersist

#include "mex.h"
#include <comedilib.h>
#include "comedipersist.h"
//#include "commat.h"
#include "matlab_arg.h"

//get_comedi_device_safe(): this code is copied from Tim's commat_util.cc
comedi_t* get_comedi_device_safe(const char *devname)
{
  if (!comedi_device_isopen(devname)) {
    mexPrintf("Error using device %s:\n",devname);
    mexErrMsgTxt("Device is not open");
  }
  return get_comedi_device(devname);
}

//IsScalar(): this code is copied from Tim's commat_util.cc
bool IsScalar(const mxArray *m)
{
  return (mxGetNumberOfElements(m) == 1);
}


void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  if(nrhs < 3)
     mexErrMsgTxt("comedi_get_maxdata(): requires 3 input arguments");

  // Device
  comedi_t * dev = get_comedi_device_safe(matlab_arg::getStringArg(prhs[0]).c_str());
  // Subdevice 
  if(!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || !IsScalar(prhs[1]))
     mexErrMsgTxt("comedi_get_maxdata(): subdevice must be a real scalar (nonneg integer)");
  unsigned int subdevice = (unsigned int) matlab_arg::getDoubleArg(prhs[1]);
  // channel
  if(!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || !IsScalar(prhs[2]))
     mexErrMsgTxt("comedi_get_maxdata(): channel must be a real scalar (nonneg integer)");
  unsigned int channel = (unsigned int) matlab_arg::getDoubleArg(prhs[2]);

  int result = comedi_get_maxdata(dev,subdevice, channel);
  if (result < 0) {
    comedi_perror("comedi_get_maxdata()");
    mexErrMsgTxt("comedi_get_maxdata(): error occured");
  }

  // Return data to matlab
  matlab_arg::setDoubleArg(plhs[0], result);
} 
