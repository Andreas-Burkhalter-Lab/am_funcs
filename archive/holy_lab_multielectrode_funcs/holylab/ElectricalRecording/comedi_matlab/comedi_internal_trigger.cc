/*
    who: comedi_interal_trigger.cc
    what: wrap comedi_internal_trigger() in <comedi>/demo/ao_waveform.cc
    when last updated: Oct 2002
    
    modified by Jason Z.S. Guo from comedi_get_n_ranges.cc

*/

#include <string.h>

#include "mex.h"
#include <comedilib.h>
#include <comedipersist.h>
#include "commat.h"

//copy from <comedi>/demo/ao_waveform.cc
static int comedi_internal_trigger(comedi_t *dev, unsigned int subd, unsigned int trignum)
{
	comedi_insn insn;
	lsampl_t data[1];

	memset(&insn, 0, sizeof(comedi_insn));
	insn.insn = INSN_INTTRIG;
	insn.subdev = subd;
	insn.data = data;
	insn.n = 1;

	data[0] = trignum;

	return comedi_do_insn(dev, &insn);
}//comedi_internal_trigge


void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int buflen;
  char *devname;
  unsigned int subdevice;
  unsigned int trigNum;
  unsigned int range;
  comedi_t *dev;
  const mxArray *curarg;
  int result;

  if (nrhs != 3)
    mexErrMsgTxt("comedi_interal_trigger: requires 3 input arguments");

  /* Device */
  curarg = prhs[0];
  buflen = matstr_len(curarg);
  devname = (char*) mxCalloc(buflen,sizeof(char));
  matstr_to_cstr(curarg,devname,buflen);
  dev = get_comedi_device_safe(devname);
  mxFree(devname);

  /* Subdevice */
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_internal_trigger: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);

  /* trigNum */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_internal_trigger: trigNum must be a real scalar (nonneg integer)");
  trigNum = (unsigned int) mxGetScalar(curarg);

  result = comedi_internal_trigger(dev,subdevice,trigNum);
  if (result < 0) {
    comedi_perror("comedi_internal_trigger");
    mexErrMsgTxt("comedi_internal_trigger: error on internal_trigger");
  }

  /* Return data to matlab */
  plhs[0] = mxCreateScalarDouble(result);
} 


