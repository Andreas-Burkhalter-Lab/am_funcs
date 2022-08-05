/*
    comedi_get_range
    Copyright (C) 2002 Timothy E. Holy <holy@pcg.wustl.edu>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation, version 2.1
    of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA.
*/

#include "mex.h"
#include <comedilib.h>
#include <comedipersist.h>
#include "commat.h"

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int buflen;
  char *devname;
  unsigned int subdevice;
  unsigned int channel;
  unsigned int range;
  comedi_t *dev;
  comedi_range *crange;
  const mxArray *curarg;
  int indx;

  if (nrhs < 4)
    mexErrMsgTxt("comedi_get_range: requires 4 input arguments");

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
    mexErrMsgTxt("comedi_get_range: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);
  /* Channel */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_get_range: channel must be a real scalar (nonneg integer)");
  channel = (unsigned int) mxGetScalar(curarg);
  /* Range */
  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_get_range: range must be a real scalar (nonneg integer)");
  range = (unsigned int) mxGetScalar(curarg);

  /* Get the range */
  crange = comedi_get_range(dev,subdevice,channel,range);
  if (crange == NULL) {
    comedi_perror("comedi_get_range");
    mexErrMsgTxt("comedi_get_range: error getting range from device");
  }

  /* Return data to matlab */
  plhs[0] = mxCreateDoubleMatrix(2,1,mxREAL);
  *mxGetPr(plhs[0]) = crange->min;
  *(mxGetPr(plhs[0])+1) = crange->max;
  if (nlhs > 1) {
    indx = intmatch(crange->unit,unitvals,n_unit);
    if (indx < 0)
      mexErrMsgTxt("comedi_get_range: error on the unit setting");
    plhs[1] = mxCreateString(unitnames[indx]);
  }
} 
