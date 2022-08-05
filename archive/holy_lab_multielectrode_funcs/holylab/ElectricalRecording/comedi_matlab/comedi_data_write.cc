/*
    comedi_data_write
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
  unsigned int aref;
  int indx;
  lsampl_t data;
  comedi_t *dev;
  const mxArray *curarg;
  int nwritten;

  if (nrhs < 6)
    mexErrMsgTxt("comedi_data_write: requires 6 input arguments");

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
    mexErrMsgTxt("comedi_data_write: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);
  /* Channel */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_data_write: channel must be a real scalar (nonneg integer)");
  channel = (unsigned int) mxGetScalar(curarg);
  /* Range */
  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_data_write: range must be a real scalar (nonneg integer)");
  range = (unsigned int) mxGetScalar(curarg);
  /* Aref */
  curarg = prhs[4];
  aref = arefvals[getvalfromstr(curarg,arefnames,n_aref)];
  /* Data */
  curarg = prhs[5];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_data_write: data must be a real scalar (nonneg integer) from 0 to 4095");
  data = (lsampl_t) mxGetScalar(curarg);

  /* mexPrintf("Writing value %d to channel %d\n",data,channel); */

  /* Write the data point */
  if ((nwritten = comedi_data_write(dev,subdevice,channel,range,aref,data)) < 0) {
    comedi_perror("comedi_data_write");
    mexErrMsgTxt("comedi_data_write: error writing data to device");
  }
} 
