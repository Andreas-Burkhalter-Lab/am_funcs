/*
    comedi_dio_bitfield
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
  int subdevtype;
  unsigned int write_mask;
  unsigned int bits;
  int indx;
  comedi_t *dev;
  const mxArray *curarg;

  if (nrhs < 4)
    mexErrMsgTxt("comedi_dio_bitfield: requires 4 input arguments");

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
    mexErrMsgTxt("comedi_dio_bitfield: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);
  /* Write_mask */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_dio_bitfield: write_mask must be a real scalar (nonneg integer)");
  write_mask = (unsigned int) mxGetScalar(curarg);
  /* Bits */
  curarg = prhs[3];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_dio_bitfield: bits must be a real scalar (nonneg integer)");
  bits = (unsigned int) mxGetScalar(curarg);

  /* Check that this is a DIO device */
  subdevtype = comedi_get_subdevice_type(dev,subdevice);
  if (subdevtype != COMEDI_SUBD_DIO && subdevtype != COMEDI_SUBD_DI)
    mexErrMsgTxt("Subdevice is not of type DI or DIO");
     
  /* Execute the function */
  if (comedi_dio_bitfield(dev,subdevice,write_mask,&bits) < 0) {
    mexErrMsgTxt("comedi_dio_bitfield: error calling comedi library");
  }

  /* Return data to matlab */
  plhs[0] = mxCreateScalarDouble((double) bits);
} 
