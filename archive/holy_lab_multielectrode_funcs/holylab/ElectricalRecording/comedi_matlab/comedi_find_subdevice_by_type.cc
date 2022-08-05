/*
    comedi_find_subdevice_by_type
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
  int subdevice_type;
  comedi_t *dev;
  int nrows,ncols;
  int ret;
  const mxArray *curarg;

  if (nrhs < 2)
    mexErrMsgTxt("comedi_find_subdevice_by_type: requires 2 or 3 input arguments");

  /* Device */
  curarg = prhs[0];
  buflen = matstr_len(curarg);
  devname = (char*) mxCalloc(buflen,sizeof(char));
  matstr_to_cstr(curarg,devname,buflen);
  dev = get_comedi_device_safe(devname);
  mxFree(devname);
  /* Type */
  curarg = prhs[1];
  subdevice_type = getvalfromstr(curarg,subdevnames,n_subdevnames);
  /* Subdevice start number */
  if (nrhs > 2) {
    curarg = prhs[2];
    if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
      mexErrMsgTxt("comedi_find_subdevice_by_type: subdevice must be a real scalar (nonneg integer)");
    subdevice = (unsigned int) mxGetScalar(curarg);
  }
  else
    subdevice = 0;


  /* Call the comedi library */
  ret = comedi_find_subdevice_by_type(dev,subdevvals[subdevice_type],subdevice);
  if (ret < 0)
    mexErrMsgTxt("comedi_find_subdevice_by_type: can't find this subdevice type");

  /* Return the subdevice number */
  plhs[0] = mxCreateScalarDouble((double) ret);
} 
