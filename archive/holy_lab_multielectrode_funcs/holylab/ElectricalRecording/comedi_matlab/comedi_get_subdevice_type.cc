/*
    comedi_get_subdevice_type
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
  comedi_t *dev;
  int nrows,ncols;
  int ret;
  const mxArray *curarg;

  if (nrhs < 2)
    mexErrMsgTxt("comedi_get_subdevice_type: requires 2 input arguments");

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
    mexErrMsgTxt("comedi_get_subdevice_type: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);


  /* Call the comedi library */
  ret = comedi_get_subdevice_type(dev,subdevice);
  if (ret < 0)
    mexErrMsgTxt("comedi_get_subdevice_type: error fetching subdevice type");

  /* Convert to a string and output*/
  ret = intmatch(ret,subdevvals,n_subdevnames);
  if (ret < 0)
    mexErrMsgTxt("comedi_get_subdevice_type: error matching subdevice type");
  plhs[0] = mxCreateString(subdevnames[ret]);
} 
