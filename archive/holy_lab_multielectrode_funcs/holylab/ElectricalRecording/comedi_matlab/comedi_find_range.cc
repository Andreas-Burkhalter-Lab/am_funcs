/*
    comedi_find_range
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
  char *unitstr;
  int indx;
  unsigned int unit;
  comedi_t *dev;
  double minv,maxv,*rangep;
  int nrows,ncols;
  int ret;
  const mxArray *curarg;

  if (nrhs < 5)
    mexErrMsgTxt("comedi_find_range: requires 5 input arguments");

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
    mexErrMsgTxt("comedi_find_range: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);
  /* Channel */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_find_range: channel must be a real scalar (nonneg integer)");
  channel = (unsigned int) mxGetScalar(curarg);
  /* Unit */
  curarg = prhs[3];
  unit = getvalfromstr(curarg,unitnames,n_unit);
  /* Range */
  curarg = prhs[4];
  nrows = mxGetM(curarg);
  ncols = mxGetN(curarg);
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || nrows*ncols != 2)
    mexErrMsgTxt("comedi_find_range: channel must be a real 2-vector [min max]");
  rangep = mxGetPr(curarg);
  minv = rangep[0];
  maxv = rangep[1];


  /* Call the comedi library */
  ret = comedi_find_range(dev,subdevice,channel,unit,minv,maxv);
  if (ret < 0)
    mexErrMsgTxt("comedi_find_range: no matching range available");

  /* Return data to matlab */
  plhs[0] = mxCreateScalarDouble(ret);
} 
