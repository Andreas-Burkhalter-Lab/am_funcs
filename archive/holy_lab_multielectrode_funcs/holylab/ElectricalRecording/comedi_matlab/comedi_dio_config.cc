/*
    comedi_dio_config
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
  unsigned int dio_direction;
  int subdevtype;
  int indx;
  comedi_t *dev;
  const mxArray *curarg;

  if (nrhs < 4)
    mexErrMsgTxt("comedi_dio_config: requires 4 input arguments");

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
    mexErrMsgTxt("comedi_dio_config: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);
  /* Channel */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedi_dio_config: channel must be a real scalar (nonneg integer)");
  channel = (unsigned int) mxGetScalar(curarg);
  /* Direction */
  curarg = prhs[3];
  dio_direction = getvalfromstr(curarg,dio_direction_names,n_dio_directions);


  /* Check to be sure this subdevice is DIO */
  subdevtype = comedi_get_subdevice_type(dev,subdevice);
  if (subdevtype != COMEDI_SUBD_DIO)
    mexErrMsgTxt("Subdevice is not of type DIO");

  /* Configure the channel */
  if (comedi_dio_config(dev,subdevice,channel,dio_direction) < 0)
    mexErrMsgTxt("comedi_dio_config: error configuring digital device");
} 
