/*
    comedi_cancel_and_flush
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
#include <unistd.h>
#include <comedilib.h>
#include <comedipersist.h>
#include "commat.h"


const int sizenull = 1024;
sampl_t devnull[sizenull];

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int buflen;
  int curlength,maxlength;  /* all in units of samples, not bytes */
  int nread;                /* in units of bytes */
  int i;
  char *devname;
  comedi_t *dev;
  int subdevice;
  const mxArray *curarg;
  int ret;

  if (nrhs < 2)
    mexErrMsgTxt("syntax: comedi_cancel_and_flush(devicename,subdevice)");

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
    mexErrMsgTxt("comedi_cancel_and_flush: subdevice must be a real scalar (nonneg integer)");
  subdevice = (unsigned int) mxGetScalar(curarg);


  /* Terminate the command */
  ret = comedi_cancel(dev,subdevice);
  if (ret < 0)
    mexErrMsgTxt("comedi_cancel_and_flush: error cancelling command");

  /* If it's of type AI, clear read buffer (is this necessary?) */
  if (comedi_get_subdevice_type(dev,subdevice) == COMEDI_SUBD_AI) {
    nread = 1;
    while (nread > 0) {
      nread = read(comedi_fileno(dev),devnull,sizenull*sizeof(sampl_t));
      /* Should do some error checking here? */
    }
  }

  return;
}
