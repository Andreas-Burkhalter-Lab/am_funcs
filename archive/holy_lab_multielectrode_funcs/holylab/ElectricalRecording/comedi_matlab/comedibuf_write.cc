/*
    comedibuf_write
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

/* 
   Syntax: ret = comedibuf_write(device,buffer)
   Data are written to the output buffer until the buffer is full.

   On exit, the value of ret is the number of _samples_ written.
   Beware of the difference from C-code, where the number of _bytes_
   is returned.  (Each sample is 2 bytes long.)
  */


typedef short int16;

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int buflen;
  int datalength;  /* all in units of samples, not bytes */
  int nwritten, maxread;       /* in units of bytes */
  int i;
  char *devname;
  unsigned int channel;
  unsigned int range;
  char *arefstr;
  int indx;
  unsigned int aref;
  lsampl_t data;
  comedi_t *dev;
  const mxArray *curarg;
  int16 *buf;

  if (nrhs < 2)
    mexErrMsgTxt("syntax: comedibuf_write(devicename,buffer)");

  /* Device */
  curarg = prhs[0];
  buflen = matstr_len(curarg);
  devname = (char*) mxCalloc(buflen,sizeof(char));
  matstr_to_cstr(curarg,devname,buflen);
  dev = get_comedi_device_safe(devname);
  mxFree(devname);
  /* Buffer */
  curarg = prhs[1];
  if (!mxIsInt16(curarg))
    mexErrMsgTxt("comedibuf_write: buffer must be an int16");
  buf = (int16 *) mxGetData(curarg);
  datalength = mxGetNumberOfElements(curarg);

  /* Write the data */
  nwritten = write(comedi_fileno(dev),buf,datalength*sizeof(int16));
  if (nwritten < 0) {
    perror("write");
    mexErrMsgTxt("Error writing data to output channel(s)");
  }
  /*
  if (nwritten < datalength) {
    mexPrintf("Number of samples written: %d; number requested: %d\n",nwritten,datalength);
    mexErrMsgTxt("Less data was written than requested");
  }
  */

  plhs[0] = mxCreateScalarDouble(nwritten/sizeof(int16));
  return;
}
