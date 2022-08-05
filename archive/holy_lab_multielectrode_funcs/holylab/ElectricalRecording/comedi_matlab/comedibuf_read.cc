/*
    comedibuf_read
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
   Syntax: ret = comedibuf_read(device,buffer,curlength)
   Data are dumped into buffer, starting at offset curlength
   and continuing until the buffer is full or the read command
   is finished.  If the read command finishes before the buffer
   is full, the buffer is resized to make the buffer full.

   On exit, the value of ret is:
     -1: read command is finished with no more data to return
     >=0: the new curlength
   You can tell when the buffer is ready by checking whether
   curlength == prod(size(buffer)).  Note that if the read
   command finishes before the buffer is full, this function will
   resize the buffer so that this equality will be valid; therefore,
   be sure to compute prod(size(buffer)) after this command returns.
   (Also, be sure to create a new buffer of the desired size if you
   do any subsequent reads!)
*/


typedef short int16;

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  int buflen;
  int curlength,maxlength;  /* all in units of samples, not bytes */
  int nread, maxread;       /* in units of bytes */
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

  if (nrhs < 3)
    mexErrMsgTxt("syntax: comedibuf_read(devicename,buffer,curlength)");

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
    mexErrMsgTxt("comedibuf_read: buffer must be an int16");
  buf = (int16 *) mxGetData(curarg);
  maxlength = mxGetNumberOfElements(curarg);
  /* Curlength */
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !IsScalar(curarg))
    mexErrMsgTxt("comedibuf_read: curlength must be a real scalar (nonneg integer)");
  curlength = (unsigned int) mxGetScalar(curarg);
  if (curlength >= maxlength)
    mexErrMsgTxt("comedibuf_read: curlength is >= the total buffer size.\nDid you forget to reset curlength = 0 after a complete buffer read?");



  /* Get the data */
  maxread = (maxlength-curlength)*sizeof(sampl_t);
  nread = read(comedi_fileno(dev),buf+curlength,maxread);
  //mexPrintf("fileno %d, base %p, readbufoffset %d, readstart %d, maxlength %d, maxread %d, nread %d\n",comedi_fileno(dev),buf,curlength*sizeof(sampl_t),buf+curlength,maxlength,maxread,nread);
  if (nread < 0) {
    perror("read");
    mexErrMsgTxt("comedibuf_read: Device reading error. Buffer overflowed before\nread? Try increasing the input buffer size with\ncomedi_set_buffer_size.");
  }
  if (nread % sizeof(sampl_t))
    mexErrMsgTxt("The number of bytes read was not divisible by sizeof(sampl_t)!");
  curlength += nread/sizeof(sampl_t);
  /* 
     Are we done with this buffer?
     We could be done either because it's full, or because the
     read command is finished.
  */
  if (curlength == maxlength || nread == 0) {
    if (curlength) {
      /* Normal exit */
      if (curlength < maxlength) {
	/* read command is finished */
	mxSetN(plhs[0],(maxlength-curlength)/mxGetM(plhs[0]));
      }
      plhs[0] = mxCreateScalarDouble(curlength);   /* read probably isn't finished */
    }
    else
      plhs[0] = mxCreateScalarDouble(-1); /* read command finished, but we're
					     still trying to get more! */
  }
  else
    plhs[0] = mxCreateScalarDouble(curlength);   /* buffer not done, read prob. not done */
  return;
}
