/*
    comedi_open
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
  char *devname;
  int buflen;

  if (nrhs < 1)
    mexErrMsgTxt("comedi_open: requires device name");

  buflen = matstr_len(prhs[0]);
  devname = (char*) mxCalloc(buflen,sizeof(char));
  matstr_to_cstr(prhs[0],devname,buflen);
  if (comedi_device_isopen(devname))
    mexErrMsgTxt("Device already open");
  open_comedi_device(devname);
  mxFree(devname);
}
