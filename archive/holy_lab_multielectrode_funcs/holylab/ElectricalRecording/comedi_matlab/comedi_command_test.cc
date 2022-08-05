/*
    comedi_command_test
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

const int dims[] = {1,1};

char *cmdstrings[]={
        "success",
        "invalid source",
        "source conflict",
        "invalid argument",
        "argument conflict",
        "invalid chanlist",
};


void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  char *devname;
  int buflen;
  const mxArray *curarg;
  comedi_t *dev;
  comedi_cmd comcmd, *cmd = &comcmd;
  int ret;

  if (nrhs < 2)
    mexErrMsgTxt("comedi_command_test: requires device name and command");

  /* Device */
  curarg = prhs[0];
  buflen = matstr_len(curarg);
  devname = (char*) mxCalloc(buflen,sizeof(char));
  matstr_to_cstr(prhs[0],devname,buflen);
  dev = get_comedi_device_safe(devname);
  mxFree(devname);
  /* Command */
  curarg = prhs[1];
  matstruct_to_cmd(curarg,cmd);

  /* Now do the test */
  /* This will cause some fields to be edited */
  ret = comedi_command_test(dev,cmd);

  /* Set up return values */
  plhs[0] = mxCreateString(cmdstrings[ret]);
  if (nlhs > 1) {
    plhs[1] = mxCreateStructArray(2,dims,n_cmdfields-3,cmdfields);
    cmd_to_matstruct(cmd,plhs[1]);
  }

  comedi_cmd_mxFree(cmd);
}
