/*
    comedi_do_insnlist
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
  const mxArray *curarg;
  comedi_t *dev;
  comedi_insnlist comilist,*ilist = &comilist;
  int i, ret;
  comedi_insn in2[2];
  lsampl_t data[1];
  double *dP;

  if (nrhs < 2)
    mexErrMsgTxt("comedi_do_insnlist: requires device name and instructions");

  /* Device */
  curarg = prhs[0];
  buflen = matstr_len(curarg);
  devname = (char*) mxCalloc(buflen,sizeof(char));
  matstr_to_cstr(prhs[0],devname,buflen);
  dev = get_comedi_device_safe(devname);
  mxFree(devname);
  /* Instruction list */
  curarg = prhs[1];
  matstruct_to_insn(curarg,ilist);

  /* Execute the instruction */
  /*
   * This is the version you'd think of first, but there seems to be a
   * comedi bug (as of 0.7.62/lib0.7.18) which prevents this from working
  ret = comedi_do_insnlist(dev,ilist);
  plhs[0] = mxCreateScalarDouble(ret);
   */
  
  // Here's a (near) substitute which works. The funny thing is,
  // comedi_do_insn calls comedi_do_insnlist! It must be in the driver code.
  plhs[0] = mxCreateDoubleMatrix(1,ilist->n_insns,mxREAL);
  dP = mxGetPr(plhs[0]);
  for (i = 0; i < ilist->n_insns; i++)
    //mexPrintf("Executing instruction %d with result %d\n",i,
    dP[i] = (double) comedi_do_insn(dev,&(ilist->insns[i]));
  //);

  mxFree(ilist->insns);
}
