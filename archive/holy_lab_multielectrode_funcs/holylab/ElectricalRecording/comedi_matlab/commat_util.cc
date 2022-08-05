/*
    commat_util
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
#include <string.h>
#include "commat.h"

int warnMsg = 1;

const char *cmdfields[] = {
  "subdev",
  "flags",
  "start_src",
  "start_arg",
  "scan_begin_src",
  "scan_begin_arg",
  "convert_src",
  "convert_arg",
  "scan_end_src",
  "scan_end_arg",
  "stop_src",
  "stop_arg",
  "packedchanlist",
  "chanlist",
  "rangelist",
  "areflist",
};
const int n_cmdfields = 16;

const char *arefnames[] = {
  "AREF_GROUND",
  "AREF_COMMON",
  "AREF_DIFF",
  "AREF_OTHER",
};
const unsigned int arefvals[] = {
  AREF_GROUND,
  AREF_COMMON,
  AREF_DIFF,
  AREF_OTHER,
};
const int n_aref = 4;

const char *unitnames[] = {
  "UNIT_volt",
  "UNIT_mA",
  "UNIT_none",
};
const unsigned int unitvals[] = {
  UNIT_volt,
  UNIT_mA,
  UNIT_none,
};
const int n_unit = 3;

const char *srcnames[] = {
  "TRIG_NONE",
  "TRIG_NOW",
  "TRIG_FOLLOW",
  "TRIG_TIME",
  "TRIG_TIMER",
  "TRIG_COUNT",
  "TRIG_EXT",
  "TRIG_INT",
  "TRIG_OTHER",
};
const unsigned int srcvals[] = {
  TRIG_NONE,
  TRIG_NOW,
  TRIG_FOLLOW,
  TRIG_TIME,
  TRIG_TIMER,
  TRIG_COUNT,
  TRIG_EXT,
  TRIG_INT,
  TRIG_OTHER,
};
const int n_src = 9;

const char *subdevnames[] = {
  "COMEDI_SUBD_UNUSED",
  "COMEDI_SUBD_AI",
  "COMEDI_SUBD_AO",
  "COMEDI_SUBD_DI",
  "COMEDI_SUBD_DO",
  "COMEDI_SUBD_DIO",
  "COMEDI_SUBD_COUNTER",
  "COMEDI_SUBD_TIMER",
  "COMEDI_SUBD_MEMORY",
  "COMEDI_SUBD_CALIB",
  "COMEDI_SUBD_PROC"};
const unsigned int subdevvals[] = {
  COMEDI_SUBD_UNUSED,
  COMEDI_SUBD_AI,
  COMEDI_SUBD_AO,
  COMEDI_SUBD_DI,
  COMEDI_SUBD_DO,
  COMEDI_SUBD_DIO,
  COMEDI_SUBD_COUNTER,
  COMEDI_SUBD_TIMER,
  COMEDI_SUBD_MEMORY,
  COMEDI_SUBD_CALIB,
  COMEDI_SUBD_PROC};
const int n_subdevnames = 10;

const char *dio_direction_names[] = {
  "COMEDI_INPUT",
  "COMEDI_OUTPUT"};
const unsigned int dio_direction_vals[] = {
  COMEDI_INPUT,
  COMEDI_OUTPUT};
const int n_dio_directions = 2;

const char *insn_fields[] = {
  "insn",
  "n",
  "data",
  "subdev",
  "chan",
  "range",
  "aref",
};
const int n_insn_fields = 8;

const char *insn_names[] = {
  "INSN_READ",
  "INSN_WRITE",
  "INSN_BITS",
  "INSN_CONFIG",
  "INSN_GTOD",
  "INSN_WAIT",
  "INSN_INTTRIG",
};
const unsigned int insn_vals[] = {
  INSN_READ,
  INSN_WRITE,
  INSN_BITS,
  INSN_CONFIG,
  INSN_GTOD,
  INSN_WAIT,
  INSN_INTTRIG};
const int n_insn = 7;
  

/*
void comedipersist_error(const char *msg)
{
  mexErrMsgTxt(msg);
}
*/

/*
FILE* stderrstream = NULL;

void capture_stderr()
{
  if (stderrstream != NULL)
    mexErrMsgTxt("capture_stderr: stderr is already being captured");
  stderrstream = freopen("stderr.out","w, type=memory", stderr);
  if (stderrstream == NULL) {
    perror("capture_stderr");
    mexErrMsgTxt("capture_stderr: error capturing stderr");
  }
}

char* stderrcpy(char *dest)
{
  int filelength;
  char *cpyresult;

  if (stderrstream == NULL)
    mexErrMsgTxt("stderrcpy: stderr is not being captured");
  if (fclose(stderrstream))
    mexErrMsgTxt("stderrcpy: error closing capture stream");
  stderrstream = fopen("stderr.out","r");
  if (stderrstream == NULL)
    mexErrMsgTxt("stderrcpy: error reopening capture file in read mode");
  if (fseek(stderrstream,0,SEEK_END))
    mexErrMsgTxt("stderrcpy: error seeking end of capture stream");
  filelength = ftell(stderrstream);
  rewind(stderrstream);
  dest = (char *) mxCalloc(filelength, sizeof(char));
  if (!dest)
    mexErrMsgTxt("stderrcpy: error allocating string");
  cpyresult = fgets(dest,filelength,stderrstream);
  if (cpyresult == NULL)
    mexErrMsgTxt("stderrcpy: error copying stream");
  if (fclose(stderrstream))
    mexErrMsgTxt("stderrcpy: error closing capture stream after read");
  stderrstream = NULL;
  return dest;
}
*/

int matstr_len(const mxArray *matstr)
{
  int nrows,ncols;

  if (!mxIsChar(matstr))
    mexErrMsgTxt("matstr_len: argument should be a string");
  nrows = mxGetM(matstr);
  ncols = mxGetN(matstr);
  if (nrows != 1 && ncols != 1)
    mexErrMsgTxt("matstr_len: argument must be a string vector, not matrix");

  return nrows*ncols*sizeof(mxChar) + 1;
}

char *matstr_to_cstr(const mxArray *matstr,char *cstr,int buflen)
{
  if (cstr == NULL)
    mexErrMsgTxt("matstr_to_cstr: error allocating space for c string");
  if (mxGetString(matstr,cstr,buflen))
    mexErrMsgTxt("matstr_to_cstr: some error copying device name");
  return cstr;
}

bool IsScalar(const mxArray *m)
{
  return (mxGetNumberOfElements(m) == 1);
}

int intmatch(unsigned int uint,const unsigned int vals[],int n_vallist)
{
  int i;

  i = 0;
  while (i < n_vallist)
    if (uint == vals[i])
      break;
    else
      i++;
  if (i < n_vallist)
    return i;
  else
    return -1;
}
      
int strmatch(const char *str,const char *strlist[],int n_strlist)
{
  int i;

  i = 0;
  while (i < n_strlist)
    if (!strcmp(str,strlist[i]))
      break;
    else
      i++;
  if (i < n_strlist)
    return i;
  else {
    if (warnMsg) {
      mexPrintf("%s is not among the following choices:\n",str);
      for (i = 0; i < n_strlist; i++)
	mexPrintf("  %s\n",strlist[i]);
    }
    return -1;
  }
}


int getvalfromstr(const mxArray *matstr,const char *strs[],int nvals)
{
  int buflen;
  char *cstr;
  int ret;

  buflen = matstr_len(matstr);
  cstr = (char*) mxCalloc(buflen,sizeof(char));
  if (!cstr)
    mexErrMsgTxt("getvalfromstr: error allocating space for string");
  matstr_to_cstr(matstr,cstr,buflen);
  ret = strmatch(cstr,strs,nvals);
  mxFree(cstr);
  if (ret < 0)
    mexErrMsgTxt("getvalfromstr: string is not in list");
  return ret;
}

comedi_t* get_comedi_device_safe(const char *devname)
{
  if (!comedi_device_isopen(devname)) {
    mexPrintf("Error using device %s:\n",devname);
    mexErrMsgTxt("Device is not open");
  }
  return get_comedi_device(devname);
}
			    

void matstruct_to_cmd(const mxArray *matcmd,comedi_cmd *comcmd)
{
  int nfields,i;
  const mxArray *curarg;
  const char *curfieldname;
  int curfieldindx;
  commat_chanspec cspec;
  double *chanp;
  int j;

  if (!mxIsStruct(matcmd))
    mexErrMsgTxt("Error matstruct_to_cmd: first input must be a MATLAB structure");
  if (mxGetNumberOfElements(matcmd) != 1)
    mexErrMsgTxt("Error matstruct_to_cmd: first input must be have only 1 structure element");

  memset(comcmd,0,sizeof(*comcmd));
  prep_chanspec(&cspec);

  nfields = mxGetNumberOfFields(matcmd);

  for (i = 0; i < nfields; i++) {
    curfieldname = mxGetFieldNameByNumber(matcmd,i);
    curfieldindx = strmatch(curfieldname, cmdfields, n_cmdfields);
    switch (curfieldindx) {
    case -1:
      mexPrintf("Fieldname %s not recognized\n",curfieldname);
      mexErrMsgTxt("matstruct_to_cmd");
      break;
    case 0:   // subdev
      comcmd->subdev = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 1:   // flags
      if (!IsScalar(mxGetFieldByNumber(matcmd,0,i)))
	mexErrMsgTxt("nonscalar flags not yet implemented");
      comcmd->flags = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 2:   // start_src
      comcmd->start_src = srcvals[getvalfromstr(mxGetFieldByNumber(matcmd,0,i),srcnames,n_src)];
      break;
    case 3:   // start_arg
      comcmd->start_arg = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 4:   // scan_begin_src
      comcmd->scan_begin_src = srcvals[getvalfromstr(mxGetFieldByNumber(matcmd,0,i),srcnames,n_src)];
      break;
    case 5:   // scan_begin_arg
      comcmd->scan_begin_arg = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 6:   // convert_src
      comcmd->convert_src = srcvals[getvalfromstr(mxGetFieldByNumber(matcmd,0,i),srcnames,n_src)];
      break;
    case 7:   // convert_arg
      comcmd->convert_arg = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 8:   // scan_end_src
      comcmd->scan_end_src = srcvals[getvalfromstr(mxGetFieldByNumber(matcmd,0,i),srcnames,n_src)];
      break;
    case 9:   // scan_end_arg
      comcmd->scan_end_arg = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 10:  // scan_src
      comcmd->stop_src = srcvals[getvalfromstr(mxGetFieldByNumber(matcmd,0,i),srcnames,n_src)];
      break;
    case 11:  // stop_arg
      comcmd->stop_arg = (unsigned int) mxGetScalar(mxGetFieldByNumber(matcmd,0,i));
      break;
    case 12:  // packedchanlist
      curarg = mxGetFieldByNumber(matcmd,0,i);
      cspec.nchans_packed = mxGetNumberOfElements(curarg);
      cspec.packedchanlist = (unsigned int *) mxCalloc(cspec.nchans_packed,sizeof(unsigned int));
      if (!cspec.packedchanlist)
	mexErrMsgTxt("matstruct_to_cmd: error allocating space for packedchanlist");
      chanp = mxGetPr(curarg);
      for (j = 0; j < cspec.nchans_packed; j++)
	cspec.packedchanlist[j] = (unsigned int) chanp[j];
      break;
    case 13:  // chanlist
      curarg = mxGetFieldByNumber(matcmd,0,i);
      cspec.nchans_chan = mxGetNumberOfElements(curarg);
      cspec.chanlist = (unsigned int *) mxCalloc(cspec.nchans_chan,sizeof(unsigned int));
      if (!cspec.chanlist)
	mexErrMsgTxt("matstruct_to_cmd: error allocating space for chanlist");
      chanp = mxGetPr(curarg);
      for (j = 0; j < cspec.nchans_chan; j++)
	cspec.chanlist[j] = (unsigned int) chanp[j];
      break;
    case 14:  // rangelist
      curarg = mxGetFieldByNumber(matcmd,0,i);
      cspec.nchans_range = mxGetNumberOfElements(curarg);
      cspec.rangelist = (unsigned int *) mxCalloc(cspec.nchans_range,sizeof(unsigned int));
      if (!cspec.rangelist)
	mexErrMsgTxt("matstruct_to_cmd: error allocating space for rangelist");
      chanp = mxGetPr(curarg);
      for (j = 0; j < cspec.nchans_range; j++)
	cspec.rangelist[j] = (unsigned int) chanp[j];
      break;
    case 15:  // areflist
      curarg = mxGetFieldByNumber(matcmd,0,i);
      cspec.nchans_aref = mxGetNumberOfElements(curarg);
      cspec.areflist = (unsigned int *) mxCalloc(cspec.nchans_aref,sizeof(unsigned int));
      if (!cspec.areflist)
	mexErrMsgTxt("matstruct_to_cmd: error allocating space for areflist");
      chanp = mxGetPr(curarg);
      for (j = 0; j < cspec.nchans_aref; j++)
	cspec.areflist[j] = arefvals[(int) chanp[j]];
      break;
    }
  }
  comcmd->chanlist_len = do_chan_packing(&(comcmd->chanlist),&cspec);
}

void cmd_to_matstruct(const comedi_cmd *comcmd, mxArray *matcmd)
{
  mxArray *packedchanlist;
  double *pchan;
  int j;

  mxSetField(matcmd,0,"subdev",mxCreateScalarDouble(comcmd->subdev));
  mxSetField(matcmd,0,"flags",mxCreateScalarDouble(comcmd->flags));
  mxSetField(matcmd,0,"start_src",mxCreateString(srcnames[intmatch(comcmd->start_src,srcvals,n_src)]));
  mxSetField(matcmd,0,"start_arg",mxCreateScalarDouble(comcmd->start_arg));
  mxSetField(matcmd,0,"scan_begin_src",mxCreateString(srcnames[intmatch(comcmd->scan_begin_src,srcvals,n_src)]));
  mxSetField(matcmd,0,"scan_begin_arg",mxCreateScalarDouble(comcmd->scan_begin_arg));
  mxSetField(matcmd,0,"convert_src",mxCreateString(srcnames[intmatch(comcmd->convert_src,srcvals,n_src)]));
  mxSetField(matcmd,0,"convert_arg",mxCreateScalarDouble(comcmd->convert_arg));
  mxSetField(matcmd,0,"scan_end_src",mxCreateString(srcnames[intmatch(comcmd->scan_end_src,srcvals,n_src)]));
  mxSetField(matcmd,0,"scan_end_arg",mxCreateScalarDouble(comcmd->scan_end_arg));
  mxSetField(matcmd,0,"stop_src",mxCreateString(srcnames[intmatch(comcmd->stop_src,srcvals,n_src)]));
  mxSetField(matcmd,0,"stop_arg",mxCreateScalarDouble(comcmd->stop_arg));
  packedchanlist = mxCreateDoubleMatrix((int) comcmd->chanlist_len,1,mxREAL);
  pchan = mxGetPr(packedchanlist);
  for (j = 0; j < (int) comcmd->chanlist_len; j++) {
    pchan[j] = comcmd->chanlist[j];
  }
  mxSetField(matcmd,0,"packedchanlist",packedchanlist);
}
	
void matstruct_to_insn(const mxArray *matinsn,comedi_insnlist *cominsnlist)
{
  int nfields,ninsns,i;
  const char *curfieldname;
  int curfieldindx;
  unsigned int chan,range,aref;
  int k,dataindex,temp;
  comedi_insn *curinsn;

  if (!mxIsStruct(matinsn))
    mexErrMsgTxt("Error matstruct_to_insn: first input must be a MATLAB structure");

  // Determine the number of instructions & allocate space
  ninsns = mxGetNumberOfElements(matinsn);
  if (ninsns < 1)
    mexErrMsgTxt("Error matstruct_to_insn: first input must be have at least 1 structure element");
  cominsnlist->n_insns = ninsns;
  cominsnlist->insns = (comedi_insn *) mxCalloc(ninsns,sizeof(comedi_insn));
  memset(cominsnlist->insns,0,ninsns*sizeof(comedi_insn));

  nfields = mxGetNumberOfFields(matinsn);

  // Loop over instructions
  for (k = 0; k < ninsns; k++) {
    dataindex = -1;    // indicate no data field processed yet
    chan = range = aref = 0;
    curinsn = &(cominsnlist->insns[k]);
    // Loop over fieldnames
    for (i = 0; i < nfields; i++) {
      curfieldname = mxGetFieldNameByNumber(matinsn,i);
      curfieldindx = strmatch(curfieldname, insn_fields, n_insn_fields);
      switch (curfieldindx) {
      case -1:
	mexPrintf("Fieldname %s not recognized\n",curfieldname);
	mexErrMsgTxt("matstruct_to_insn");
	break;
      case 0:   // insn
	temp = getvalfromstr(mxGetFieldByNumber(matinsn,k,i),insn_names,n_insn);
	curinsn->insn = insn_vals[getvalfromstr(mxGetFieldByNumber(matinsn,k,i),insn_names,n_insn)];
	break;
      case 1:   // n
	curinsn->n = (unsigned int) mygetscalar(mxGetFieldByNumber(matinsn,k,i));
	break;
      case 2:   // data
	if (mxGetFieldByNumber(matinsn,k,i))  // if it's allocated...
	  dataindex = i; // remember index, deal with this later
	break;
      case 3:   // subdev
	curinsn->subdev = (unsigned int) mygetscalar(mxGetFieldByNumber(matinsn,k,i));
	break;
      case 4:  // chan
	chan = (unsigned int) mygetscalar(mxGetFieldByNumber(matinsn,k,i));
	break;
      case 5:  // range
	range = (unsigned int) mygetscalar(mxGetFieldByNumber(matinsn,k,i));
	break;
      case 6:  // aref
	aref = (unsigned int) mygetscalar(mxGetFieldByNumber(matinsn,k,i));
	break;
      }
    }
    // Pack the channel list
    curinsn->chanspec = CR_PACK(chan,range,aref);
    // Now we know enough about the process to deal with the data field
    if (dataindex >= 0) {
      switch (curinsn->insn) {
      case INSN_WAIT:
	curinsn->data = (lsampl_t *) mxCalloc(1,sizeof(lsampl_t));
	*(curinsn->data) = (lsampl_t) mxGetScalar(mxGetFieldByNumber(matinsn,k,dataindex));
	break;
      default:
	mexErrMsgTxt("matstruct_to_insn: data field not implemented for all types yet");
      }
    }
    // Enforce only supported types
    if (curinsn->insn != INSN_WAIT &&
	curinsn->insn != INSN_INTTRIG)
      mexErrMsgTxt("matstruct_to_insn: only WAIT and INTTRIG are currently supported");
  }
}

void prep_chanspec(commat_chanspec *cspec)
{
  cspec->chanlist = cspec->rangelist = cspec->areflist = cspec->packedchanlist = 0;
  cspec->nchans_chan = cspec->nchans_range = cspec->nchans_aref = cspec->nchans_packed = 0;
}

int do_chan_packing(unsigned int **pclist,commat_chanspec *cspec)
{
  // Generic channel-packing function.  But then I realized that
  // comedi_insn only allows scalar channels, so comedi_cmd is the
  // only thing that uses this!

  // Check to see if channels are defined; if so, do CR_PACK
  int j;
  if (cspec->nchans_chan && cspec->nchans_packed == 0) {
    if (!cspec->nchans_range) {
      cspec->rangelist = (unsigned int *) mxCalloc(cspec->nchans_chan,sizeof(unsigned int));
      if (!cspec->rangelist)
	mexErrMsgTxt("do_chan_packing: error allocating space for rangelist");
      for (j = 0; j < cspec->nchans_chan; j++)
	cspec->rangelist[j] = 0;
    }
    else if (cspec->nchans_range != cspec->nchans_chan)
      mexErrMsgTxt("do_chan_packing: the number of channels and ranges do not match");
    if (!cspec->nchans_aref) {
      cspec->areflist = (unsigned int *) mxCalloc(cspec->nchans_chan,sizeof(unsigned int));
      if (!cspec->areflist)
	mexErrMsgTxt("do_chan_packing: error allocating space for areflist");
      for (j = 0; j < cspec->nchans_chan; j++)
	cspec->areflist[j] = AREF_GROUND;
    }
    else if (cspec->nchans_aref != cspec->nchans_chan)
      mexErrMsgTxt("do_chan_packing: the number of channels and arefs do not match");
    cspec->packedchanlist = (unsigned int *) mxCalloc(cspec->nchans_chan,sizeof(unsigned int));
    if (!cspec->packedchanlist)
      mexErrMsgTxt("do_chan_packing: error allocating space for packedchanlist");
    for (j = 0; j < cspec->nchans_chan; j++) {
      cspec->packedchanlist[j] = CR_PACK(cspec->chanlist[j],cspec->rangelist[j],cspec->areflist[j]);
    }
    mxFree(cspec->chanlist);
    mxFree(cspec->rangelist);
    mxFree(cspec->areflist);
    *pclist = cspec->packedchanlist;
    return cspec->nchans_chan;
  }
  else {
    if (cspec->nchans_chan)
      mexErrMsgTxt("do_chan_packing: Can't have both chanlist and packedchanlist");
    *pclist = cspec->packedchanlist;
    return cspec->nchans_packed;
  }
}

void comedi_cmd_mxFree(comedi_cmd *comcmd)
{
  mxFree(comcmd->chanlist);
}


double mygetscalar(const mxArray *m)
{
  if (!m)
    return 0;
  if (mxIsEmpty(m))
    return 0;
  else
    return mxGetScalar(m);
}

/* 
 * The following output functions were copied directly from
 * David A. Schleef's common.c file in the comedilib demo directory.
 */

char *cmd_src(int src,char *buf)
{
        buf[0]=0;
 
        if(src&TRIG_NONE)strcat(buf,"none|");
        if(src&TRIG_NOW)strcat(buf,"now|");
        if(src&TRIG_FOLLOW)strcat(buf, "follow|");
        if(src&TRIG_TIME)strcat(buf, "time|");
        if(src&TRIG_TIMER)strcat(buf, "timer|");
        if(src&TRIG_COUNT)strcat(buf, "count|");
        if(src&TRIG_EXT)strcat(buf, "ext|");
        if(src&TRIG_INT)strcat(buf, "int|");
#ifdef TRIG_OTHER
        if(src&TRIG_OTHER)strcat(buf, "other|");
#endif
 
        if(strlen(buf)==0){
                sprintf(buf,"unknown(0x%08x)",src);
        }else{
                buf[strlen(buf)-1]=0;
        }
 
        return buf;
}
 
void dump_cmd(FILE *out,comedi_cmd *cmd)
{
        char buf[100];
 
        fprintf(out,"start:      %-8s %d\n",
                cmd_src(cmd->start_src,buf),
                cmd->start_arg);
 
        fprintf(out,"scan_begin: %-8s %d\n",
                cmd_src(cmd->scan_begin_src,buf),
                cmd->scan_begin_arg);
 
        fprintf(out,"convert:    %-8s %d\n",
                cmd_src(cmd->convert_src,buf),
                cmd->convert_arg);
 
        fprintf(out,"scan_end:   %-8s %d\n",
                cmd_src(cmd->scan_end_src,buf),
                cmd->scan_end_arg);
 
        fprintf(out,"stop:       %-8s %d\n",
                cmd_src(cmd->stop_src,buf),
                cmd->stop_arg);
}
