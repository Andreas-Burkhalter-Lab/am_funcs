#include "mat.h"
#include "mat_variables.h"

/* 
 * mat_load_variables: load particular variables from a .mat file
 * mat_save_variables: save variables to a .mat file
 *
 * These are useful functions particularly for testing, debugging, and
 * profiling MEX files.  With them you can build stand-alone C
 * programs that use Matlab for all the file I/O: just save a .mat
 * file with all the input variables (named appropriately) and call
 * your mexFunction.
 *
 */
 
void mat_load_variables(const char *filename,const char *varnames[],int n_variables,mxArray *vars[])
{
  MATFile *matfp;
  int varIterator;

  printf("About to open file %s:",filename);
  matfp = matOpen(filename,"r");
  if (matfp == NULL)
    error(1,ENOENT,"Can't open file %s for reading",filename);
  else
    printf("done\n");

  for (varIterator = 0; varIterator < n_variables; varIterator++) {
    printf("Working on variable %s:\n",varnames[varIterator]);
    vars[varIterator] = matGetVariable(matfp,varnames[varIterator]);
    if (vars[varIterator] == NULL)
      error(1,0,"Variable %s is not defined",varnames[varIterator]);
    else
      printf(" done\n");
  }

  if (matClose(matfp) != 0)
    error(1,0,"Error closing file %s",filename);
}

void mat_save_variables(const char *filename,const char *varnames[],int n_variables,mxArray *vars[])
{
  MATFile *matfp;
  int varIterator, status;

  matfp = matOpen(filename,"w");
  if (matfp == NULL)
    error(1,ENOENT,"Can't open file %s for writing",filename);

  for (varIterator = 0; varIterator < n_variables; varIterator++) {
    status = matPutVariable(matfp,varnames[varIterator],vars[varIterator]);
    if (status != 0)
      error(1,0,"Error saving variable %s",varnames[varIterator]);
  }

  if (matClose(matfp) != 0)
    error(1,0,"Error closing file %s",filename);
}

void mexErrMsgTxt(const char *msg)
{
  error(1,0,"%s\n",msg);
}

void mexWarnMsgTxt(const char *msg)
{
  printf("Warning:  %s\n",msg);
}

/*
int mexPrintf(const char *msg, ...)
{
  return printf(msg,...);
}

void mexErrMsgIdAndTxt(const char *errorid,const char *errmsg,...)
{
  error(1,0,errmsg,...);
}
*/
