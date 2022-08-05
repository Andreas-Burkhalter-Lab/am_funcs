#include <error.h>
#include <errno.h>

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

#ifdef __cplusplus
extern "C" {
#endif
void mat_load_variables(const char *filename,const char *varnames[],int n_variables,mxArray *vars[]);
void mat_save_variables(const char *filename,const char *varnames[],int n_variables,mxArray *vars[]);
void mexErrMsgTxt(const char *msg);
void mexWarnMsgTxt(const char *msg);
/*
int mexPrintf(const char *msg, ...);
void mexErrMsgIdAndTxt(const char *errorid,const char *errmsg,...);
*/
#ifdef __cplusplus
}
#endif
#define mexPrintf printf
/*
#define mexErrMsgIdAndTxt(a, b, rest...) \
           error(1, 0, b , ## rest)
*/
#define mexErrMsgIdAndTxt(a, b, ...) error (1, 0, b, ##__VA_ARGS__)
#define mexWarnMsgIdAndTxt(a, b, ...) printf (b, ##__VA_ARGS__)
#define mexFunctionName() "<function name>"
