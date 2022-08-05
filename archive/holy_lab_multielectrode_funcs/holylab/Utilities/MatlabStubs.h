#ifndef MATLABSTUBS_H
#define MATLABSTUBS_H

#include <error.h>
#include <errno.h>

/* 
 * These are stubs for the common mex function calls, resolving them
 * to C-library functions.  These assist in compiling a standalone
 * version of a MEX file, which can be a very helpful step in
 * debugging MEX files (e.g., allow setting breakpoints, valgrind
 * memcheck, etc.).
 *
 */

#define mexPrintf printf
#define mexErrMsgTxt(a) error(1, 0, a)
#define mexErrMsgIdAndTxt(a, b, ...) error (1, 0, b, ##__VA_ARGS__)
#define mexWarnMsgIdAndTxt(a, b, ...) printf (b, ##__VA_ARGS__)
#define mexFunctionName() "<function name>"
#define mexEvalString printf

#endif // MATLABSTUBS_H
