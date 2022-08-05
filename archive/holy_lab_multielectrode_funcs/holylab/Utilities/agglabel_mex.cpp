/*=================================================================
  Aggregate entries with the same label
 *=================================================================*/

#include "mex.h"
#include <vector>

using namespace std;


//---------------------------------------------------------------------

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  long mrows,ncols,nelem,i,labelsize;
  double *label,*labelend,*clabelsize;
  int clabel;
  vector<vector<double> > labeltable;  // This is where the labels are stored
  mxArray *tArray;

  // Check for proper number of input and output arguments
  if (nrhs !=1) {
    mexErrMsgTxt("One input argument required.");
  } 
  if(nlhs != 2)
    mexErrMsgTxt("Two output arguments required.");

  // input must be a real vector
  if( !mxIsDouble(prhs[0]) ) mexErrMsgTxt("Input must be double.");
  if( mxIsComplex(prhs[0]) ) mexErrMsgTxt("Input cannot be complex.");
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  //mexPrintf("mrows: %d, ncols: %d\n",mrows,ncols);
  if (mrows > 1 && ncols > 1) mexErrMsgTxt("Input must be a vector.");
  nelem = mrows*ncols;
  label = mxGetPr(prhs[0]);

  // Now start aggregating. Don't separate this out to a work
  // function, because it requires a lot of interaction with Matlab
  // memory allocation to report the results back to Matlab.
  // Using STL here, with growing containers, allows O(N logN) behavior

  for (i = 0; i < nelem; i++) {
    clabel = int(label[i]);   // current label
    if (clabel > 0) {
      // Do we have a new label we've never seen?
      // Check first if we need to grow the labeltable
      if (clabel > labeltable.capacity())
	labeltable.reserve(2*clabel);  // So don't realloc often
      if (clabel > labeltable.size())
	labeltable.resize(clabel);
      
      // Now add this index onto the end of the list
      labeltable[clabel-1].push_back(double(i+1));
    }
  }
  //mexPrintf("Processed %d elements, with the largest label being %d\n",nelem,labeltable.size());

  // Now convert results to MATLAB format
  plhs[0] = mxCreateCellMatrix(1,labeltable.size());  // the label list
  plhs[1] = mxCreateDoubleMatrix(1,labeltable.size(),mxREAL); // the lengths of lists
  clabelsize = mxGetPr(plhs[1]);
  for (i = 0; i < labeltable.size(); i++) {
    labelsize = labeltable[i].size();
    // Store the length of each label list
    clabelsize[i] = double(labelsize);
    // Allocate a matlab vector for each label list
    tArray = mxCreateDoubleMatrix(1,labelsize,mxREAL);
    mxSetCell(plhs[0],i,tArray);
    // Copy the values
    copy(labeltable[i].begin(),labeltable[i].end(),mxGetPr(tArray));
  }
}
