/*=================================================================

 *=================================================================*/
#include <string>
#include <vector>
using std::string;
using std::vector;

#include "mex.h"
#include "matlab_arg.h"
using namespace matlab_arg;

//syntax doublearray=split_dbl(astring, delimiter)
//e.g. a=split_dbl('1 2 3', ' '); 
//    we got a=[1 2 3]


//another implementation of splitString() that uses strtok().
//but this one has problem that it skips empty items:
//      e.g. splitString2("a,,c",",") output "a" "c" instead of "a" "" "c"
vector<string> splitStringSkipEmptyItem(string const & rStr, string const delimiter)
{
   char* p=new char[rStr.size()+1];
   strcpy(p, rStr.c_str());
   vector<string> result;
   p=strtok(p, delimiter.c_str());
   while(p){
      result.push_back(p);
      p = strtok(NULL, delimiter.c_str());
   }

   return result;
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    /* Check for proper number of input and  output arguments */    
    if (nrhs !=2) {
        mexErrMsgTxt("two input argument (string_to_split, delimiter) required.");
    } 
    if(nlhs >1 ){
        mexErrMsgTxt("0~1 output arguments required.");
    }

    /* input must be a string */
    if( mxIsChar(prhs[0]) != 1)  mexErrMsgTxt("Input must be a string.");

    //handle empty input string:
    //mexPrintf("matrix dimensions: %d x %d \n", mxGetM(prhs[0]), mxGetN(prhs[0]));
    if(mxGetM(prhs[0])*mxGetN(prhs[0])==0){
       plhs[0] =mxCreateDoubleMatrix(1, 0, mxREAL);
       //or: plhs[0] =mxCreateDoubleMatrix(0, 0, mxREAL); ?
       return;
    }

    /* input must be a row vector */
    if(mxGetM(prhs[0])!=1)  mexErrMsgTxt("Input must be a row vector.");

    if( mxIsChar(prhs[1]) != 1)  mexErrMsgTxt("Input must be a string.");

    /* input must be a row vector */
    if(mxGetM(prhs[1])!=1)  mexErrMsgTxt("Input must be a row vector.");

    string toSplit=getStringArg(prhs[0]);
    string delimiter=getStringArg(prhs[1]);
    vector<string> stringList=splitStringSkipEmptyItem(toSplit, delimiter);

    vector<double> dblList;
    for(int i=0; i<stringList.size(); ++i)dblList.push_back(atof(stringList[i].c_str()));
                    //for long integer in stringList, atof() may lose precision

    //we can always assume one LHS arg there---if user does not supply one explicitly, it is variable ans:
    plhs[0] = mxCreateDoubleMatrix(1, dblList.size(), mxREAL);
    double* p=mxGetPr(plhs[0]);
    copy(dblList.begin(), dblList.end(), p); 

}

