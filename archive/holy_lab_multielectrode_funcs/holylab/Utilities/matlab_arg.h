//function set used to ease 

#ifndef MATLAB_ARG_H
#define MATLAB_ARG_H


#include <string>
using std::string;


#include "mex.h"

namespace matlab_arg {

template<typename T>struct TMatrixParam {
   T * p;  //the buf holds items
   int m;  //# of rows
   int n;  //# of cols

};

typedef TMatrixParam<double> TDoubleMatrixParam;


//usage sample:
//  getDoubleArg(prhs[0]);
double getDoubleArg(const mxArray * arg);

//usage sample:
//plhs[0]=setDoubleArg(10);
//or:  setDoubleArg(plhs[0], 10);
mxArray * setDoubleArg(double value);
mxArray * & setDoubleArg(mxArray * & rArg, double value);


string getStringArg(const mxArray * arg);
mxArray * setStringArg(string& value);
mxArray * & setStringArg(mxArray * & rArg, string& value);


TDoubleMatrixParam getDoubleMatrixArg(const mxArray * arg);

} //namespace, matlab_arg


#endif //MATLAB_ARG_H
