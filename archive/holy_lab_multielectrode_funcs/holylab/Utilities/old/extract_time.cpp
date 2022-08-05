//compile
//for c++ stand-alone executable:
//g++ -I ~/commong -o extract_time extract_time.cpp ~/commong/miscpp.cpp
//for matlab mex file:
// mex -DMEX -I/home/jason/commong extract_time.cpp ~/commong/miscpp.cpp matlab_arg.cpp
//run:
// ./extract_time "30s 10m and 2 min"

#ifdef MEX
#   include "mex.h"
#   include "matlab_arg.h"
#endif

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

#include <ctype.h>

#include "miscpp.h"

//--------------------------------------------------------------------------
//return time in second
int extract_time(string aTimeString, string* pSimplifiedTimeString=0)
{
   vector<string> tokens;
   tokens=splitStringSkipEmptyItem(aTimeString, " ,\n\t");
   for(int i=0; i<tokens.size(); ++i){
      tokens[i]=strTrim(tokens[i]);
   }

   char* p[]={"s", "sec", "sec.", "second", "seconds",
              "m", "min", "min.", "minute", "minutes",
              "h", "hr", "hr.", "hour", "hours", "hrs", "hrs."
   };
   vector<string> units(p, p+sizeof(p)/sizeof(char*));
   //cout<<units[0]<<endl;
   //cout<<units[units.size()-1]<<endl;
   //cout<<units.size()<<endl;

   vector<string> times;
   int curIdx=0;
   while(curIdx<tokens.size()){
      if(is_all_digits(tokens[curIdx])){
         if(curIdx!=tokens.size()-1 && find_first_of(units, tokens[curIdx+1], false)!=-1){
            times.push_back(tokens[curIdx]+" "+tokens[curIdx+1]);
            ++curIdx; //advance one more than other cases, i.e. advance by 2 for next iteration
         }
      }
      else if(tokens[curIdx].size()>0 && isdigit(tokens[curIdx].at(0))) {
         string digits, after_digits;
         separate_digits(tokens[curIdx], digits, after_digits);
         if(find_first_of(units, after_digits, false)!=-1){
            times.push_back(tokens[curIdx]);
         }
      }//else if, start w/ digit
      ++curIdx;
   }

   // now calc total time and extract time string:
   //copy(times.begin(), times.end(), ostream_iterator<string>(cout, " ") );
   if(pSimplifiedTimeString){
      string simplified;
      for(int i=0; i<times.size(); ++i) simplified += times[i]+" ";
      *pSimplifiedTimeString=simplified;
   }

   int sum=0; //unit is second
   for(int i=0; i<times.size(); ++i){
      string digits, after_digits;
      separate_digits(times[i], digits, after_digits);
      after_digits=to_upper(strTrim(after_digits));
      int unit;
      switch(after_digits.at(0)){
         case 'S': unit=1; break;
         case 'M': unit=60; break;
         case 'H': unit=3600; break;
         default: cerr<<"you caught a bug\n"; exit(1);
      }//switch,
      sum+=string2int(digits)*unit;
   }
   
   return sum;
}
//--------------------------------------------------------------------------
#if !defined(MEX)
int main(int argc, char* argv[])
{
   if(argc!=2) {
      cout<<"usage: "<<argv[0]<<" a_time_string"<<endl; return 1;
   }
   string simplifiedTimeString;
   int totaltime=extract_time(argv[1], &simplifiedTimeString);
   cout<<simplifiedTimeString<<"="<<totaltime<<"s"<<endl;

   return 0;
}

#else
void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
   /* Check for proper number of input and  output arguments */    
   if (nrhs !=1) {
      mexErrMsgTxt("one input argument required.");
   }   
   if(nlhs >2 ){
      mexErrMsgTxt("0~2 output arguments required.");
   }
   
   /* input must be a string */
   if( mxIsChar(prhs[0]) != 1)  mexErrMsgTxt("Input must be a string.");
   
   //handle empty input string:
   //mexPrintf("matrix dimensions: %d x %d \n", mxGetM(prhs[0]), mxGetN(prhs[0]));
   if(mxGetM(prhs[0])*mxGetN(prhs[0])==0){
      mexErrMsgTxt("Input cannot be empty");
   }
   
   /* input must be a row vector */
   if(mxGetM(prhs[0])!=1)  mexErrMsgTxt("Input must be a row vector.");
   
   string input=matlab_arg::getStringArg(prhs[0]);
   string simplifiedTimeString;
   int totaltime=extract_time(input, &simplifiedTimeString);
   
   matlab_arg::setDoubleArg(plhs[0], totaltime);
   if(nlhs==2){
      matlab_arg::setStringArg(plhs[1], simplifiedTimeString);
   }
}

#endif
