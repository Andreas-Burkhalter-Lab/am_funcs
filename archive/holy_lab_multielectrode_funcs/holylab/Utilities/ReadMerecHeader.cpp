/*=================================================================
  This is a func similar to ReadAIHeader.
  It reads a header created by Merec and pass a struct back to matlab

  Note that: Not all fields in merec header are passed back---
     only the fields that appear in ReadAIHeader are passed back.

 *=================================================================*/

#include "parsemerecheader.h"

#include "mex.h"

//arg: int fd

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    /* Check for proper number of input and  output arguments */    
    if (nrhs !=1) {
        mexErrMsgTxt("one input argument (fd) required.");
    } 
    if(nlhs > 2){
        mexErrMsgTxt("0~2 output arguments required.");
    }

    //get fd:
    /* The input must be a noncomplex scalar double.*/
    int mrows = mxGetM(prhs[0]);
    int ncols = mxGetN(prhs[0]);
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !(mrows==1 && ncols==1) ) {
       mexErrMsgTxt("Input must be a noncomplex scalar double.");
    }
  
    /* get pointers to input arg */
    int fd = *mxGetPr(prhs[0]);

    
    TMerecHeader header;
    if(!parseMerecHeader(fd, header)){mexErrMsgTxt("error when open the merec data file");}
    
    /* Create a 1-by-1 array of structs. */ 
    int nfields=sizeof(fields)/sizeof(char*);
    plhs[0] = mxCreateStructMatrix(1, 1, nfields, fields);


    for(int ifield=0; ifield<nfields; ifield++) {
       mxArray *value;
       switch(ifield){
          case 0:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.headersize;
             break;
          }
          case 1:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.nscans;
             break;
          }
          case 2:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.numch;
             break;
          }
          case 3:{
             value = mxCreateDoubleMatrix(1,header.numch,mxREAL);
             double* p=mxGetPr(value);
             copy(header.channels.begin(), header.channels.end(), p);
             break;
          }
          case 4:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.scanrate;
             break;
          }
          case 5:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.scalemult;
             break;
          }
          case 6:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.scaleoff;
             break;
          }
          case 7:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.voltageMin;
             break;
          }
          case 8:{
             value = mxCreateDoubleMatrix(1,1,mxREAL);
             *mxGetPr(value) = header.voltageMax;
             break;
          }
          case 9:{
             value = mxCreateString(header.date.c_str());
             break;
          }
          case 10:{
             value = mxCreateString(header.time.c_str());
             break;
          }
          case 11:{
             value = mxCreateString(header.usrhdr.c_str());
             break;
          }
       }//switch, ifield

       mxSetFieldByNumber(plhs[0], 0, ifield, value);
       //                       //0: only 1 struct (struct 0)
    }

    if(nlhs==2){
       plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
       *mxGetPr(plhs[1]) = header.headersize;
    }

}

