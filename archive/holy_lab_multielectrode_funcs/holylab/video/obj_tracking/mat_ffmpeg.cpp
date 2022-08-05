#include "errno.h"

#if defined(WIN32) || defined(__WIN32) || defined(WIN)
#   define _sys_errlist sys_errlist
#endif

#include <stdlib.h>

#include "mex.h"

#include "matlab_arg.h"
using namespace matlab_arg;
#include <vector>
using namespace std;

#include "ffmpeg.hpp"

//---------------------------------------------------------------------
enum {eOpen=1, eClose,
      eGetCurTime, 
      eNextFrame, eSeekFrame, eSeekPerfectFrame, 
      eGetHeight, eGetWidth,
      eGetDuration,
      eGetFrame,
      eGetErrorString,
      eGetFrameRate,
      eForceFrameRate,
};

vector<Ffmpeg*> mpgs;

   //	   mexPrintf("dmirror: initialized\n");
   //   mexErrMsgTxt(err.what());

void cleanup()
{
   for(unsigned i=0; i<mpgs.size(); ++i) delete mpgs[i];
   mpgs.clear();
}

//---------------------------------------------------------------------
void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
   if(nrhs<1) mexErrMsgTxt("At least 1 input argument required.");

   mexAtExit(cleanup);

   int action=getDoubleArg(prhs[0]);
   int id;
   switch(action){
      case eOpen: {
         string filename=getStringArg(prhs[1]);
         Ffmpeg * mpg=new Ffmpeg(filename.c_str());
         int result=-1;
         if(mpg->isOk) {
            mpgs.push_back(mpg);
            result=mpgs.size()-1;
         }
         //else: delete mpg (TODO: fixed Ffmpeg's dtor)
         setDoubleArg(plhs[0], result);
         break;
      }//case, eOpen
      case eClose:
         id=getDoubleArg(prhs[1]);
         //todo: check id is in range
         delete mpgs[id];
         mpgs.erase(mpgs.begin()+id);
         break;
      case eGetCurTime:{
         id=getDoubleArg(prhs[1]);
         double result=mpgs[id]->curTime();
         setDoubleArg(plhs[0], result);
         break;
      }
      case eNextFrame:
         id=getDoubleArg(prhs[1]);
         mpgs[id]->nextFrame();
         setDoubleArg(plhs[0], mpgs[id]->isOk);
         break;
      case eSeekFrame: {
         id=getDoubleArg(prhs[1]);
         double time=getDoubleArg(prhs[2]);
         mpgs[id]->seekFrame(time);
         setDoubleArg(plhs[0], mpgs[id]->isOk);
         break;
      }
      case eSeekPerfectFrame: {
         id=getDoubleArg(prhs[1]);
         double time=getDoubleArg(prhs[2]);
         mpgs[id]->seekFrame_perfect(time);
         setDoubleArg(plhs[0], mpgs[id]->isOk);
         break;
      }
      case eGetHeight:
         id=getDoubleArg(prhs[1]);
         setDoubleArg(plhs[0], mpgs[id]->height);
         break;
      case eGetWidth:
         id=getDoubleArg(prhs[1]);
         setDoubleArg(plhs[0], mpgs[id]->width);
         break;
      case eGetDuration:
         id=getDoubleArg(prhs[1]);
         setDoubleArg(plhs[0], mpgs[id]->duration);
         break;
      case eGetFrame: {
         id=getDoubleArg(prhs[1]);
         int nbytes=3*mpgs[id]->width*mpgs[id]->height;
         plhs[0]=mxCreateNumericMatrix(1, nbytes, mxUINT8_CLASS, mxREAL);
         memcpy(mxGetPr(plhs[0]), mpgs[id]->data, nbytes);
         break;
      }
      case eGetErrorString:
         id=getDoubleArg(prhs[1]);
         setStringArg(plhs[0], mpgs[id]->errMsg);
         break;
      case eGetFrameRate:{
         id=getDoubleArg(prhs[1]);
         double result=mpgs[id]->fps;
         setDoubleArg(plhs[0], result);
         break;
      }
      case eForceFrameRate:{
         id=getDoubleArg(prhs[1]);
         double fps=getDoubleArg(prhs[2]);
         mpgs[id]->fps=fps;
         break;
      }
   }//switch, action
   //      if(nrhs!=2)mexErrMsgTxt("2 input arguments required.");
}

