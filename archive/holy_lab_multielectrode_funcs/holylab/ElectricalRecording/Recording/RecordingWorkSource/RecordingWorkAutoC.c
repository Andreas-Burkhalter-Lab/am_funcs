static char mc_version[] = "MATLAB Compiler 1.2 Jan  7 1998 infun";
/*
 *  MATLAB Compiler: 1.2
 *  Date: Jan  7 1998
 *  Arguments: -r RecordingWorkTryNostruct 

	The compiler-generated code was extensively modified, and facilties added
	to drive the National Instruments Card
	Tim Holy, 4/20/99
 */
#ifndef ARRAY_ACCESS_INLINING
#error You must use the -inline option when compiling MATLAB compiler generated code with MEX or MBUILD
#endif
#ifndef MATLAB_COMPILER_GENERATED_CODE
#define MATLAB_COMPILER_GENERATED_CODE
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "nidaq.h"
#include "nidaqcns.h"
#include "nidaqerr.h"
//#include "iotypes.h"
//#include "FileHeaders.cp"

// defines...
#define     kTrue             1
#define     kFalse            0

#define     errcheck( err, msg )                                  \
               {                                                  \
                  if( err )                                       \
                  {                                               \
                     printf( "Error %d from %s.\n", err, msg );   \
                     if( err < 0 )                                \
                        return;                                   \
                  }                                               \
               }

static void RecordingWork_fillmm2(mxArray *, mxArray *, mxArray *, mxArray *);
/* static array S0_ (1 x 20) text, line 21: '%d seconds remaining' */
static unsigned short S0__r_[] =
{
         37,  100,   32,  115,  101,   99,  111,  110,
        100,  115,   32,  114,  101,  109,   97,  105,
        110,  105,  110,  103,
};
static mxArray S0_ = mccCINIT( mccTEXT,  1, 20, S0__r_, 0);
/* static array S1_ (1 x 6) text, line 21: 'String' */
static unsigned short S1__r_[] =
{
         83,  116,  114,  105,  110,  103,
};
static mxArray S1_ = mccCINIT( mccTEXT,  1, 6, S1__r_, 0);
/* static array S2_ (1 x 6) text, line 27: 'String' */
static unsigned short S2__r_[] =
{
         83,  116,  114,  105,  110,  103,
};
static mxArray S2_ = mccCINIT( mccTEXT,  1, 6, S2__r_, 0);
/* static array S3_ (1 x 9) text, line 27: 'All done!' */
static unsigned short S3__r_[] =
{
         65,  108,  108,   32,  100,  111,  110,  101,
         33,
};
static mxArray S3_ = mccCINIT( mccTEXT,  1, 9, S3__r_, 0);
/* static array S4_ (1 x 3) text, line 29: 'Tag' */
static unsigned short S4__r_[] =
{
         84,   97,  103,
};
static mxArray S4_ = mccCINIT( mccTEXT,  1, 3, S4__r_, 0);
/* static array S5_ (1 x 13) text, line 29: 'StopRecButton' */
static unsigned short S5__r_[] =
{
         83,  116,  111,  112,   82,  101,   99,   66,
        117,  116,  116,  111,  110,
};
static mxArray S5_ = mccCINIT( mccTEXT,  1, 13, S5__r_, 0);
/* static array S6_ (1 x 6) text, line 29: 'Enable' */
static unsigned short S6__r_[] =
{
         69,  110,   97,   98,  108,  101,
};
static mxArray S6_ = mccCINIT( mccTEXT,  1, 6, S6__r_, 0);
/* static array S7_ (1 x 3) text, line 29: 'off' */
static unsigned short S7__r_[] =
{
        111,  102,  102,
};
static mxArray S7_ = mccCINIT( mccTEXT,  1, 3, S7__r_, 0);
mxArray *mlfFlipdim( mxArray *, mxArray * );
/* static array S8_ (1 x 1) text, line 35: 'b' */
static unsigned short S8__r_[] =
{
         98,
};
static mxArray S8_ = mccCINIT( mccTEXT,  1, 1, S8__r_, 0);
/* static array S9_ (1 x 9) text, line 35: 'EdgeColor' */
static unsigned short S9__r_[] =
{
         69,  100,  103,  101,   67,  111,  108,  111,
        114,
};
static mxArray S9_ = mccCINIT( mccTEXT,  1, 9, S9__r_, 0);
/* static array S10_ (1 x 1) text, line 35: 'b' */
static unsigned short S10__r_[] =
{
         98,
};
static mxArray S10_ = mccCINIT( mccTEXT,  1, 1, S10__r_, 0);
/* static array S11_ (1 x 1) text, line 35: 'b' */
static unsigned short S11__r_[] =
{
         98,
};
static mxArray S11_ = mccCINIT( mccTEXT,  1, 1, S11__r_, 0);
/* static array S12_ (1 x 9) text, line 35: 'EdgeColor' */
static unsigned short S12__r_[] =
{
         69,  100,  103,  101,   67,  111,  108,  111,
        114,
};
static mxArray S12_ = mccCINIT( mccTEXT,  1, 9, S12__r_, 0);
/* static array S13_ (1 x 1) text, line 35: 'b' */
static unsigned short S13__r_[] =
{
         98,
};
static mxArray S13_ = mccCINIT( mccTEXT,  1, 1, S13__r_, 0);
/***************** Compiler Assumptions ****************
 *
 *       B0_         	boolean scalar temporary
 *       BM0_        	boolean vector/matrix temporary
 *       I0_         	integer scalar temporary
 *       IM0_        	integer vector/matrix temporary
 *       IM1_        	integer vector/matrix temporary
 *       IM2_        	integer vector/matrix temporary
 *       R0_         	real scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       RM1_        	real vector/matrix temporary
 *       RM2_        	real vector/matrix temporary
 *       RecordingWork	<function being defined>
 *       RecordingWork/fillmm2	<function>
 *       TM0_        	string vector/matrix temporary
 *       axH         	real vector/matrix
 *       axes        	<function>
 *       buffersize  	real vector/matrix
 *       channels    	real vector/matrix
 *       counter     	integer scalar
 *       delete      	<function>
 *       drawnow     	<function>
 *       filename    	real vector/matrix
 *       findobj     	<function>
 *       gIsRecording	global real vector/matrix
 *       gcbf        	<function>
 *       i           	integer scalar
 *       minmax      	real vector/matrix
 *       npix        	real scalar
 *       nscans      	real vector/matrix
 *       patchH      	real vector/matrix
 *       rand        	<function>
 *       scanrate    	real vector/matrix
 *       set         	<function>
 *       sprintf     	<function>
 *       textH       	real vector/matrix
 *       xcoord      	integer vector/matrix
 *       zeros       	<function>
 *******************************************************/

void
mexFunction(
    int nlhs_,
    mxArray *plhs_[],
    int nrhs_,
    const mxArray *prhs_[]
)
{
   mxArray *Mplhs_[1];
   mxArray *Mprhs_[3];
   

   if (nrhs_ > 8 )
   {
      mexErrMsgTxt( "Too many input arguments." );
   }
   if (nrhs_ < 7 )
   {
      mexErrMsgTxt( "Too few input arguments." );
   }

   if (nlhs_ > 0 )
   {
      mexErrMsgTxt( "Too many output arguments." );
   }

   {
      mxArray channels;
      double scanrate;
      long nscans;
      long buffersize;
      int npix = 0;
      int scanrate_set,nscans_set,buffersize_set,npix_set;
      mxArray axH;
      mxArray textH;
      mxArray filenameMx;
      char *filename;
      int counter = 0;
      mxArray xcoord;
      mxArray patchH;
      int i = 0;
      mxArray gIsRecording;
      mxArray minmax;
      int I0_ = 0;
      mxArray IM0_;
      mxArray IM1_;
      mxArray RM0_;
      double R0_ = 0.0;
      mxArray BM0_;
      unsigned short B0_ = 0;
      mxArray IM2_;
      mxArray RM1_;
      mxArray RM2_;
      mxArray TM0_;
      int saving = 0;
      int status;
      //AIHeader aih;
      
      mccRealInit(channels);
      mccImport(&channels, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      scanrate = mccImportReal(&scanrate_set, 0, prhs_[1], 0, "scanrate");
      nscans = mccImportReal(&nscans_set, 0, prhs_[2], 0, "nscans");
      buffersize = mccImportReal(&buffersize_set, 0, prhs_[3], 0, "buffersize");
      npix = mccImportReal(&npix_set, 0, prhs_[4], 0, "npix");
      mccRealInit(axH);
      mccImport(&axH, ((nrhs_>5) ? prhs_[5] : 0), 0, 0);
      mccRealInit(textH);
      mccImport(&textH, ((nrhs_>6) ? prhs_[6] : 0), 0, 0);
      mccRealInit(filenameMx);
      mccImport(&filenameMx, ((nrhs_>7) ? prhs_[7] : 0), 0, 0);
      mccIntInit(xcoord);
      mccRealInit(patchH);
      mccRealInit(gIsRecording);
      mccRealInit(minmax);
      mccIntInit(IM0_);
      mccIntInit(IM1_);
      mccRealInit(RM0_);
      mccBoolInit(BM0_);
      mccIntInit(IM2_);
      mccRealInit(RM1_);
      mccRealInit(RM2_);
      mccTextInit(TM0_);
      
      /* % RecordingWorkTry(recparams,npix,axH,textH) */
      /* global gIsRecording */


	  // First check the input parameters
      if (!scanrate_set)
         mexErrMsgTxt( "variable scanrate undefined" );
      if (scanrate <= 0)
      	mexErrMsgTxt("Scanrate is <= 0");

      if (!nscans_set)
         mexErrMsgTxt( "variable nscans undefined" );
      if (nscans <= 0)
         mexErrMsgTxt( "nscans is <= 0" );


      if (!buffersize_set)
         mexErrMsgTxt( "variable buffersize undefined" );
      if (buffersize <= 0)
         mexErrMsgTxt( "buffersize is <= 0" );

      if (!npix_set)
         mexErrMsgTxt( "variable npix undefined" );
      if (npix <= 0)
         mexErrMsgTxt( "npix is <= 0" );
      
/*
	// Set up file if saving
	saving = 0;
	if (nrhs_ > 7) {
		if (!mxIsChar(filenameMx))
			mexErrMsgTxt("Filename must be a string");
		saving = 1;
		// First get length of filename
		filenamelen = (mxGetM(filenameMx)*mxGetN(filenameMx))+1;
		filename = mxCalloc(filenamelen,sizeof(char));
		status = mxGetString(filenameMx,filename,filenamelen);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, filename is truncated");
		// CONTINUE HERE
	}
*/
	// Set up handles for plotting
     /* xcoord = 1:npix; */
      mccColon2(&xcoord, (double)1, npix);
      /* patchH = zeros(1,64); */
      mccZerosMN(&patchH, 1, 64);
      /* for i = 1:64 */
      for (I0_ = 1; I0_ <= 64; I0_ = I0_ + 1)
      {
         i = I0_;
         // So delete(patchH) doesn't give an error, have to fill with a dummy set
         /* patchH(i) = fillmm2(zeros(1,npix),zeros(1,npix),xcoord); */
         mccZerosMN(&IM0_, 1, ((int)npix));
         mccZerosMN(&IM1_, 1, ((int)npix));
         RecordingWork_fillmm2(&RM0_, &IM0_, &IM1_, &xcoord);
         R0_ = (mccGetRealVectorElement(&RM0_, mccRint(1)));
         mccSetRealVectorElement(&patchH, mccRint(i), R0_);
         /* end */
      }
      
      // Now get the card ready to go
      
      
      
      /* while (gIsRecording & counter < nscans) */
      while (1)
      {
         mccGetGlobal(&gIsRecording, "gIsRecording");
         B0_ = mccIfCondition(&gIsRecording);
         if ((double)B0_)
         {
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               unsigned short *p_BM0_;
               int I_BM0_=1;
               double *p_nscans;
               int I_nscans=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&nscans), mccN(&nscans));
               mccAllocateMatrix(&BM0_, m_, n_);
               I_BM0_ = (mccM(&BM0_) != 1 || mccN(&BM0_) != 1);
               p_BM0_ = mccSPR(&BM0_);
               I_nscans = (mccM(&nscans) != 1 || mccN(&nscans) != 1);
               p_nscans = mccPR(&nscans);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_BM0_+=I_BM0_, p_nscans+=I_nscans)
                     {
                        *p_BM0_ = ( (counter < *p_nscans) && !mccREL_NAN(*p_nscans) );
                     }
                  }
               }
            }
            B0_ = mccIfCondition(&BM0_);
         }
         if (!((double)B0_))
            break;
         /* delete(patchH); */
         Mprhs_[0] = &patchH;
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "delete", 13);
         /* patchH = zeros(1,64); */
         mccZerosMN(&patchH, 1, 64);
         /* minmax = rand(64,2*npix); */
         Mprhs_[0] = mccTempMatrix(64, 0., mccINT, 0 );
         Mprhs_[1] = mccTempMatrix((2 * (double) npix), 0., mccINT, 0 );
         Mplhs_[0] = &minmax;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "rand", 15);
         /* minmax(:,xcoord) = -minmax(:,xcoord); */
         mccFindIndex(&IM0_, &xcoord);
         mccFindIndex(&IM1_, &xcoord);
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_minmax;
            int I_minmax=1;
            int *p_IM1_;
            int I_IM1_=1;
            double *p_1minmax;
            int I_1minmax=1;
            int *p_IM0_;
            int I_IM0_=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&minmax), (mccM(&IM0_) * mccN(&IM0_)));
            mccGrowMatrix(&minmax, m_, mccGetMaxIndex(&IM1_ ,mccN(&minmax)));
            mccCheckMatrixSize(&minmax, mccM(&minmax), mccGetMaxIndex(&IM0_ ,mccN(&minmax)));
            I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
            p_IM1_ = mccIPR(&IM1_);
            I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
            p_IM0_ = mccIPR(&IM0_);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_, p_IM1_ += I_IM1_, p_IM0_ += I_IM0_)
               {
                  p_minmax = mccPR(&minmax) + mccM(&minmax) * ((int)(*p_IM1_ - .5)) + 0;
                  p_1minmax = mccPR(&minmax) + mccM(&minmax) * ((int)(*p_IM0_ - .5)) + 0;
                  for (i_=0; i_<m_; ++i_, p_minmax+=I_minmax, p_1minmax+=I_1minmax)
                  {
                     *p_minmax = (-*p_1minmax);
                  }
               }
            }
         }
         /* for i = 1:64 */
         for (I0_ = 1; I0_ <= 64; I0_ = I0_ + 1)
         {
            i = I0_;
            /* axes(axH(i)); */
            if(mccNOTSET(&axH))
            {
               mexErrMsgTxt( "variable axH undefined, line 18" );
            }
            Mprhs_[0] = mccTempVectorElement(&axH, i);
            Mplhs_[0] = 0;
            mccCallMATLAB(0, Mplhs_, 1, Mprhs_, "axes", 18);
            /* patchH(i) = fillmm2(minmax(i,xcoord),minmax(i,xcoord+npix),xcoord); */
            mccFindIndex(&IM0_, &xcoord);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM0_;
               int I_RM0_=1;
               double *p_minmax;
               int I_minmax=1;
               int *p_IM0_;
               int I_IM0_=1;
               m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM0_) * mccN(&IM0_)));
               mccAllocateMatrix(&RM0_, m_, n_);
               mccCheckMatrixSize(&minmax, i, mccGetMaxIndex(&IM0_ ,mccN(&minmax)));
               I_RM0_ = (mccM(&RM0_) != 1 || mccN(&RM0_) != 1);
               p_RM0_ = mccPR(&RM0_);
               I_IM0_ = (mccM(&IM0_) != 1 || mccN(&IM0_) != 1);
               p_IM0_ = mccIPR(&IM0_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_IM0_ += I_IM0_)
                  {
                     p_minmax = mccPR(&minmax) + mccM(&minmax) * ((int)(*p_IM0_ - .5)) + (i-1);
                     for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_minmax+=I_minmax)
                     {
                        *p_RM0_ = *p_minmax;
                     }
                  }
               }
            }
            mccSTRING(&RM0_) = 0;
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               int *p_IM1_;
               int I_IM1_=1;
               int *p_xcoord;
               int I_xcoord=1;
               m_ = mcmCalcResultSize(m_, &n_, mccM(&xcoord), mccN(&xcoord));
               mccAllocateMatrix(&IM1_, m_, n_);
               I_IM1_ = (mccM(&IM1_) != 1 || mccN(&IM1_) != 1);
               p_IM1_ = mccIPR(&IM1_);
               I_xcoord = (mccM(&xcoord) != 1 || mccN(&xcoord) != 1);
               p_xcoord = mccIPR(&xcoord);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_)
                  {
                     for (i_=0; i_<m_; ++i_, p_IM1_+=I_IM1_, p_xcoord+=I_xcoord)
                     {
                        *p_IM1_ = (((int)*p_xcoord) + npix);
                     }
                  }
               }
            }
            mccFindIndex(&IM2_, &IM1_);
            {
               int i_, j_;
               int m_=1, n_=1, cx_ = 0;
               double *p_RM1_;
               int I_RM1_=1;
               double *p_minmax;
               int I_minmax=1;
               int *p_IM2_;
               int I_IM2_=1;
               m_ = mcmCalcResultSize(m_, &n_, 1, (mccM(&IM2_) * mccN(&IM2_)));
               mccAllocateMatrix(&RM1_, m_, n_);
               mccCheckMatrixSize(&minmax, i, mccGetMaxIndex(&IM2_ ,mccN(&minmax)));
               I_RM1_ = (mccM(&RM1_) != 1 || mccN(&RM1_) != 1);
               p_RM1_ = mccPR(&RM1_);
               I_IM2_ = (mccM(&IM2_) != 1 || mccN(&IM2_) != 1);
               p_IM2_ = mccIPR(&IM2_);
               if (m_ != 0)
               {
                  for (j_=0; j_<n_; ++j_, p_IM2_ += I_IM2_)
                  {
                     p_minmax = mccPR(&minmax) + mccM(&minmax) * ((int)(*p_IM2_ - .5)) + (i-1);
                     for (i_=0; i_<m_; ++i_, p_RM1_+=I_RM1_, p_minmax+=I_minmax)
                     {
                        *p_RM1_ = *p_minmax;
                     }
                  }
               }
            }
            mccLOG(&RM1_) = 0;
            mccSTRING(&RM1_) = 0;
            RecordingWork_fillmm2_1(&RM2_, &RM0_, &RM1_, &xcoord);
            R0_ = (mccGetRealVectorElement(&RM2_, mccRint(1)));
            mccSetRealVectorElement(&patchH, mccRint(i), R0_);
            /* end */
         }
         /* set(textH,'String',sprintf('%d seconds remaining',nscans-counter)); */
         if(mccNOTSET(&textH))
         {
            mexErrMsgTxt( "variable textH undefined, line 21" );
         }
         {
            int i_, j_;
            int m_=1, n_=1, cx_ = 0;
            double *p_RM2_;
            int I_RM2_=1;
            double *p_nscans;
            int I_nscans=1;
            m_ = mcmCalcResultSize(m_, &n_, mccM(&nscans), mccN(&nscans));
            mccAllocateMatrix(&RM2_, m_, n_);
            I_RM2_ = (mccM(&RM2_) != 1 || mccN(&RM2_) != 1);
            p_RM2_ = mccPR(&RM2_);
            I_nscans = (mccM(&nscans) != 1 || mccN(&nscans) != 1);
            p_nscans = mccPR(&nscans);
            if (m_ != 0)
            {
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM2_+=I_RM2_, p_nscans+=I_nscans)
                  {
                     *p_RM2_ = (*p_nscans - counter);
                  }
               }
            }
         }
         Mprhs_[0] = &S0_;
         Mprhs_[1] = &RM2_;
         Mplhs_[0] = &TM0_;
         mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "sprintf", 21);
         Mprhs_[0] = &textH;
         Mprhs_[1] = &S1_;
         Mprhs_[2] = &TM0_;
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 3, Mprhs_, "set", 21);
         /* counter = counter + 1; */
         counter = (counter + 1);
         /* %if (mod(counter,100) == 0) */
         /* drawnow */
         Mplhs_[0] = 0;
         mccCallMATLAB(0, Mplhs_, 0, Mprhs_, "drawnow", 24);
         /* %end */
         /* end */
      }
      /* set(textH,'String','All done!'); */
      if(mccNOTSET(&textH))
      {
         mexErrMsgTxt( "variable textH undefined, line 27" );
      }
      Mprhs_[0] = &textH;
      Mprhs_[1] = &S2_;
      Mprhs_[2] = &S3_;
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 3, Mprhs_, "set", 27);
      /* gIsRecording = 0; */
      {
         double tr_ = 0;
         mccAllocateMatrix(&gIsRecording, 1, 1);
         *mccPR(&gIsRecording) = tr_;
      }
      mccSetGlobal("gIsRecording", &gIsRecording);
      /* set(findobj(gcbf,'Tag','StopRecButton'),'Enable','off'); % enable Stop button */
      Mplhs_[0] = &RM2_;
      mccCallMATLAB(1, Mplhs_, 0, Mprhs_, "gcbf", 29);
      Mprhs_[0] = &RM2_;
      Mprhs_[1] = &S4_;
      Mprhs_[2] = &S5_;
      Mplhs_[0] = &RM1_;
      mccCallMATLAB(1, Mplhs_, 3, Mprhs_, "findobj", 29);
      Mprhs_[0] = &RM1_;
      Mprhs_[1] = &S6_;
      Mprhs_[2] = &S7_;
      Mplhs_[0] = 0;
      mccCallMATLAB(0, Mplhs_, 3, Mprhs_, "set", 29);
      /* return */
      goto MretR;
      
      MretR: ;
   }
   return;
}
/*
R = RecordingWork/fillmm2(RRI)
*/
/***************** Compiler Assumptions ****************
 *
 *       RM0_        	real vector/matrix temporary
 *       RecordingWork/fillmm2	<function being defined>
 *       fill        	<function>
 *       flipdim     	<function>
 *       h           	real vector/matrix
 *       vmax        	real vector/matrix
 *       vmin        	real vector/matrix
 *       x           	integer vector/matrix
 *       xx          	real vector/matrix
 *       yy          	real vector/matrix
 *******************************************************/
static void RecordingWork_fillmm2_1(mxArray *h_PP_, mxArray *vmin_P_, mxArray *vmax_P_, mxArray *x_P_ )
{
   mxArray *Mplhs_[1];
   mxArray *Mprhs_[5];
   mxArray h;
   mxArray yy;
   mxArray xx;
   mxArray RM0_;
   
   mccRealInit(h);
   mccRealInit(yy);
   mccRealInit(xx);
   mccRealInit(RM0_);
   
   /* yy = [vmin flipdim(vmax,2)]; */
   Mprhs_[0] = vmax_P_;
   Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
   Mplhs_[0] = &RM0_;
   mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "flipdim", 33);
   mccCatenateColumns(&yy, vmin_P_, &RM0_);
   /* xx = [x flipdim(x,2)]; */
   Mprhs_[0] = x_P_;
   Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
   Mplhs_[0] = &RM0_;
   mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "flipdim", 34);
   mccCatenateColumns(&xx, x_P_, &RM0_);
   /* h = fill(xx,yy,'b','EdgeColor','b'); */
   Mprhs_[0] = &xx;
   Mprhs_[1] = &yy;
   Mprhs_[2] = &S8_;
   Mprhs_[3] = &S9_;
   Mprhs_[4] = &S10_;
   Mplhs_[0] = &h;
   mccCallMATLAB(1, Mplhs_, 5, Mprhs_, "fill", 35);
   /* return */
   goto MretH;
   
   MretH:
   mccCopy(h_PP_, &h);
   mccFreeMatrix(&h);
   mccFreeMatrix(&yy);
   mccFreeMatrix(&xx);
   mccFreeMatrix(&RM0_);
   return;
}
/*
R = RecordingWork/fillmm2(III)
*/
/***************** Compiler Assumptions ****************
 *
 *       RM0_        	real vector/matrix temporary
 *       RecordingWork/fillmm2	<function being defined>
 *       fill        	<function>
 *       flipdim     	<function>
 *       h           	real vector/matrix
 *       vmax        	integer vector/matrix
 *       vmin        	integer vector/matrix
 *       x           	integer vector/matrix
 *       xx          	real vector/matrix
 *       yy          	real vector/matrix
 *******************************************************/
static void RecordingWork_fillmm2(mxArray *h_PP_, mxArray *vmin_P_, mxArray *vmax_P_, mxArray *x_P_ )
{
   mxArray *Mplhs_[1];
   mxArray *Mprhs_[5];
   mxArray h;
   mxArray yy;
   mxArray xx;
   mxArray RM0_;
   
   mccRealInit(h);
   mccRealInit(yy);
   mccRealInit(xx);
   mccRealInit(RM0_);
   
   /* yy = [vmin flipdim(vmax,2)]; */
   Mprhs_[0] = vmax_P_;
   Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
   Mplhs_[0] = &RM0_;
   mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "flipdim", 33);
   mccCatenateColumns(&yy, vmin_P_, &RM0_);
   /* xx = [x flipdim(x,2)]; */
   Mprhs_[0] = x_P_;
   Mprhs_[1] = mccTempMatrix(2, 0., mccINT, 0 );
   Mplhs_[0] = &RM0_;
   mccCallMATLAB(1, Mplhs_, 2, Mprhs_, "flipdim", 34);
   mccCatenateColumns(&xx, x_P_, &RM0_);
   /* h = fill(xx,yy,'b','EdgeColor','b'); */
   Mprhs_[0] = &xx;
   Mprhs_[1] = &yy;
   Mprhs_[2] = &S11_;
   Mprhs_[3] = &S12_;
   Mprhs_[4] = &S13_;
   Mplhs_[0] = &h;
   mccCallMATLAB(1, Mplhs_, 5, Mprhs_, "fill", 35);
   /* return */
   goto MretH;
   
   MretH:
   mccCopy(h_PP_, &h);
   mccFreeMatrix(&h);
   mccFreeMatrix(&yy);
   mccFreeMatrix(&xx);
   mccFreeMatrix(&RM0_);
   return;
}
