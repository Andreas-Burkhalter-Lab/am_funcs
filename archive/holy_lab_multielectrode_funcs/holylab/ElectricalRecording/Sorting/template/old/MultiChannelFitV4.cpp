

/*
 * =============================================================
 * xtimesy.c - example found in API guide
 *
 * Multiplies an input scalar times an input matrix and outputs a
 * matrix. 
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 1984-2000 The MathWorks, Inc.
 * =============================================================
 */

/* $Revision: 1.9 $ */
#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* To transfer a matrix A from Matlab where [m,n]=size(A) then:
   A(i,j)=*(ptr+m*(j-1)+(i-1))
   
   ptr is a pointer
*/





#define MaximumChannelNumber 30
#define MaximumNumberOfSpikes 40000
#define MaximumNumberOfOverlap 20

const int cMaxNumOfShiftMatrixRows=20521; //that is, at most 360 templates if shift is in [-28 28] 

float unarrangedT[cMaxNumOfShiftMatrixRows][MaximumChannelNumber*64];

//double ShiftMatrix[2500][4];
//double T[2500][MaximumChannelNumber*64];






// ST=MultiChannelFitV1(Data,T,ShiftMatrix,Lambda,NumberOfOverlap,I);
/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   float *ShiftMatrix,*Data,*I;
   double *ST;
   int WaveformNumber,NumberOfShifts,b,NumberOfOverlap,NumberOfChannels,NumberOfSamples;
   float Lambda;
   
   float *InT;
   float **rearrangedT;
   float **T;
   
   //  float *TransformedTemplates;
   int i,j,k,count,Index=1;
   // mexErrMsgTxt("Function beginning");
   
   
   //  mexPrintf("MEX-file executing %d %d %d %d\n",NumberOfShifts,NumberOfChannels,NumberOfSamples,MaximumChannelNumber*64);
   /*  Check for proper number of arguments. */
   /* NOTE: You do not need an else statement when using
      mexErrMsgTxt within an if statement. It will never
      get to the else statement if mexErrMsgTxt is executed.
      (mexErrMsgTxt breaks you out of the MEX-file.) 
   */
   if (nrhs != 6) 
      mexErrMsgTxt("Six inputs required.");
   
   
   
   
   /* Create a pointer to the input matrix y. */
   Data = (float*)mxGetPr(prhs[0]);    
   NumberOfChannels = mxGetM(prhs[0]);
   NumberOfSamples = mxGetN(prhs[0]);
   
   
   InT = (float*)mxGetPr(prhs[1]); 
   
   //  TransformedTemplates=malloc(sizeof(float)*mxG
   
   
   ShiftMatrix = (float*)mxGetPr(prhs[2]);    
   NumberOfShifts = mxGetM(prhs[2]);
   b = mxGetN(prhs[2]);
   
   
   //  mexPrintf("MEX-file executing %d %d %d %d\n",NumberOfShifts,NumberOfChannels,NumberOfSamples,MaximumChannelNumber*64);
   
   rearrangedT=(float **)malloc(sizeof(float *)*NumberOfShifts);
   
   for(i=0;i<NumberOfShifts;i++)
      rearrangedT[i]=(float *)malloc(sizeof(float)*NumberOfChannels*64);
   
   for (i=1;i<=NumberOfShifts;i++)
      for (j=1;j<=NumberOfChannels*64;j++)
         unarrangedT[i-1][j-1]=(*(InT+NumberOfShifts*(j-1)+(i-1)));
   
   // Now rearrange T
   for(i=0;i<NumberOfShifts;i++)
      for(j=0;j<NumberOfChannels;j++)
         for(k=0;k<64;k++)
            rearrangedT[i][k*NumberOfChannels+j]=unarrangedT[i][j*64+k];
   
   // Choose which ordering of the templates to use
   // T=unarrangedT;
   T=rearrangedT;
   
   //mexPrintf("MEX-file executing %f\n",Lambda);
   
   
   
   Lambda = (float)mxGetScalar(prhs[3]);
   
   NumberOfOverlap = (int)mxGetScalar(prhs[4]);
   
   I = (float*)mxGetPr(prhs[5]);
   WaveformNumber = mxGetN(prhs[5]);
   
   
   // Set the output pointer to the output matrix. 
   plhs[0] = mxCreateDoubleMatrix(MaximumNumberOfSpikes,2, mxREAL);
   //plhs[1] = mxCreateDoubleMatrix(WaveformNumber,NumberOfOverlap+1, mxREAL);
   // Create a C pointer to a copy of the output matrix. 
   ST = (double *)mxGetPr(plhs[0]);
   
   
   
   //  mexPrintf("MEX-file executing %f %d %d %d %d\n",Lambda,WaveformNumber,NumberOfShifts,NumberOfChannels,NumberOfSamples);
   
   
   // Currently the data formats are:
   // Data: waveform data as [all channels at time 1, all channels at time 2, ...]
   // T: template data as [all times at channel 1, all times at channel 2,...] (all times being the 64 time samples from -32:31)
   // Note: zero template is the LAST template in T
   
   
   
   // *** Loop through all possible spikes
   for (i=1;i<=WaveformNumber;i++)
   {
      float Waveform[MaximumChannelNumber*64];
      float Error[MaximumNumberOfOverlap],WaveformToRemove[MaximumChannelNumber*64];
      float BestFitScore,Score,MinFitValue;
      int k,n,m,BestFit,NN,MinFit;
      int Fit[MaximumNumberOfOverlap];
      
      
      // *** Initialize fit result arrays (Fit stores which waveform, error the error)
      for (k=0;k<NumberOfOverlap;k++) 
      {
         Fit[k]=0;
         Error[k]=0.;
      }
      
      m=0;
      
      // This option is to transform the data
      // *** Transform waveform data from [all channels at time 1, all channels at time 2,...] to [all times for channel 1, all times for channel 2,...]
      /*		for (n=1;n<=NumberOfChannels;n++)
			for (k=-32;k<32;k++)
			{
                        Waveform[m]=*(Data+((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1));
                        
                        m++;
			}*/
      
      // This option just copies the data
      for (k=-32;k<32;k++)
         for (n=1;n<=NumberOfChannels;n++)
         {
            Waveform[m]=*(Data+((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1));
            
            m++;
         }
      
      
      // We start trying the first fit
      NN=0;
      
      do
      {            
         // *** Calculate mean squared error between the zero template and the waveform under investigation
         BestFitScore=0.0;
         
         for (n=1;n<=NumberOfChannels*64;n++)
            BestFitScore+=(Waveform[n-1]-T[NumberOfShifts-1][n-1])*(Waveform[n-1]-T[NumberOfShifts-1][n-1]);
         
         BestFit=NumberOfShifts-1;
         
         
         // *** Loop through all possible templates, checking if the mean squared error is less than the current best
         for (k=0;k<(NumberOfShifts-1);k++)
         {
            Score=0.0;
            
            for (n=NumberOfChannels*24;n<NumberOfChannels*64;)
            {
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
               
               if(Score>BestFitScore)
                  n=NumberOfChannels*64+1;
            }
            
            if(n==NumberOfChannels*64)
               for (n=0;n<NumberOfChannels*24;)
               {
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  Score+=(Waveform[n]-T[k][n])*(Waveform[n]-T[k][n]); n++;
                  
                  if(Score>BestFitScore)
                     n=NumberOfChannels*64+1;
               }
            
            if (Score<BestFitScore)
            {
               BestFit=k;
               BestFitScore=Score;
            }
         }
         
         // *** Subtract best fit from waveform, and add it to the fits array
         for (n=1;n<=NumberOfChannels*64;n++)
            Waveform[n-1]=Waveform[n-1]-T[BestFit][n-1];
         
         Fit[NN]=BestFit+1;
         Error[NN]=sqrt(BestFitScore)+Lambda*(NN+1);
         
         NN++;
      }
      while( (NN<NumberOfOverlap) && (Fit[NN-1]!=NumberOfShifts) );
      // keep repeating until we have reached NumberOfOverlap fits, or have fit the zero template
      
      // If we broke out by fitting the zero template, fill in the rest of the 'fits' with large error zero fits.
      if(NN<NumberOfOverlap)
         for(j=NN;j<NumberOfOverlap;j++)
         {
            Fit[j]=NumberOfShifts;
            Error[j]=Lambda*1000; // Massive error; these fits are definitely not used later on
         }
      
      /* **************************************     
         
      //for (k=0;k<NumberOfChannels*64;k++)
      //	mexPrintf("%f ",Waveform[k]);
      
      // *** Calculate mean squared error between the first template and the waveform under investigation
      BestFitScore=0.;
      for (n=1;n<=NumberOfChannels*64;n++)
      BestFitScore+=(Waveform[n-1]-T[0][n-1])*(Waveform[n-1]-T[0][n-1]);
      BestFitScore=sqrt(BestFitScore);
      
      //	mexPrintf("\n %f \n",BestFitScore);
      
      //for (k=1;k<=NumberOfChannels*64;k++)
      //	mexPrintf("%f ",*(T+(k-1)*NumberOfShifts));
      
      
      // *** Loop through all possible templates, checking if the mean squared error is less than the current best
      for (k=2;k<=NumberOfShifts;k++)
      {
      Score=0.;
      for (n=1;n<=NumberOfChannels*64;n++)
      Score+=(Waveform[n-1]-T[k-1][n-1])*(Waveform[n-1]-T[k-1][n-1]);
      Score=sqrt(Score);
      
      if (Score<BestFitScore)
      {
      BestFit=k;
      BestFitScore=Score;
      //mexPrintf(" %d %f %f\n",BestFit,Score,BestFitScore);
      }
      
      //mexPrintf(" %f ",Score);
      }
      
      //	mexPrintf("\n %f  %d \n ",BestFitScore,BestFit);
      
      
      // *** Subtract best fit from waveform, and add it to the fits array
      for (n=1;n<=NumberOfChannels*64;n++)
      Waveform[n-1]=Waveform[n-1]-T[BestFit-1][n-1];
      
      Fit[0]=BestFit;
      Error[0]=BestFitScore+Lambda;
      
      //     mexPrintf(" %d %f %f\n",Fit[0],Score,Error[0]);
      
      //mexPrintf(" %d    %f    %f    \n",BestFit,Score,BestFitScore);
      
      
      // *** If the best fit WASN'T the zero template...
      if (BestFit<NumberOfShifts)
      {
      
      // *** Loop through all earlier stuff for overlapping waveforms
      // Note: apparently we can't break out of this - even if we are constantly matching zero we just keep going until we reach NumberOfOverlap fits
      // In fact we don't even TRY and match the zero template (k only ranges to NumberOfShifts-1). Why not? 
      for (NN=2;NN<=NumberOfOverlap;NN++)
      {
      BestFit=1;
      BestFitScore=0.;
      for (n=1;n<=NumberOfChannels*64;n++)
      BestFitScore+=(Waveform[n-1]-T[0][n-1])*(Waveform[n-1]-T[0][n-1]);
      BestFitScore=sqrt(BestFitScore);
      for (k=2;k<=NumberOfShifts;k++)
      {
      Score=0.;
      for (n=1;n<=NumberOfChannels*64;n++)
      Score+=(Waveform[n-1]-T[k-1][n-1])*(Waveform[n-1]-T[k-1][n-1]);
      Score=sqrt(Score);
      
      if (Score<BestFitScore)
      {
      BestFit=k;
      BestFitScore=Score;
      //mexPrintf(" %d %f %f\n",BestFit,Score,BestFitScore);
      }
      
      //mexPrintf(" %f ",Score);
      }
      
      if(BestFit==NumberOfShifts) // matched zero!
      {
      // Fill in fits array for remaining fits to do
      for(j=NN;j<=NumberOfOverlap;j++)
      {
      Fit[NN-1]=NumberOfShifts;
      Error[NN-1]=Lambda*1000; // Massive error; these fits are definitely not used later on
      }
      
      mexPrintf("Breaking out after %d\n",NN);
      // Break out of loop
      NN=NumberOfOverlap+1;
      
      }
      
      for (n=1;n<=NumberOfChannels*64;n++)
      Waveform[n-1]=Waveform[n-1]-T[BestFit-1][n-1];				
      
      Fit[NN-1]=BestFit;
      
      // *** The error for each fit is the root mean square difference between the waveform and the template, plus a term corresponding to how many fits have so far been done
      // If we have a perfectly clean waveform that is the sum of two of our templates, we could identify both of them but report a large error, despite the fact the waveforms sum together to fit the data perfectly. Maybe not the best method to calculate an error - maybe it should be recalculated after all the fits are done, for combinations of fits.
      Error[NN-1]=BestFitScore+NN*Lambda;
      }
      
      ************************** */
      
      //mexPrintf("########### %d %f %f %d %f \n",BestFit,Score,BestFitScore,Fit[NN-1],Error[NN-1]);
      
      if (NN>1) // if we fit at least one waveform
      {
         
         for (n=0;n<NumberOfChannels*64;n++)
            WaveformToRemove[n]=0;
         
         // *** MinFit is really a MaxFit - the last (maximum index) fit that 'fits' (by 'fits' we mean it has less error than any previous fit)
         MinFit=0;
         MinFitValue=Error[0];
         for (n=1;n<NumberOfOverlap;n++)
            if (Error[n]<MinFitValue)
            {
               MinFit=n;
               MinFitValue=Error[n];
            }
         //mexPrintf("########### %d \n",MinFit);
         
         // *** Now for all the fits up to MinFit (so all the hopefully good fits) we look in the second column of the ShiftMatrix
         // for the corresponding shift (and require it to be less than 20 and greater than -20 - why 20? max(SM(:,2)) reports 28;
         //  are we just ignoring 30% of the possible shifts?). If it is less than 20 then we add the fitted waveform to the WaveformToRemove.
         
         for (k=0;k<=MinFit;k++)				
            if (fabs(*(ShiftMatrix+NumberOfShifts+(Fit[k]-1)))<20)
            {
               for (n=1;n<=NumberOfChannels*64;n++)
                  WaveformToRemove[n-1]=WaveformToRemove[n-1]+T[Fit[k]][n-1];
               
               // *** ST (presumable meaning SpikeTimes) is the output list of spike times. Index is initialized to one at the beginning
               // of the code and incremented by one each time here. ST is a 2*MaximumNumberOfSpikes array. MaximumNumberOfSpikes is
               // hard coded (a #define at the top of the code). Could fail for large files? Each time here we store the time and the fit
               // (the fit being out of numberOfCells*numberOfShifts)
               
               *(ST+(Index-1))=*(I+(i-1));
               *(ST+MaximumNumberOfSpikes+(Index-1))=Fit[k];
               Index=Index+1;
               //mexPrintf(" %d %f ",Fit[k],Error[k]);
            }
         
         // *** Finally remove WaveformToRemove from the original data
         m=0;
         
         // for unarrangedT
         /*	for (n=1;n<=NumberOfChannels;n++)
                for (k=-32;k<32;k++)
                {
                //mexPrintf("%d\n ",((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1));
                *(Data+((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1))=*(Data+((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1))-WaveformToRemove[m];
                m++;
         */
         
         // for rearrangedT
         for (k=-32;k<32;k++)
            for (n=1;n<=NumberOfChannels;n++)
            {
               //mexPrintf("%d\n ",((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1));
               *(Data+((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1))=*(Data+((int)(*(I+i-1)+k-1))*NumberOfChannels+(n-1))-WaveformToRemove[m];
               m++;
            }
         
      }
      
   }
   
   for(i=0; i<NumberOfShifts; i++)
      free(rearrangedT[i]);
   
   free(rearrangedT);
   
}  

