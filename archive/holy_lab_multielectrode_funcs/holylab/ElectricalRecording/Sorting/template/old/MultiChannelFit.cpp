/* MultiChannelFit: fit templates to waveforms for multichannel
   electrical recording data

   Note: this version breaks compatibility with the Berry lab
   software! It is, however, a significant cleanup :-).

   Syntax:
      st = MultiChannelFit(waveforms,templates,template_shift_matrix,lambda,maxNumberOfOverlap,event_times)
      [st,residual_waveform] = MultiChannelFit(waveforms,templates,template_shift_matrix,lambda,maxNumberOfOverlap,event_times,template_offset)

   where

     waveforms is a n_channels-by-n_scans single-precision matrix
       containing a block of recorded data. (Note: unlike the Berry
       lab software, this is not overwritten, according to Matlab conventions)

     templates is a n_channels-by-len-by-n_templates single-precision
       3D array (len = # of scans in a template). You pretty much need
       to insure that the last template is an all-zeros template (the
       final template doesn't "count" and is used in termination
       criteria)

     template_shift_matrix is a n_templates-by-2 matrix, where
       template_shift_matrix(i,1) is the numeric tag associated with
       template i (usually its true "template #", before temporal
       shifting occurred) and template_shift_matrix(i,2) is the amount
       of temporal shift applied.

     lambda (a scalar) is the effective voltage penalty per added
       template (used to avoid over-fitting)

     maxNumberOfOverlap (a scalar) is upper bound on total number of
       templates that can be used to describe an event (?)

     event_times is a double-precision vector of event times (measured
       in scans) to be described in terms of the templates. Note this
       is unit-offset, so an event time of 1 corresponds to the first
       scan.
    
     template_offset (default 0) is a scalar, giving the number of
       samples into each template until the peak is reached (this
       keeps the final event times close to the input event times)

     st is a nspikes-by-2 matrix, where st(i,1) is the template
       numerical tag (taken from shift_matrix) of the ith spike, and
       st(i,2) is its time of occurrance (relative to the first scan
       in the waveform)
    
     residual is a single-precision n_channels-by-n_scans matrix of
       the residual waveform.
    
     Note: waveforms and templates must be of type single.  Note that
     event_times is double (less
     likely to get integer-overflow from long waveform samples).
*/

    // Copyright 2007 by Timothy E. Holy.  Based in part on MultiChannelFitV4 from the Berry lab, although it's essentially a total rewrite.
/* $Revision: 1.0 $ */
//#include "mex.h"


#ifdef MAIN
#include "mat.h"
#include "mat_variables.h"
#define mexPrintf printf
#include <error.h>
#include <errno.h>
#else
#include "mex.h"
#endif

//#include <math.h>
#include <vector>

void mcf_work(float *Data,int NumberOfChannels,int NumberOfSamples,const float *templates,int NumberOfTemplates,int LengthOfTemplates,const float *ShiftMatrix,int b,float Lambda,int NumberOfOverlap,const double *event_times,int NumberOfEvents,int template_offset,std::vector<int> &template_number,std::vector<double> &template_time);


/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  float *wave_out;
  const float *Data, *ShiftMatrix, *templates;
  const double *event_times;
  int NumberOfShifts,b,NumberOfOverlap,NumberOfChannels,NumberOfSamples;
  int NumberOfTemplates,LengthOfTemplates,NumberOfEvents;
  int template_offset;
  float Lambda;
  double *output_tmp;
  
  const mxArray *curarg;
  const mwSize *dims;
  
  if (nrhs < 6 || nrhs > 7) 
    mexErrMsgTxt("Six or seven inputs required.");
  if (nlhs < 1)
    mexErrMsgTxt("One or two outputs required.");
  
  // Parse the inputs
  // waveform data
  curarg = prhs[0];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg) ||
      mxGetNumberOfDimensions(curarg) != 2)
    mexErrMsgTxt("MultiChannelFit: waveform must be a real single-precision matrix");
  Data = (float *) mxGetData(curarg);
  NumberOfChannels = mxGetM(curarg);
  NumberOfSamples = mxGetN(curarg);
  
  // templates
  curarg = prhs[1];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg))
    mexErrMsgTxt("MultiChannelFit: templates must be a real single-precision 3D-array");
  if (mxGetNumberOfDimensions(curarg) != 3)
    mexErrMsgTxt("MultiChannelFit: template must be 3-dimensional, n_channels-by-len-by-n_templates");
  dims = mxGetDimensions(curarg);
  if (NumberOfChannels != dims[0])
    mexErrMsgTxt("MultiChannelFit: number of channels in templates does not agree with waveform");
  templates = (float*)mxGetData(curarg);
  NumberOfTemplates = dims[2];
  LengthOfTemplates = dims[1];
  
  // shift matrix
  curarg = prhs[2];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsSingle(curarg) ||
      mxGetNumberOfDimensions(curarg) != 2)
    mexErrMsgTxt("MultiChannelFit: shift_matrix must be a real single-precision matrix");
  if (mxGetM(curarg) != NumberOfTemplates)
    mexErrMsgTxt("MultiChannelFit: shift_matrix does not have one row per template");
  if (mxGetN(curarg) != 2)
    mexErrMsgTxt("MultiChannelFit: shift_matrix must be a n_templates-by-2 matrix");
  ShiftMatrix = (float*)mxGetData(curarg);
  
  // Lambda
  curarg = prhs[3];
  if (mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("MultiChannelFit: lambda must be a scalar");
  Lambda = (float)mxGetScalar(curarg);
  
  // NumberOfOverlap
  curarg = prhs[4];
  if (mxGetNumberOfElements(curarg) != 1)
    mexErrMsgTxt("MultiChannelFit: NumberOfOverlap must be a scalar");
  NumberOfOverlap = (int)mxGetScalar(curarg);
  
  // Event times
  curarg = prhs[5];
  if (!mxIsNumeric(curarg) || mxIsComplex(curarg) || !mxIsDouble(curarg) ||
      mxGetNumberOfDimensions(curarg) != 2)
    mexErrMsgTxt("MultiChannelFit: event_times must be a real double-precision vector");
  NumberOfEvents = mxGetNumberOfElements(curarg);
  event_times = mxGetPr(curarg);
  
  // template_offset
  if (nrhs > 6) {
    curarg = prhs[6];
    if (mxGetNumberOfElements(curarg) != 1)
      mexErrMsgTxt("MultiChannelFit: template_offset must be a positive integer scalar");
    template_offset = (int)mxGetScalar(curarg);
    if (template_offset < 0)
      mexErrMsgTxt("MultiChannelFit: template_offset must be a positive integer scalar");
  } else
    template_offset = 1;  // unit offset
  
  
  // Set up the output: we'll append spikes as we go
  std::vector<int> template_number;
  std::vector<double> template_time;
  
  //make it more efficient on vector usage
  size_t nPrealloc=size_t(NumberOfEvents*1.2); //20% extra for overlapping cases
  template_number.reserve(nPrealloc);
  template_time.reserve(nPrealloc);
  
  // Copy over the waveform data, to avoid overwriting the input
  if (nlhs > 1) {
    plhs[1] = mxCreateNumericMatrix(NumberOfChannels,NumberOfSamples,mxSINGLE_CLASS,mxREAL);
    wave_out = (float *) mxGetData(plhs[1]);
  } else
    wave_out = (float *) mxMalloc(NumberOfChannels*NumberOfSamples*sizeof(float));
  std::copy(Data,Data+NumberOfChannels*NumberOfSamples,wave_out);
  
  // Do the work
  mcf_work(wave_out,NumberOfChannels,NumberOfSamples,templates,NumberOfTemplates,LengthOfTemplates,ShiftMatrix,b,Lambda,NumberOfOverlap,event_times,NumberOfEvents,template_offset,template_number,template_time);
  
  // Copy over the results to a matlab array
  //mexPrintf("Number of spikes: %d\n",template_number.size());
  plhs[0] = mxCreateDoubleMatrix(template_number.size(),2, mxREAL);
  output_tmp = (double *) mxGetPr(plhs[0]);
  std::copy(template_number.begin(),template_number.end(),output_tmp);
  output_tmp += template_number.size();
  std::copy(template_time.begin(),template_time.end(),output_tmp);
  
  if (nlhs < 2)
    mxFree(wave_out);
}

//
// This function does the actual work of template fitting
//
void mcf_work(float *Data,int NumberOfChannels,int NumberOfSamples,const float *templates,int NumberOfTemplates,int LengthOfTemplates,const float *ShiftMatrix,int b,float Lambda,int NumberOfOverlap,const double *event_times,int NumberOfEvents,int template_offset,std::vector<int> &template_number,std::vector<double> &template_time) {
  // Currently the data formats are:
  // Data: waveform data as [all channels at time 1, all channels at time 2, ...]
  // Note: zero template is the LAST template
   
  int eventIndex, overlapIndex, templateIndex, best_template, n_overlaps;
  float err, smallest_err, err_tmp;
  float *residual_tmp, *scan_start, *scan_end;
  const float *template_tmp;
  bool keep_fitting;
  
  float *Residual, *Error, *Residual_end;
  int *Fit;
  double t_end,t_current;

  Residual = (float *)mxMalloc(NumberOfChannels*LengthOfTemplates*sizeof(float));
  Residual_end = Residual+NumberOfChannels*LengthOfTemplates;
  Error = (float *) mxMalloc(NumberOfOverlap*sizeof(float));
  Fit = (int *) mxMalloc(NumberOfOverlap*sizeof(int));
   
  // *** Loop through all possible spikes
  for (eventIndex=0;eventIndex<NumberOfEvents;eventIndex++) {
    // Copy over the waveform data into Residual
    scan_start = Data + size_t((event_times[eventIndex]-template_offset)*NumberOfChannels);
    scan_end = scan_start + NumberOfChannels*LengthOfTemplates;
    if (scan_start < Data || scan_end >= Data + NumberOfSamples*NumberOfChannels)
      continue;  // make sure we don't go beyond the edge of the data
    std::copy(scan_start,scan_end,Residual);
    
    // *** Initialize fit result arrays (Fit stores which waveform, error the error)
    for (overlapIndex=0;overlapIndex<NumberOfOverlap;overlapIndex++) {
      Fit[overlapIndex]=0;
      Error[overlapIndex]=0.;
    }
    
    // Loop over the # of overlaps, aborting if we fit the zero template
    overlapIndex = 0;
    keep_fitting = true;
    while (keep_fitting) {
      smallest_err = mxGetInf();
      best_template = -1;
      // Determine the template with the smallest fitting error
      for (templateIndex = 0, template_tmp = templates;
	   templateIndex < NumberOfTemplates; templateIndex++) {
	err = 0;
	// Calculate the fitting error
	for (residual_tmp = Residual; residual_tmp < Residual_end; residual_tmp++,template_tmp++) {
	  err_tmp = *residual_tmp - *template_tmp;
	  err += err_tmp*err_tmp;
	}
	if (err < smallest_err) {
	  smallest_err = err;
	  best_template = templateIndex;
	}
      }
      // Add the best template to the Fits array
      Fit[overlapIndex] = best_template;
      Error[overlapIndex] = smallest_err + 
	Lambda*Lambda*(overlapIndex+1)*(overlapIndex+1);
      if (best_template == NumberOfTemplates-1)
	keep_fitting = false;
      else if (overlapIndex > 0 &&
	       Error[overlapIndex] > Error[overlapIndex-1])
	keep_fitting = false;
      else {
	// *** Subtract best fit from waveform residual
	template_tmp = templates +
	  best_template*NumberOfChannels*LengthOfTemplates;
	for (residual_tmp = Residual; residual_tmp < Residual_end;
	     residual_tmp++, template_tmp++)
	  *residual_tmp -= *template_tmp;
	overlapIndex++;
	if (overlapIndex == NumberOfOverlap)
	  keep_fitting = false;
      }
    }  // while (keep_fitting)
    
    // For the templates that were centered at a time before the next
    // scheduled event, keep these as spike times. Note because the
    // fitting might assign a slightly different spike time than the
    // user supplied, we subtract 1 from the next event time to
    // prevent it from being "consumed" by this one. (If it does,
    // we've recorded it anyway; the only problem would be a small
    // "tail" that doesn't get properly subtracted because it's over
    // the right boundary of the shifted template.)
    n_overlaps = overlapIndex;
    if (eventIndex < NumberOfEvents-1)
      t_end = event_times[eventIndex+1]-1;  // t_end is at next event-1
    else
      t_end = NumberOfSamples;            // t_end is at end of recording
    if (n_overlaps > 0) {    // if we fit at least one template
      for (overlapIndex = 0; overlapIndex < n_overlaps; overlapIndex++) {
	t_current = ShiftMatrix[Fit[overlapIndex]+NumberOfTemplates] + 
	  event_times[eventIndex];
	if (t_current < t_end) {
	  // It occurred before the next event, keep it
	  template_number.push_back(1+int(Fit[overlapIndex]));
	  template_time.push_back(t_current);
	  // subtract this template permanently from the waveform, so
	  // we don't double-count it
	  template_tmp = templates +
	    Fit[overlapIndex]*NumberOfChannels*LengthOfTemplates;
	  for (residual_tmp = scan_start; residual_tmp < scan_end;
	       residual_tmp++,template_tmp++)
	    *residual_tmp -= *template_tmp;
	}
      }
    }  // if any templates were fit
  }  // loop over events
  
  mxFree(Residual);
  mxFree(Error);
  mxFree(Fit);
}


#ifdef MAIN   // this will compile only when we're creating the standalone MAT file for debugging
int main(int argc,char *argv[])
{
  char *filein,*fileout;
  const int n_inputs = 6;
  const int n_outputs = 1;
  bool use_offset = true;
  const mxArray *input[n_inputs+1];
  mxArray *output[n_outputs];
  const char *input_names[] = {
     "allWave", "T", "ShiftMatrix", "lambda", "maxNumberOfOverlap", "eventTime", "template_offset"
  };
  const char *output_names[] = {
    "ST"
  };

  if (argc < 3) {
    printf("Usage:\n  %s infile outfile\nwhere infile is a .mat file with variables g, psi1, psi2, psig, and sqrtdetJ\n",argv[0]);
    return EXIT_FAILURE;
  }

  filein = argv[1];
  fileout = argv[2];

  // Load the inputs
  mat_load_variables(filein,input_names,n_inputs,input);
  if (use_offset) {
    MATFile *matfp;

    matfp = matOpen(filein,"r");
    if (matfp == NULL)
      error(1,ENOENT,"Can't open file %s for reading",filein);
    input[n_inputs] = matGetVariable(matfp,"template_offset");
    if (input[n_inputs] == NULL) {
      printf("template_offset not defined");
      use_offset = false;
    }
    if (matClose(matfp) != 0)
      error(1,0,"Error closing file %s",filein);
  }

  // Call the C program, using the MEX interface function
  mexFunction(n_outputs,output,n_inputs+use_offset,input);

  // Save the outputs
  mat_save_variables(fileout,output_names,n_outputs,output);
  
  return EXIT_SUCCESS;
}
#endif

