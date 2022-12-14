/* $Revision: 1.1 $ */
// cedFunction.cpp
// Automatically generated by Matlab AppWizard version 1.0
//
// This is the gateway routine for a MATLAB Math/Graphics Library-based
// C MATLAB MEX File.

#include "main.h"


void son_getEventData(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
	long numEvents,			// Number of returned events
		 max;				// Max # of events to return.
	short handle;			// SON file handle
	WORD chan;				// channel
	TSTime sTime,			// start time
		   eTime;			// end time
	TpSTime evtDataP;		// Pointer to event data buffer
	BOOLEAN tpb;
	int dims[2];
	TpSTime ptr;
	int i;
	char **fieldNames;

	// Check the number of right-hand-side arguments.
	if (nrhs != 6) {
		mexErrMsgTxt("Incorrect Number of Inputs for SonGetEventData.");
	}
	// Check the number of left-hand-side arguments.
	else if (nlhs != 2) {
		mexErrMsgTxt("Incorrect number of outputs for SonGetEventData.");
	}
	else {
		// Get inputs
		handle = (short)mxGetScalar(prhs[1]);
		chan = (WORD)mxGetScalar(prhs[2]);
		max = (long)mxGetScalar(prhs[3]);
		sTime = (TSTime)mxGetScalar(prhs[4]);
		eTime = (TSTime)mxGetScalar(prhs[5]);

		// Initialize fieldNames, which will describe the the structure entries in the mxArray of
		// structures returned to Matlab.
		fieldNames = (char**)malloc(1 * sizeof(char*));
		fieldNames[0] = (char*)malloc(25 * sizeof(char));
		strcpy(fieldNames[0], "time");

		// Grab event data in the channel.
		if ((evtDataP = (TpSTime)malloc(max * sizeof(TSTime))) == NULL) {
			mexErrMsgTxt("Malloc Failure in GetEventData!\n");
			return;
		}
		numEvents = SONGetEventData(handle, chan, evtDataP, max, sTime, eTime, &tpb, NULL);

		// Do error checking on the function.
		if (numEvents < 0) {
			son_error(numEvents);
			return;
		}

		// Copy the data into a size-optimized structure array.
		ptr = (TpSTime)malloc(numEvents * sizeof(TSTime));
		memcpy(ptr, evtDataP, numEvents * sizeof(TSTime));

		// Store the number of events we found.
		dims[0] = 1; dims[1] = 1;
		plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		*mxGetPr(plhs[0]) = (double)numEvents;

		// Create the mxArray that will hold the structure array values.
		dims[1] = numEvents;
		plhs[1] = mxCreateStructArray(2, dims, 1, fieldNames);

		// Copy all the data over.  I do it this way because I had problems with casting
		// a long to a double.
		for (i = 0; i < numEvents; i++) {
			mxArray *eventTime = mxCreateDoubleMatrix(1, 1, mxREAL);
			mxGetPr(eventTime)[0] = ptr[i];
			mxSetFieldByNumber(plhs[1], i, 0, eventTime);
		}

		free(evtDataP);
		free(ptr);
	}
}


void son_chanKind(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs) {

	short handle;	// file handle
	WORD chan;		// channel
	int dims[2] = {1,1};

	if (nrhs != 3) {
		mexErrMsgTxt("Incorrect Number of Inputs for SonChanKind.");
		return;
	}
	else if (nlhs != 1) {
		mexErrMsgTxt("Incorrect number of outputs for SonChanKind.");
		return;
	}
	else {
		// Get inputs.
		handle = (short)mxGetScalar(prhs[1]);
		chan = (WORD)mxGetScalar(prhs[2]);

		// Call SonChanKind
		plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		*mxGetPr(plhs[0]) = (double)SONChanKind(handle, chan);
	}
}





/*************************************************************************/
/*	void son_getAllADCData(int, mxArray**, int, const mxArray**)		 */
/*																		 */
/*************************************************************************/
void son_getAllADCData(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {

	if (nrhs == 7)
	{	
		short handle = (short)mxGetScalar(prhs[1]);	// file handle
		WORD chan =  (WORD)mxGetScalar(prhs[2]);	// channel to be read in
		long max = (long)mxGetScalar(prhs[3]);		// max number of pts to be read in
		long sTime = (long)mxGetScalar(prhs[4]);	// start time for reading in data
		long eTime = (long)mxGetScalar(prhs[5]);	// end time to read in data
		int ticks = (int)mxGetScalar(prhs[6]);		// number of time units per ADC pt
		int dims[2];
		short *AdcData, *data_calloc, *data_ptr;
		int StartCounter;
		long num_pts_read;
		double* np_ptr;
		TpSTime pbTime = (TpSTime)malloc(sizeof(TSTime));
		
		//create an mxArray to hold the maximum number of data pts
		dims[0] = 1;
		dims[1] = eTime/ticks;
		plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS,mxREAL);
		AdcData = (short*)mxGetPr(plhs[0]);
		
		//initialize a temp array to hold numpts worth of data
		data_calloc = (short*)calloc(max, sizeof(short));
		
		data_ptr = data_calloc;
		
		//counter variables to be used within the for loop
		StartCounter = 0;
		num_pts_read = 0;
		
		//while the number of points read is less than the total number of pts
		while (num_pts_read < max)
		{
			int counter, i;
			
			long np = SONGetADCData(handle, chan, data_ptr, max, sTime, eTime, pbTime, NULL);
			
			// advance the starting point for data acquisition by the appropriate 
			// increment
			sTime= *pbTime + np * ticks;
			
			//insert data into output array
			counter = 0;
			StartCounter = *pbTime / ticks;
			for (i=StartCounter;i<StartCounter + np; i++)
			{
				AdcData[i] = data_calloc[counter];
				counter++;
			}
			
			num_pts_read += np;
		}
		
		//create an mxArray to hold the number of data pts retrieved
		dims[0] = 1;
		dims[1] = 1;
		plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		np_ptr = mxGetPr(plhs[1]);
		np_ptr[0] = (double)num_pts_read;
		
		free((short *)data_calloc);
		data_calloc = NULL;
	}
	else {
		mexErrMsgTxt("Incorrect Number of Inputs for SonGetAllADCData.");
	}
}

//**********************************************************************************//
//	void son_getusPerTime(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)	//
//**********************************************************************************//
void son_getusPerTime(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
	int dims[2] = {1,1};
	double* uspertime_ptr;
	short handle;
	WORD uspertime;

	if (nrhs != 2) {
		mexErrMsgTxt("Incorrect Number of Inputs for SONGetusPerTime.");
	}
	else if (nlhs != 1) {
		mexErrMsgTxt("Incorrect Number of Outputs for SONGetusPerTime.");
	}
	else {
		// File handle
		handle = (short)mxGetScalar(prhs[1]);

		uspertime = SONGetusPerTime(handle);
		
		// Create an mxArray to hold the number of microsecs per time inter
		plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		*mxGetPr(plhs[0]) = (double)uspertime;
	}
}


//**********************************************************************************//
// void son_getNumADCPts(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)	//
//**********************************************************************************//
void son_getNumADCPts(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
	if (nrhs == 7)
	{
		TpAdc data;
		int dims[2] = {1,1};
		double* numpts_ptr;
		short handle = (short)mxGetScalar(prhs[1]);		//file handle
		WORD chan = (WORD)mxGetScalar(prhs[2]);			//channel to be read in
		long max = (long)mxGetScalar(prhs[3]);			//max number of pts to be read in
		TSTime sTime = (TSTime)mxGetScalar(prhs[4]);	//start time for reading in data
		TSTime eTime = (TSTime)mxGetScalar(prhs[5]);	//end time to read in data
		TSTime ticks = (TSTime)mxGetScalar(prhs[6]);	//number of time units per ADC pt
		TpSTime pbTime = (TpSTime)malloc(sizeof(TSTime));
		TSTime np, numpts = 0;

		data = (TpAdc)malloc(max * sizeof(short));

		while (sTime < eTime)
		{	
			np = SONGetADCData(handle, chan, data, max, sTime, eTime, pbTime, NULL);
			sTime = *pbTime + (np * ticks);
			numpts += np;
		}
		
		// create an mxArray to hold the number of data pts
		plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		numpts_ptr = mxGetPr(plhs[0]);
		numpts_ptr[0] = (double)numpts;
	}
	else
	{
		mexErrMsgTxt("Incorrect Number of Inputs for SONGetNumADCPts.");
	}
}



/*************************************************************************/
/*	void son_getADCData(int, mxArray**, int, const mxArray**)			 */
/*																		 */
/*************************************************************************/
void son_getADCData(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {

	if (nrhs == 6)
	{ 
		long dims[2], np;
		double *np_ptr, *pbTime_ptr;
		short handle = (short)mxGetScalar(prhs[1]);		// file handle
		WORD chan = (WORD)mxGetScalar(prhs[2]);			// channel to be read in
		long max = (long)mxGetScalar(prhs[3]);			// max number of pts to be read in
		TSTime sTime = (TSTime)mxGetScalar(prhs[4]);	// start time for reading in data
		TSTime eTime = (TSTime)mxGetScalar(prhs[5]);	// end time to read in data
		TpSTime pbTime = (TpSTime)malloc(sizeof(TSTime));
		TpAdc AdcData;
		
		//create an mxArray to hold the maximum number of data pts
		dims[0] = 1;
		dims[1] = max;
		plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS,mxREAL);
		AdcData = (TpAdc)mxGetPr(plhs[0]);
		
		np = SONGetADCData(handle, chan, AdcData, max, sTime, eTime, pbTime, NULL);
		
		// create an mxArray to hold the number of data pts retrieved
		dims[0] = 1;
		dims[1] = 1;
		plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		np_ptr = mxGetPr(plhs[1]);
		np_ptr[0] = (double)np;
		
		// create an mxArray to hold the time of the first data point
		dims[0] = 1;
		dims[1] = 1;
		plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		pbTime_ptr = mxGetPr(plhs[2]);
		pbTime_ptr[0] = (double)(*pbTime);
	}
	else
	{
		mexErrMsgTxt("Incorrect Number of Inputs for SONGetAdcData.");
	}
}

/*********************************************************************/
/* void son_chanDivide(int, mxArray**, int, const mxArray**)		 */
/*																	 */
/*********************************************************************/
void son_chanDivide(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {

	if (nrhs == 3)
	{
		double* output;
		short handle = (short)mxGetScalar(prhs[1]);		// file handle
		WORD chan = (WORD)mxGetScalar(prhs[2]);			// channel
		
		TSTime chan_interval = SONChanDivide(handle, chan);

		/* create an mxArray for output of basic time units per waveform conversion */
		int dims[2] = {1,1};
		plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		output = mxGetPr(plhs[0]);
		output[0] = (double)chan_interval;
	}
	else {
		mexErrMsgTxt("Incorrect Number of Inputs for SONChanDivide.");
	}
}


/*******************************************************************************/
/* void son_chanMaxTime(int, mxArray**, int, const mxArray**)				   */
/*																			   */
/*******************************************************************************/
void son_chanMaxTime(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {

	if (nrhs == 3)
	{
		short handle = (short)mxGetScalar(prhs[1]);		// file handle
		WORD chan = (WORD)mxGetScalar(prhs[2]);			// channel
		int dims[2] = {1,1};
		double* output;
		TSTime maxTime;

		// Get the max time from the smr file.
		maxTime = SONChanMaxTime(handle, chan);

		/* create an mxArray for output of maximum time of channel */
		plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		output = mxGetPr(plhs[0]);
		output[0] = (double)maxTime;
	}
	else
	{
		mexErrMsgTxt("Incorrect Number of inputs for SONChamMaxTime.");
	}
}


//**********************************************************************************//
//	void son_error(int error_no)													//
//**********************************************************************************//
void son_error(int error_no)
{
	switch (error_no)
	{
	case SON_NO_FILE:
		mexPrintf("SON File error: unknown file\n");
		break;
	case SON_NO_DOS_FILE:
		mexPrintf("SON File error: file doesn't exist\n");
		break;
	case SON_NO_PATH:
		mexPrintf("SON File error: path to file doesn't exist\n");
		break;
	case SON_NO_HANDLES:
		mexPrintf("SON File error: too many open files\n");
		break;
	case SON_NO_ACCESS:
		mexPrintf("SON File error: access denied, authorities on their way\n");
		break;
	case SON_READ_ONLY:
		mexPrintf("SON File write error: file is read-only\n");
		break;
	case SON_BAD_HANDLE:
		mexPrintf("SON File access error: handle is invalid\n");
		break;
	case SON_INVALID_DRIVE:
		mexPrintf("SON File operation failure: non-existent drive\n");
		break;
	case SON_MEMORY_ZAP:
		mexPrintf("SON System failure: son library memory allocation failed\n");
		break;
	case SON_OUT_OF_MEMORY:
		mexPrintf("SON System failure: couldn't allocate memory\n");
		break;
	case SON_NO_CHANNEL:
		mexPrintf("SON Channel error: requested channel out of valid range\n");
		break;
	case SON_CHANNEL_USED:
		mexPrintf("SON Channel error: requested channel already in use\n");
		break;
	case SON_CHANNEL_UNUSED:
		mexPrintf("SON Channel error: requested channel is not in use\n");
		break;
	case SON_PAST_EOF:
		mexPrintf("SON File error: attempt made to move past end of file\n");
		break;
	case SON_PAST_SOF:
		mexPrintf("SON File error: attempt made to move past start of file\n");
		break;
	case SON_WRONG_FILE:
		mexPrintf("SON File open error: file name is not a SON library file\n");
		break;
	case SON_NO_EXTRA:
		mexPrintf("SON File error: attempt to read/write past end of extra data, or no extra data\n");
		break;
	case SON_BAD_READ:
		mexPrintf("SON File read error: attempt to read from disk failed\n");
		break;
	case SON_BAD_WRITE:
		mexPrintf("SON File write error: attempt to write to disk failed\n");
		break;
	case SON_CORRUPT_FILE:
		mexPrintf("SON File error: file is corrupt\n");
		break;
	case SON_OUT_OF_HANDLES:
		mexPrintf("SON File error: SON library has run out of file handles because our library sucks\n");
		break;
	default:
		mexPrintf("SON Unknown error\n");
		break;
	}
} // End void son_error(int error_no)


void son_getMarkData(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
	if (nrhs == 6) {  /* Good, we have the correct # of input variables. */
		
		long transferredMarkers;		// Number of transferred markers.
		TpMarker markBuffer;			// Buffer which holds data read from file.
		long max;						// Maximum number of markers to be read.
		TSTime startTime, endTime;		// Start and end time
		int dimensions[2] = {1,1};
		char **fieldNames;
		int i;
		
		max = (long)mxGetScalar(prhs[3]);			// Get the maximum number of markers to return.
		startTime = (TSTime)mxGetScalar(prhs[4]);	// Get the start time of the search.
		endTime = (TSTime)mxGetScalar(prhs[5]);		// Get the end time of the search.
		
		// Initialize fieldNames, which will describe the the structure entries in the mxArray of
		// structures returned to Matlab.
		fieldNames = (char**)malloc(2 * sizeof(char*));
		fieldNames[0] = (char*)malloc(25 * sizeof(char));
		fieldNames[1] = (char*)malloc(25 * sizeof(char));
		strcpy(fieldNames[0], "mark");
		strcpy(fieldNames[1], "mvals");
		
		// Inititalize the buffer area where the Son32 function stores all event marker
		// data found.
		markBuffer = (TpMarker)malloc(max * sizeof(TMarker));
		
		// Get the number of markers.
		transferredMarkers = SONGetMarkData((short)mxGetScalar(prhs[1]), // handle
											(WORD)mxGetScalar(prhs[2]),	 // channel
											markBuffer,			 // holds marker data
											max,				 // maximum markers
											startTime,			 // start time
											endTime,			 // end time
											NULL                 /* filter */);
		
		// Check for error codes
		if (transferredMarkers == SON_BAD_READ || transferredMarkers == SON_NO_CHANNEL) {
			son_error(transferredMarkers);
			return;
		}

		free(markBuffer);
		markBuffer = (TpMarker)malloc(transferredMarkers * sizeof(TMarker));

		transferredMarkers = SONGetMarkData((short)mxGetScalar(prhs[1]), // handle
											(WORD)mxGetScalar(prhs[2]),	 // channel
											markBuffer,			 // holds marker data
											transferredMarkers,	 // maximum markers
											startTime,			 // start time
											endTime,			 // end time
											NULL                 /* filter */);

		// Check for error codes
		if (transferredMarkers < 0) {
			son_error(transferredMarkers);
			return;
		}
		
		// Create an mxArray to hold the number of retrieved markers and an mxStructure
		// array to hold the actual data.
		plhs[0] = mxCreateNumericArray(2, dimensions, mxDOUBLE_CLASS, mxREAL);
		mxGetPr(plhs[0])[0] = transferredMarkers;
		dimensions[1] = transferredMarkers;
		plhs[1] = mxCreateStructArray(2, dimensions, 2, fieldNames); 
		
		// Store all the data for Matlab.
		for (i = 0; i < transferredMarkers; i++) {
			mxArray* markData = mxCreateDoubleMatrix(1,1,mxREAL);
			mxArray* mvalData = mxCreateDoubleMatrix(1,1,mxREAL);
			mxGetPr(markData)[0] = markBuffer[i].mark;
			mxGetPr(mvalData)[0] = markBuffer[i].mvals[0];
			
			// Set the values.
			mxSetFieldByNumber(plhs[1], i, 0, markData);
			mxSetFieldByNumber(plhs[1], i, 1, mvalData);
		}
		
		// Free up memory.
		free(markBuffer);

	} // End if statement
	else { // Bad!!!
		mexErrMsgTxt("Incorrect number of inputs for SonGetMarkData.");
	}
}



/*********************************************************************************/
/*	Opens an existing SMR data file.											 */
/*********************************************************************************/
void son_openOldFile(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs) 
{
	if (nrhs == 3) {
		TpStr filename = (TpStr)mxArrayToString(prhs[1]);
		int openmode = (int)mxGetScalar(prhs[2]);
		short* output;
		int dims[2] = {1, 1};

		// Try to open the file.
		short handle = SONOpenOldFile(filename, openmode);

		// Print any errors.
		if (handle < 0) {
			son_error(handle);
		}

		// Create an mxArray for output of the handle.
		plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS,mxREAL);
		output = (short*)mxGetPr(plhs[0]);
		output[0] = handle;
	}
	else {
		mexErrMsgTxt("Incorrect Number of Inputs for SonOpenOldFile");
	}
}




/*********************************************************************************/
/*		Closes an open SMR file.					 */
/*********************************************************************************/
void son_closeFile(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {

	if (nrhs == 2) {
		// File handle
		short handle = (short)mxGetScalar(prhs[1]);
		
		// Try to close the SON file.
		short close_flag = SONCloseFile(handle);
		
		if (close_flag < 0) {
			son_error(close_flag);
		}
	}
	else
	{
		mexErrMsgTxt("Incorrect Number of Inputs for SONCloseFile");
	}
}



void mexFunction(
	int nlhs,              // Number of left hand side (output) arguments
	mxArray *plhs[],       // Array of left hand side arguments
	int nrhs,              // Number of right hand side (input) arguments
	const mxArray *prhs[]) // Array of right hand side arguments
{
	char *cedFunctionName;

	if(nrhs < 1) {
		mexErrMsgTxt("Need at least one input.");
	}

	if(mxIsChar(prhs[0]) == 1)
	{
		cedFunctionName = mxArrayToString(prhs[0]);

		if (strcmp("SonOpenOldFile", cedFunctionName) == STRING_MATCH)
		{
			son_openOldFile(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonCloseFile", cedFunctionName) == STRING_MATCH)
		{
			son_closeFile(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonChanMaxTime", cedFunctionName) == STRING_MATCH)
		{
			son_chanMaxTime(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonChanDivide", cedFunctionName) == STRING_MATCH)
		{
			son_chanDivide(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonGetADCData", cedFunctionName) == STRING_MATCH)
		{
			son_getADCData(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonGetNumADCPts", cedFunctionName) == STRING_MATCH)
		{
			son_getNumADCPts(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonGetusPerTime", cedFunctionName) == STRING_MATCH)
		{
			son_getusPerTime(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonGetAllADCData", cedFunctionName) == STRING_MATCH)
		{
			son_getAllADCData(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonGetMarkData", cedFunctionName) == STRING_MATCH) {
			son_getMarkData(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonChanKind", cedFunctionName) == STRING_MATCH) {
			son_chanKind(nlhs, plhs, nrhs, prhs);
		}
		else if (strcmp("SonGetEventData", cedFunctionName) == STRING_MATCH) {
			son_getEventData(nlhs, plhs, nrhs, prhs);
		}
		else {
			mexErrMsgTxt("Invalid function");
		}
	}
	else {
		mexErrMsgTxt("First argument must be a string.");
	}

} // End mexFunction