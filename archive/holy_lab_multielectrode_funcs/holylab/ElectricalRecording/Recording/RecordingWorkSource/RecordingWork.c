/*
	Acquire from National Instruments card, save to disk, and send
	graphical data back to MATLAB.
	Tim Holy, 4/20/99
 */

#include "mex.h"
#include "nidaq.h"
#include "nidaqcns.h"
#include "nidaqerr.h"
#include "iotypes.h"
#include "FileHeaders.cp"
#include <stdio.h>

extern "C" int16 Select_Signal(int16,uint32,uint32,uint32);

// Globals for driving the card
const int16 device = 1;
const int16 dbModeOff = 0;

void THMessage(char *s)
{
	DAQ_Clear( device );
	DAQ_DB_Config( device, dbModeOff );
	Timeout_Config( device, -1 /* disable timeout */ );

	mexErrMsgTxt(s);
}

//#define     kTrue             1
//#define     kFalse            0
char gtempstr[256];
#define     errcheck( err, msg )                                  \
               {                                                  \
                  if( err )                                       \
                  {  sprintf( gtempstr,"Error %d from %s.\n", err, msg );  \
					aih.nscans = loopCount*buffersize;			\
					aih.nscans.update(fpout);					\
					fclose(fpout);									\
               	 	THMessage(gtempstr);					  \
                  }                                               \
               }


int RecordingWork(int16 nchan,int16 *channels,float64 scanrate,uint32 nscans,uint32 buffersize,int trigger,int gain,
				int npix,const mxArray *axHP,const mxArray *textHP,const char *usrheader,const char *filename);

const int printing = 0;

// Gateway routine: do the argument parsing
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	const mxArray *axHP;
	const mxArray *textHP;

	int16 *channels;
	int nchan;
	double *dblP;
	float64 scanrate;
	uint32 nscans;
	uint32 buffersize;
	int npix;
	char *usrheader = "",*filename = "";
	int status;
	int mrows,ncols;
	int trigger;
	int gain;
	int i;
	int cbufflen;
	
	if (nrhs > 11 )
		mexErrMsgTxt( "Too many input arguments." );
	if (nrhs < 9 )
		mexErrMsgTxt( "Too few input arguments." );
	if (nlhs > 1 )
		mexErrMsgTxt( "Too many output arguments." );

	// Argument order: channels,scanrate,nscans,buffersize,trigger,gain,npix,axH,textH,usrhdr,filename
	// (usrhdr & filename are optional)	
	// Channels: noncomplex double vector, turn into array of ints
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
		mexErrMsgTxt("channels input must be non-complex double");
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	if (mrows != 1 && ncols != 1)
		mexErrMsgTxt("channels input must be a vector");
	nchan = mrows*ncols;
	channels = (int16*) mxMalloc(nchan*sizeof(int16));
	dblP = mxGetPr(prhs[0]);
	for (i = 0; i < nchan; i++) {
		channels[i] = dblP[i];
		if (channels[i] < 0 || channels[i] > 63)
			mexErrMsgTxt("channel numbers must be between 0 and 63");
	}
	// scanrate: noncomplex double scalar, turn into double
	mrows = mxGetM(prhs[1]);
	ncols = mxGetN(prhs[1]);
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mrows*ncols != 1)
		mexErrMsgTxt("scanrate input must be non-complex scalar double");
	scanrate = mxGetScalar(prhs[1]);
	// nscans: noncomplex double scalar, turn into uint32
	mrows = mxGetM(prhs[2]);
	ncols = mxGetN(prhs[2]);
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mrows*ncols != 1)
		mexErrMsgTxt("nscans input must be non-complex scalar double");
	if (mxGetScalar(prhs[2]) <= 0)
		mexErrMsgTxt("nscans must be positive");
	nscans = (uint32) mxGetScalar(prhs[2]);
	// buffersize: noncomplex double scalar, turn into uint32
	mrows = mxGetM(prhs[3]);
	ncols = mxGetN(prhs[3]);
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mrows*ncols != 1)
		mexErrMsgTxt("buffersize input must be non-complex scalar double");
	if (mxGetScalar(prhs[3]) <= 0)
		mexErrMsgTxt("buffersize must be positive");
	buffersize = (uint32) mxGetScalar(prhs[3]);
	// trigger: noncomplex double scalar, turn into int
	mrows = mxGetM(prhs[4]);
	ncols = mxGetN(prhs[4]);
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mrows*ncols != 1)
		mexErrMsgTxt("trigger input must be non-complex scalar double");
	trigger = (int) mxGetScalar(prhs[4]);
	if (trigger < 0 || trigger > 1)
		mexErrMsgTxt("trigger must be 0 or 1");
	// gain: noncomplex double scalar, turn into int
	mrows = mxGetM(prhs[5]);
	ncols = mxGetN(prhs[5]);
	if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mrows*ncols != 1)
		mexErrMsgTxt("gain input must be non-complex scalar double");
	gain = (int) mxGetScalar(prhs[5]);
	if (gain != -1 && gain != 1 && gain != 2 && gain != 5 && gain != 10 && gain != 20 && gain != 50 && gain != 100)
		mexErrMsgTxt("Valid gains are -1 (for 0.5), 1, 2, 5, 10, 20, 50, and 100");
	// npix: noncomplex double scalar, turn into int
	mrows = mxGetM(prhs[6]);
	ncols = mxGetN(prhs[6]);
	if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mrows*ncols != 1)
		mexErrMsgTxt("npix input must be non-complex scalar double");
	npix = (int) mxGetScalar(prhs[6]);
	// axH: noncomplex double vector
	mrows = mxGetM(prhs[7]);
	ncols = mxGetN(prhs[7]);
	if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]))
		mexErrMsgTxt("axH input must be non-complex double");
	if (mrows != 1 && ncols != 1)
		mexErrMsgTxt("axH must be a vector");
	axHP = prhs[7];
	// textH: noncomplex double scalar
	mrows = mxGetM(prhs[8]);
	ncols = mxGetN(prhs[8]);
	if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || mrows*ncols != 1)
		mexErrMsgTxt("textH input must be non-complex scalar double");
	textHP = prhs[8];

	if (nrhs > 9) {
		// usrheader: row vector of chars, convert to array of chars
		mrows = mxGetM(prhs[9]);
		ncols = mxGetN(prhs[9]);
/*
		if (!mxIsChar(prhs[9]) || mrows != 1)
			mexErrMsgTxt("usrheader input must be a row vector of characters");
*/
		cbufflen = mrows*ncols*sizeof(mxChar)+1;
		if (printing)
			mexPrintf("usrheader cbufflen = %d, [rows cols] = [%d %d]\n",cbufflen,mrows,ncols);
		usrheader = (char*) mxCalloc(cbufflen,sizeof(char));
		status = mxGetString(prhs[9],usrheader,cbufflen);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, usrheader string is truncated");
	}
	if (nrhs == 11) {
		// filename: row vector of chars, convert to array of chars
		mrows = mxGetM(prhs[10]);
		ncols = mxGetN(prhs[10]);
		if (!mxIsChar(prhs[10]) || mrows != 1)
			mexErrMsgTxt("filename input must be a row vector of characters");
		cbufflen = mrows*ncols*sizeof(mxChar)+1;
		filename = (char*) mxCalloc(cbufflen,sizeof(char));
		status = mxGetString(prhs[10],filename,cbufflen);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, filename string is truncated");
	}
	
	// Check the values, to make sure nothing bad will happen
	if (scanrate <= 0)
		mexErrMsgTxt("scanrate is <= 0");
	if (npix <= 0)
		 mexErrMsgTxt( "npix is <= 0" );

	if (printing) {
		mexPrintf("%d channels: ",nchan);
		for (i = 0; i < nchan; i++)
			mexPrintf(" %d",channels[i]);
		mexPrintf("\nscanrate %g, %d scans, buffersize %d, trigger %d\n",scanrate,nscans,buffersize,trigger);
		mexPrintf("gain %d, npix %d, usrheader %s,\nfilename %s\n",gain,npix,usrheader,filename);
	}
	
	status = RecordingWork(nchan,channels,scanrate,nscans,buffersize,trigger,gain,npix,axHP,textHP,usrheader,filename);

	mxFree(channels);
	if (nrhs > 9)
		mxFree(usrheader);
	if (nrhs == 11)
		mxFree(filename);

	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	dblP = mxGetPr(plhs[0]);
	dblP[0] = (double) status;
}

int RecordingWork(int16 nchan,int16 *channels,float64 scanrate,uint32 nscans,uint32 buffersize,int trigger,int gain,
				int npix,const mxArray *axHP,const mxArray *textHP,const char *usrheader,const char *filename)
{
	// Variables for interface with MATLAB
	mxArray *xcoordP;		// For graphing
	mxArray *patchHP;		// For graphing
	mxArray *scalarMxP;
	mxArray *lhs[1];		// For calling MATLAB interpreter
	mxArray *rhs[3];	// For calling MATLAB interpreter
	const mxArray *gIsRecordingMxP;	// Global flag, to check if user has pressed "Stop"
	mxArray *minmaxbuff;	// To hold decimated data
	mxArray *colorstring;
	mxArray *stringStrMxP;
	int gIsRecording;	// conversion of gIsRecordingMxP to integer value
	double *minmaxD,*maxcD,*mincD;		// Access to elements of mxArrays
	double *scalarMxD,*patchHD,*axHD,*dblD;

	// Variables for driving the card
	int16 *buffer,*halfbuffer;
	const int32   timeout = 180;  //ticks (10 secs)
	const int16 units = 0;	// scans/sec
    const int16 dbModeOn = 1;
    int16 sampleTimebase,scanTimebase;
    uint16 sampleInterval,scanInterval;
	int16 isDAQstopped;
    int16 isHalfReady;
    uint32 pointsTransferred,pointsWritten;
    const float timeBaseConvert[] = {50e-9,0,200e-9,0,1e-6,1e-5,1e-4,1e-3,1e-2};

	// Miscellaneous
	int saving;
	FilePtr fpout;
	AIHeader aih;
	int i,j,k;
	int status;
	uint32 decimate;
	int16 *dataP;
	float secsPerHalfBuff;
	char cTemp[256];
	int npixc;
	int loopCount = 0;


	// Set up acquisition buffer
	uint32 count = 2*buffersize*nchan;
	buffer = (int16*) mxMalloc(count*sizeof(int16));
	halfbuffer = (int16*) mxMalloc(count*sizeof(int16)/2);
	int numberOfHalfBuffersToRead = ceil((double)nscans/buffersize);
	uint32 remainder = nscans - buffersize*floor((double)nscans/buffersize);
	secsPerHalfBuff = (float)buffersize/scanrate;
	if (printing) {
		mexPrintf("secsPerHalfBuff = %g,buffersize %d,input scanrate %g\n",secsPerHalfBuff,buffersize,scanrate);
		mexPrintf("remainder = %d\n",remainder);
	}
	
	// Set up decimation
	decimate = buffersize/npix;
	if (npix*decimate != buffersize)
		mexErrMsgTxt("buffersize must be an integer multiple of npix");
	minmaxbuff = mxCreateDoubleMatrix(1,2*npix,mxREAL);	// Storage order: min(1),...,min(npix),max(npix),...,max(1)
	minmaxD = mxGetPr(minmaxbuff);
	
	// Set up handles for plotting
	// xcoord = [1:npix,npix:-1:1]
	xcoordP = mxCreateDoubleMatrix(1,2*npix,mxREAL);
	dblD = mxGetPr(xcoordP);
	for (i = 0; i < npix; i++)
		dblD[i] = dblD[2*npix-i-1] = i+1;
	patchHP = mxCreateDoubleMatrix(1,nchan,mxREAL);
	patchHD = mxGetPr(patchHP);
	mxSetName(patchHP,"patchH");			// We'll put it in the caller workspace, so give it a name
	colorstring = mxCreateString("b");
	scalarMxP = mxCreateDoubleMatrix(1,1,mxREAL);
	scalarMxD = mxGetPr(scalarMxP);
/*
	// So delete(patchH) doesn't give an error, have to fill with a dummy set of handles
	*scalarMxD = 0;
	rhs[0] = rhs[1] = scalarMxP;
	rhs[2] = colorstring;
	for (i = 1; i <= 64; i++) {
		mexCallMATLAB(1,lhs,3,rhs,"fill");
		patchHD[i-1] = mxGetScalar(lhs[0]);
		mxDestroyArray(lhs[0]);		// ?? Should I really do this??
	}
*/
	axHD = mxGetPr(axHP);
	stringStrMxP = mxCreateString("String");
	
	// Now get the card ready to go
	// set a timeout limit of 10 seconds (180 ticks) so that if something is
	//    wrong, the program won't hang on the DAQ_DB_Transfer call
/*
	if (!trigger) {
		status = Select_Signal(device,ND_PFI_0,ND_HIGH,ND_DONT_CARE);	// output trigger on PFI0
		status = DAQ_Clear(device);
	}	
*/
	status = Timeout_Config( device, timeout );
	errcheck( status, "Timeout_Config" );
	// Set up triggering: trigger = 0 means output the trigger signal (start automatically),
	//   while trigger = 1 means wait for trigger signal to do DAQ
	if (trigger) {
		status = Select_Signal(device,ND_PFI_0,ND_NONE,ND_DONT_CARE);	// disable output on PFI0
		errcheck(status,"Select_Signal1a");
		status = Select_Signal(device,ND_IN_START_TRIGGER,ND_PFI_0,ND_LOW_TO_HIGH);	// trigger DAQ off PFI0
		errcheck(status,"Select_Signal2a");
	}
	else {
		status = Select_Signal(device,ND_IN_START_TRIGGER,ND_AUTOMATIC,ND_DONT_CARE);	// trigger DAQ in software
		errcheck(status,"Select_Signal1b");
		status = Select_Signal(device,ND_PFI_0,ND_IN_START_TRIGGER,ND_LOW_TO_HIGH);	// output trigger on PFI0
		errcheck(status,"Select_Signal2b");
	}
	// convert the scan rate, which is in units of scans per second, to
	//    the appropriate timebase and sample interval values
	sampleTimebase = scanTimebase = -3;		// Use the fastest clock (50ns)
	float dt = timeBaseConvert[sampleTimebase+3];
	scanInterval = int(round(1.0/(dt*scanrate)));
	sampleInterval = int(scanInterval/nchan);
	scanrate = 1.0/(dt*scanInterval);		// Store the real scanrate, now
	if (printing)
		mexPrintf("sampleInterval %d, scanInterval %d, real scanrate %g\n",sampleInterval,scanInterval,scanrate);
	// Specify the channels and gains
	int16 *gains = (int16*) mxMalloc(nchan*sizeof(int16));
	for (i = 0; i < nchan; i++)
		gains[i] = gain;
	status = SCAN_Setup(device,nchan,channels,gains);
	errcheck(status,"SCAN_Setup");
	mxFree(gains);
	// enable double-buffered mode
	status = DAQ_DB_Config( device, dbModeOn );
	errcheck( status, "DAQ_DB_Config" );

	// Set up file, if saving
	// We needed to wait to do this until we got the real scanrate from DAQ_Rate
	saving = 0;
	if (strlen(filename) > 0) {
		saving = 1;
		
		fpout.open(filename,"wb");
		aih.nscans = nscans;
		aih.numCh = nchan;
		for (i = 0; i < nchan; i++)
			aih.channel.push_back(channels[i]);
		aih.scanrate = scanrate;
		float gainValue = gain;
		if (gainValue < 0)
			gainValue = 0.5;
		aih.scalemult = 5.0/(2048*gainValue);	// Assumes bipolar mode
		aih.scaleoff = 0;							// The card could also have non-ideal behavior, could try to correct
		// Get date & time
		mexCallMATLAB(1,lhs,0,NULL,"date");
		char *ctemp = (char*) mxCalloc(256,sizeof(char));
		status = mxGetString(lhs[0],ctemp,256);
		mxDestroyArray(lhs[0]);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, date string is truncated");
		aih.date = string(ctemp);
		mexCallMATLAB(1,lhs,0,NULL,"clock");
		dblD = mxGetPr(lhs[0]);
		sprintf(ctemp,"%d:%d:%4.2f",int(round(dblD[3])),int(round(dblD[4])),round(100*dblD[5])/100);
		mxDestroyArray(lhs[0]);
		aih.time = string(ctemp);
		mxFree(ctemp);
		aih.usrheader = string(usrheader);
		fpout << aih;
	}
    // Start scanning!
    status = SCAN_Start(device,buffer,count,sampleTimebase,sampleInterval,scanTimebase,scanInterval);
    errcheck(status,"SCAN_Start");
/*
	if (!trigger) {
		status = Select_Signal(device,ND_PFI_0,ND_HIGH,ND_DONT_CARE);	// output trigger on PFI0
		//errcheck(status,"Select_Signal3");
	}    
*/
	// Scanning loop: acquire half-buffers, save to disk (if applicable),
	// decimate data, and send back to MATLAB for plotting.
	// Must check to see if user has hit the "Stop" button, and, if so,
	// terminate early.
	//int loopCount = 0;
	gIsRecordingMxP = mexGetArrayPtr("gIsRecording","global");
	gIsRecording = mxGetScalar(gIsRecordingMxP);

	//numberOfHalfBuffersToRead = 2;
	while(( loopCount < numberOfHalfBuffersToRead ) && ( status == noError ) && gIsRecording) {
		status = DAQ_DB_HalfReady( device, &isHalfReady, &isDAQstopped );	// Is half buffer ready?

		if(( isHalfReady == 1 ) && ( status == noError )) {
			if (printing)
				mexPrintf("Transferring halfbuffer %d out of %d\n",loopCount+1,numberOfHalfBuffersToRead);
			// transfer the half-buffer of data into halfBuffer
			status = DAQ_DB_Transfer( device, halfbuffer, &pointsTransferred,&isDAQstopped );
			errcheck( status, "DAQ_DB_Transfer" );
/*
			// If this is our first time through, reset the triggering on the card
			if (loopCount == 0) {
				status = Select_Signal(device,ND_IN_START_TRIGGER,ND_AUTOMATIC,ND_DONT_CARE);	// trigger DAQ in software
				status = Select_Signal(device,ND_PFI_0,ND_NONE,ND_DONT_CARE);	// disconnect PFI0
			}
*/
			// Increment loop counter, and adjust pointsTransferred if this is our
			// last time through the loop
			loopCount++;
			npixc = npix;
			if (loopCount == numberOfHalfBuffersToRead) {
				pointsTransferred = remainder*nchan;
				npixc = ceil((double)remainder/decimate);
			}
			
			// Save the data, if appropriate
			if (saving) {
				pointsWritten = fwrite(halfbuffer,sizeof(int16),pointsTransferred,fpout);
				if (pointsWritten != pointsTransferred)
					errcheck(writeFileError,"File writing error");
			}
			
			// Graphical output:
			// First, clear the old data
			if (loopCount > 1) {
				if (mexPutArray(patchHP,"caller"))
					THMessage("Failed to put array patchH in caller workspace");
				sprintf(cTemp,"delete(patchH)");
				if (mexEvalString(cTemp))
					THMessage("deleting patchH failed!");
			}
			// Next, decimate the data for graphical output
			for (i = 0; i < nchan; i++) {
				dataP = halfbuffer+i;
				for (j = 0; j < npixc; j++) {
					mincD = minmaxD+j;
					maxcD = minmaxD+2*npix-j-1;
					*maxcD = *mincD = *dataP;
					for (k = 1; k < decimate; k++) {
						dataP += nchan;			// Move to next acquired point on this channel
						if (*dataP > *maxcD)
							*maxcD = *dataP;
						if (*dataP < *mincD)
							*mincD = *dataP;
					}
					dataP += nchan;
				}
				for (j = npixc; j < npix; j++) {	// Fill rest with zeros (happens only on final halfbuffer)
					mincD = minmaxD+j;
					maxcD = minmaxD+2*npix-j-1;
					*maxcD = *mincD = 0;
				}
				// Create plot data
				*scalarMxD = axHD[channels[i]];
				rhs[0] = scalarMxP;
				mexCallMATLAB(0,lhs,1,rhs,"axes");
				rhs[0] = xcoordP;
				rhs[1] = minmaxbuff;
				rhs[2] = colorstring;
				mexCallMATLAB(1,lhs,3,rhs,"fill");
				patchHD[i] = mxGetScalar(lhs[0]);
				mxDestroyArray(lhs[0]);	// ??
			}
			// Inform user of # of secs remaining
			sprintf(cTemp,"%d seconds remaining",(int)round(secsPerHalfBuff*(numberOfHalfBuffersToRead-loopCount)));
			rhs[0] = (mxArray*) textHP;
			rhs[1] = stringStrMxP;
			rhs[2] = mxCreateString(cTemp);
			mexCallMATLAB(0,lhs,3,rhs,"set");
			mxDestroyArray(rhs[2]);
			// Force plotting and check for user input (Stop button)
			mexCallMATLAB(0,lhs,0,NULL,"drawnow");
			gIsRecording = mxGetScalar(gIsRecordingMxP);
		}	// if (isHalfReady == 1)...
		else {
			errcheck( status, "DAQ_DB_HalfReady" );
			mexCallMATLAB(0,lhs,0,NULL,"drawnow");	// Give user opportunity to cancel, esp. important for TrigIn
		}
	}
	// If user hit the "Stop" button, figure out how many scans were taken and
	// update disk file accordingly
	if (saving && !gIsRecording) {
		aih.nscans = loopCount*buffersize;
		aih.nscans.update(fpout);
	}
	if (status != noError)
		THMessage("Aborting acquisition because of error.");

	// Clean up card; deliberately don't check for errors
/*
	if (!trigger) {
		status = Select_Signal(device,ND_PFI_0,ND_IN_START_TRIGGER,ND_LOW_TO_HIGH);	// output trigger on PFI0
		//errcheck(status,"Select_Signal4");
	}    
*/
	status = DAQ_Clear( device );
	status = DAQ_DB_Config( device, dbModeOff );
	status = Timeout_Config( device, -1); // disable timeout
/*
	if (!trigger) {
		status = Select_Signal(device,ND_PFI_0,ND_LOW,ND_DONT_CARE);	// output trigger on PFI0
		status = DAQ_Clear(device);
	}	
*/
	//int16 deviceNumberCode;
	//status = Init_DA_Brds(device,&deviceNumberCode);

	mxDestroyArray(minmaxbuff);
	mxDestroyArray(xcoordP);
	mxDestroyArray(patchHP);
	mxDestroyArray(colorstring);
	mxDestroyArray(scalarMxP);
	mxDestroyArray(stringStrMxP);
	mxFree(buffer);
	mxFree(halfbuffer);

	if (!gIsRecording)
		return 1;			// If user hit stop

	return 0;				// Normal return, everything went well
}