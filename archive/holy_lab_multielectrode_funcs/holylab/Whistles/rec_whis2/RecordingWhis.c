/*
	Acquire from National Instruments card, save to disk, and send
	graphical data back to MATLAB.
	Also computes sonogram on the fly.
	Tim Holy, 4/20/99
	Modified 4/4/00 to drive the Sony camcorder.
	
	Normalization of output sonogram: the mean across frequencies and time is
	the mean total power. Units are in volts^2.
 */


#include <MacHeaders.h>
#include "mex.h"
#include "nidaq.h"
#include "nidaqcns.h"
#include "nidaqerr.h"
#include "iotypes.h"
#include "FileHeaders.cp"
#include <stdio.h>
#include <string>
#include "SonySerial.h"

using namespace std;

extern "C" int16 Select_Signal(int16,uint32,uint32,uint32);
extern "C" void spctrmD(float d[], double p[], int m, int k, int ovrlap);


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


mxArray* RecordingWhis(int16 channel,float64 scanrate,uint32 nfreq,int16 navg,int16 nperblock,uint32 nblocks,int gain,
				int npix,mxArray *imHP,const mxArray *textHP,int16 video,const char *usrheader,const char *filename);
void MyTextOut(char *s);
void MyWait();

const int printing = 0;

// Gateway routine: do the argument parsing
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	if (printing) mexPrintf("Version 3\n");
	const mxArray *imHP;
	const mxArray *textHP;
	mxArray *sngP;

	int16 channel;
	uint32 nfreq;
	int16 navg;
	int16 nperblock;
	uint32 nblocks;
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
	int16 video;
	
	if (nrhs > 13 )
		mexErrMsgTxt( "Too many input arguments." );
	if (nrhs < 11 )
		mexErrMsgTxt( "Too few input arguments." );
	if (nlhs > 1 )
		mexErrMsgTxt( "Too many output arguments." );

	// Argument order: channels,scanrate,nfreq,navg,nperblock,nblocks,gain,npix,imH,textH,usrhdr,filename
	// (usrhdr & filename are optional)	
	// Channels: noncomplex double vector, turn into array of ints
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mrows*ncols != 1)
		mexErrMsgTxt("channel input must be non-complex scalar double");
	channel = (int16) mxGetScalar(prhs[0]);
	// scanrate: noncomplex double scalar, turn into double
	mrows = mxGetM(prhs[1]);
	ncols = mxGetN(prhs[1]);
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mrows*ncols != 1)
		mexErrMsgTxt("scanrate input must be non-complex scalar double");
	scanrate = mxGetScalar(prhs[1]);
	// nfreq: noncomplex double scalar, turn into uint32
	mrows = mxGetM(prhs[2]);
	ncols = mxGetN(prhs[2]);
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mrows*ncols != 1)
		mexErrMsgTxt("nfreq input must be non-complex scalar double");
	if (mxGetScalar(prhs[2]) <= 0)
		mexErrMsgTxt("nfreq must be positive");
	nfreq = (uint32) mxGetScalar(prhs[2]);
	// navg: noncomplex double scalar, turn into int16
	mrows = mxGetM(prhs[3]);
	ncols = mxGetN(prhs[3]);
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mrows*ncols != 1)
		mexErrMsgTxt("navg input must be non-complex scalar double");
	if (mxGetScalar(prhs[3]) <= 0)
		mexErrMsgTxt("navg must be positive");
	navg = (int16) mxGetScalar(prhs[3]);
	if (navg % 2 != 0)
		mexErrMsgTxt("navg must be an even integer");
	// trigger: noncomplex double scalar, turn into int16
	mrows = mxGetM(prhs[4]);
	ncols = mxGetN(prhs[4]);
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mrows*ncols != 1)
		mexErrMsgTxt("nperblock input must be non-complex scalar double");
	nperblock = (int16) mxGetScalar(prhs[4]);
	// nblocks: noncomplex double scalar, turn into uint32
	mrows = mxGetM(prhs[5]);
	ncols = mxGetN(prhs[5]);
	if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mrows*ncols != 1)
		mexErrMsgTxt("nblocks input must be non-complex scalar double");
	if (mxGetScalar(prhs[5]) <= 0)
		mexErrMsgTxt("blocks must be positive");
	nblocks = (uint32) mxGetScalar(prhs[5]);
	// gain: noncomplex double scalar, turn into int
	mrows = mxGetM(prhs[6]);
	ncols = mxGetN(prhs[6]);
	if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mrows*ncols != 1)
		mexErrMsgTxt("gain input must be non-complex scalar double");
	gain = (int) mxGetScalar(prhs[6]);
	if (gain != -1 && gain != 1 && gain != 2 && gain != 5 && gain != 10 && gain != 20 && gain != 50 && gain != 100)
		mexErrMsgTxt("Valid gains are -1 (for 0.5), 1, 2, 5, 10, 20, 50, and 100");
	// npix: noncomplex double scalar, turn into int
	mrows = mxGetM(prhs[7]);
	ncols = mxGetN(prhs[7]);
	if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) || mrows*ncols != 1)
		mexErrMsgTxt("npix input must be non-complex scalar double");
	npix = (int) mxGetScalar(prhs[7]);
	// imH: noncomplex double vector
	mrows = mxGetM(prhs[8]);
	ncols = mxGetN(prhs[8]);
	if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]))
		mexErrMsgTxt("imH input must be non-complex double");
	if (mrows != 1 && ncols != 1)
		mexErrMsgTxt("imH must be a vector");
	imHP = prhs[8];
	// textH: noncomplex double scalar
	mrows = mxGetM(prhs[9]);
	ncols = mxGetN(prhs[9]);
	if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || mrows*ncols != 1)
		mexErrMsgTxt("textH input must be non-complex scalar double");
	textHP = prhs[9];
	// video: noncomplex double scalar
	mrows = mxGetM(prhs[10]);
	ncols = mxGetN(prhs[10]);
	if (!mxIsDouble(prhs[10]) || mxIsComplex(prhs[10]) || mrows*ncols != 1)
		mexErrMsgTxt("videoFlag input must be non-complex scalar double");
	video = (int) mxGetScalar(prhs[10]);

	if (nrhs > 11) {
		// usrheader: row vector of chars, convert to array of chars
		mrows = mxGetM(prhs[11]);
		ncols = mxGetN(prhs[11]);
/*
		if (!mxIsChar(prhs[11]) || mrows != 1)
			mexErrMsgTxt("usrheader input must be a row vector of characters");
*/
		cbufflen = mrows*ncols*sizeof(mxChar)+1;
		if (printing)
			mexPrintf("usrheader cbufflen = %d, [rows cols] = [%d %d]\n",cbufflen,mrows,ncols);
		usrheader = (char*) mxCalloc(cbufflen,sizeof(char));
		status = mxGetString(prhs[11],usrheader,cbufflen);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, usrheader string is truncated");
	}
	if (nrhs == 13) {
		// filename: row vector of chars, convert to array of chars
		mrows = mxGetM(prhs[12]);
		ncols = mxGetN(prhs[12]);
		if (!mxIsChar(prhs[12]) || mrows != 1)
			mexErrMsgTxt("filename input must be a row vector of characters");
		cbufflen = mrows*ncols*sizeof(mxChar)+1;
		filename = (char*) mxCalloc(cbufflen,sizeof(char));
		status = mxGetString(prhs[12],filename,cbufflen);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, filename string is truncated");
	}
	
	// Check the values, to make sure nothing bad will happen
	if (scanrate <= 0)
		mexErrMsgTxt("scanrate is <= 0");
	if (npix <= 0)
		 mexErrMsgTxt( "npix is <= 0" );

	if (printing) {
		mexPrintf("channel: %d",channel);
		mexPrintf("\nscanrate %g, %d scans, buffersize %d\n",scanrate,2*nfreq*navg*nperblock*nblocks,2*nfreq*navg*nperblock);
		mexPrintf("gain %d, npix %d, usrheader %s,\nfilename %s\n",gain,npix,usrheader,filename);
	}
	
	sngP = RecordingWhis(channel,scanrate,nfreq,navg,nperblock,nblocks,gain,npix,(mxArray*) imHP,textHP,video,usrheader,filename);

	if (nrhs > 11)
		mxFree(usrheader);
	if (nrhs == 13)
		mxFree(filename);

	plhs[0] = sngP;
}

void MyTextOut(char *s)
{
	mexPrintf("%s\n",s);
	mexEvalString("drawnow");
}

void MyWait()
{
	//SystemTask();
}

mxArray* RecordingWhis(int16 channel,float64 scanrate,uint32 nfreq,int16 navg,int16 nperblock,uint32 nblocks,int gain,
				int npix,mxArray *imHP,const mxArray *textHP,int16 video,const char *usrheader,const char *filename)
{
	// Variables for interface with MATLAB
	//mxArray *xcoordP;		// For graphing
	//mxArray *patchHP;		// For graphing
	mxArray *scalarMxP;
	mxArray *sngP,*logsngP;
	mxArray *lhs[1];		// For calling MATLAB interpreter
	mxArray *rhs[3];	// For calling MATLAB interpreter
	const mxArray *gIsRecordingMxP;	// Global flag, to check if user has pressed "Stop"
	//mxArray *minmaxbuff;	// To hold decimated data
	//mxArray *colorstring;
	mxArray *stringStrMxP,*string2StrMxP;
	int gIsRecording;	// conversion of gIsRecordingMxP to integer value
	double *minmaxD,*maxcD,*mincD;		// Access to elements of mxArrays
	double *scalarMxD,*patchHD,*dblD,*sngD,*logsngD;

	// Variables for driving the card
	uint32 buffersize;
	int16 *buffer,*halfbuffer;
	float *fourbuf;
	const int32   timeout = 180;  //ticks (10 secs)
	const int16 units = 0;	// scans/sec
    const int16 dbModeOn = 1;
    int16 sampleTimebase,scanTimebase;
    uint16 sampleInterval,scanInterval;
	int16 isDAQstopped;
    int16 isHalfReady;
    uint32 pointsTransferred,pointsWritten;
    const float timeBaseConvert[] = {50e-9,0,200e-9,0,1e-6,1e-5,1e-4,1e-3,1e-2};
	int nchan = 1;
	int16 *channels = &channel;
	uint32 nscans = 2*nfreq*navg*nperblock*nblocks;

	// Controlling the camcorder
	long finishTime,tc0;
	VTRTimeCode tc1,tc2;
	long filepos;
	SonySerial serport(0);	// Use modem port
	if (printing)
		mexPrintf("Got past serport initialization\n");
		
	// Miscellaneous
	int saving;
	FilePtr fpout;
	AIVidHeader aih;
	int i,j,k;
	int status;
	uint32 decimate;
	int16 *dataP;
	float secsPerHalfBuff;
	char cTemp[256];
	int npixc;
	int loopCount = 0;
	float scalemult,scaleoff;

	// First do the video (if required) so can exit quickly if there is trouble
	if (video) {
		if (!serport.Stop())	// Get started recording immediately, so tape can get moving
			mexErrMsgTxt("Camcorder is not responding");
		//serport.ReadTimeCode(tc1);	// Record the timecode at current spot, to be able to check that tape is really moving
	}

	// Set up acquisition buffer
	buffersize = 2*nfreq*navg*nperblock;
	uint32 count = 2*buffersize;
	buffer = (int16*) mxMalloc(count*sizeof(int16));
	halfbuffer = (int16*) mxMalloc(count*sizeof(int16)/2);
	int numberOfHalfBuffersToRead = nblocks;
	uint32 remainder = 0;
	secsPerHalfBuff = (float)buffersize/scanrate;
	if (printing) {
		mexPrintf("secsPerHalfBuff = %g,buffersize %d,input scanrate %g\n",secsPerHalfBuff,buffersize,scanrate);
		mexPrintf("remainder = %d\n",remainder);
	}
	
	// Set up stuff for computing sonogram
	fourbuf = (float*) mxMalloc(2*nfreq*navg*sizeof(float));
	sngP = mxCreateDoubleMatrix(nfreq,nperblock*nblocks,mxREAL);
	sngD = mxGetPr(sngP);
	//logsngP = mxCreateDoubleMatrix(nfreq,nperblock*nblocks,mxREAL);
	
/*
	// Set up decimation
	decimate = buffersize/npix;
	if (npix*decimate != buffersize)
		mexErrMsgTxt("buffersize must be an integer multiple of npix");
	minmaxbuff = mxCreateDoubleMatrix(1,2*npix,mxREAL);	// Storage order: min(1),...,min(npix),max(npix),...,max(1)
	minmaxD = mxGetPr(minmaxbuff);
*/
	
	// Set up plotting
	stringStrMxP = mxCreateString("CData");
	rhs[0] = imHP;
	rhs[1] = stringStrMxP;
	mexCallMATLAB(1,lhs,2,rhs,"get");	// MATLAB command: logsng = get(imageh,'CData');
	logsngP = lhs[0];
	logsngD = mxGetPr(logsngP);
			//rhs[0] = imHP;
			//rhs[1] = stringStrMxP;
			//rhs[2] = logsngP;
			//mexCallMATLAB(0,lhs,3,rhs,"set");	// MATLAB command: set(imageh,'CData',logsng)
/*
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
	imHD = mxGetPr(imHP);
*/
	string2StrMxP = mxCreateString("String");
	
	// Now get the card ready to go
	// set a timeout limit of 10 seconds (180 ticks) so that if something is
	//    wrong, the program won't hang on the DAQ_DB_Transfer call
	status = Timeout_Config( device, timeout );
	errcheck( status, "Timeout_Config" );
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

	float gainValue = gain;
	if (gainValue < 0)
		gainValue = 0.5;
	scalemult = 5.0/(2048*gainValue);	// Assumes bipolar mode
	scaleoff = 0;		// The card could also have non-ideal behavior, could try to correct

/*
	// Check & see if the tape is ready (if recording video)
	if (video) {
		long finishTime;
		int diff,err;
		
		diff = 0;
		err = 0;
		finishTime = TickCount() + 5*60;
		if (printing)
			mexPrintf("Base: %d %d %d %d\n",tc1[0],tc1[1],tc1[2],tc1[3]);
		
		while (diff == 0 && TickCount() < finishTime && !err) {
			err = !serport.ReadTimeCode(tc2);
			diff = tc2[0]-tc1[0];
			if (printing)
				mexPrintf("%d %d %d %d\n",tc2[0],tc2[1],tc2[2],tc2[3]);
			MyWait();
		}
		if (TickCount() >= finishTime || err) {
			status = DAQ_Clear( device );
			status = DAQ_DB_Config( device, dbModeOff );
			status = Timeout_Config( device, -1); // disable timeout
			mxDestroyArray(minmaxbuff);
			mxDestroyArray(xcoordP);
			mxDestroyArray(patchHP);
			mxDestroyArray(colorstring);
			mxDestroyArray(scalarMxP);
			mxDestroyArray(stringStrMxP);
			mxFree(buffer);
			mxFree(halfbuffer);
			if (err)
				mexErrMsgTxt("Timecode reading error");
			else
				mexErrMsgTxt("Tape never got rolling");
		}
		for (i = 0; i < 4; i++)
			aih.tc[i] = tc2[i];
		for (i = 0; i < 4; i++)
			aih.tc[i] = tc1[i];
		
		if (printing)
			mexPrintf("Written timecode: %d %d %d %d\n",aih.tc[0],aih.tc[1],aih.tc[2],aih.tc[3]);
	}
	else {
		for (i = 0; i < 4; i++)
			aih.tc[i] = 0;
	}	// if (video)
*/

	// Set up file, if saving
	// We needed to wait to do this until we got the real scanrate
	//	and the video timecode
	saving = 0;
	if (strlen(filename) > 0) {
		saving = 1;
		
		fpout.open(filename,"wb");
		aih.nscans = nscans;
		aih.numCh = nchan;
		for (i = 0; i < nchan; i++)
			aih.channel.push_back(channels[i]);
		aih.scanrate = scanrate;
		aih.scalemult = scalemult;
		aih.scaleoff = scaleoff;
		// Get date & time
		mexCallMATLAB(1,lhs,0,NULL,"date");
		char *ctemp = (char*) mxCalloc(256,sizeof(char));
		status = mxGetString(lhs[0],ctemp,256);
		mxDestroyArray(lhs[0]);
		if (status != 0)
			mexWarnMsgTxt("Not enough space, date string is truncated");
		aih.date = LenString(ctemp);
		mexCallMATLAB(1,lhs,0,NULL,"clock");
		dblD = mxGetPr(lhs[0]);
		sprintf(ctemp,"%d:%d:%4.2f",int(round(dblD[3])),int(round(dblD[4])),round(100*dblD[5])/100);
		mxDestroyArray(lhs[0]);
		aih.time = LenString(ctemp);
		mxFree(ctemp);
		aih.usrheader = LenString(usrheader);
		fpout << aih;
		filepos = fpout.tell();
	}

	// Start the tape recording
	if (video) {
		if (!serport.ToggleRecord())
			MyTextOut("Camcorder not responding as expected. Continuing anyway");
		//serport.ReadTimeCode(AIVidHeader::VTRTimeCode(aih.tc));
		serport.ReadTimeCode(aih.tc);
		mexPrintf("Saved time code: %d %d %d %d\n",aih.tc[0],aih.tc[1],aih.tc[2],aih.tc[3]);
	}
	
    // Start scanning!
    status = SCAN_Start(device,buffer,count,sampleTimebase,sampleInterval,scanTimebase,scanInterval);
    errcheck(status,"SCAN_Start");
	// Scanning loop: acquire half-buffers, save to disk (if applicable),
	// decimate data, and send back to MATLAB for plotting.
	// Must check to see if user has hit the "Stop" button, and, if so,
	// terminate early.
	//int loopCount = 0;
	gIsRecordingMxP = mexGetArrayPtr("gIsRecording","global");
	gIsRecording = mxGetScalar(gIsRecordingMxP);

	//numberOfHalfBuffersToRead = 1;
	while(( loopCount < numberOfHalfBuffersToRead ) && ( status == noError ) && gIsRecording) {
		status = DAQ_DB_HalfReady( device, &isHalfReady, &isDAQstopped );	// Is half buffer ready?

		if(( isHalfReady == 1 ) && ( status == noError )) {
			if (printing)
				mexPrintf("Transferring halfbuffer %d out of %d\n",loopCount+1,numberOfHalfBuffersToRead);
			// transfer the half-buffer of data into halfBuffer
			status = DAQ_DB_Transfer( device, halfbuffer, &pointsTransferred,&isDAQstopped );
			errcheck( status, "DAQ_DB_Transfer" );
			// Increment loop counter, and adjust pointsTransferred if this is our
			// last time through the loop
			loopCount++;
			npixc = npix;
			//if (loopCount == numberOfHalfBuffersToRead) {
			//	pointsTransferred = remainder*nchan;
			//	npixc = ceil((double)remainder/decimate);
			//}
			
			// Save the data, if appropriate
			if (saving) {
				pointsWritten = fwrite(halfbuffer,sizeof(int16),pointsTransferred,fpout);
				if (pointsWritten != pointsTransferred)
					errcheck(writeFileError,"File writing error");
			}
			// Compute sonogram for this stretch of data
			int16 *hbi = halfbuffer;
			uint32 nptsfour = 2*nfreq*navg;
			float *fbi;
			for (i = 0; i < nperblock; i++) {
				fbi = fourbuf;
				for (j = 0; j < nptsfour; j++,fbi++,hbi++)			// Copy the data into float buffer,
					*fbi = *hbi * scalemult + scaleoff;	// converting to volts
				spctrmD(fourbuf-1,sngD-1,nfreq,navg/2,0);		// PSD with no overlap (for speed)
				for (j = 0; j < nfreq; j++,sngD++,logsngD++) {
					*sngD *= nfreq;							// Now the mean across freqs is the mean total power
															// ("sum squared amplitude" normalization in NR parlance)
					*logsngD = log10(*sngD);
				}
			}
			// Graphical output of sonogram
			//rhs[0] = stringStrMxP;
			//mexCallMATLAB(0,lhs,1,rhs,"disp");
			rhs[0] = imHP;
			rhs[1] = stringStrMxP;
			rhs[2] = logsngP;
			//mexPrintf("imHP %g\n",*mxGetPr(imHP));
			mexCallMATLAB(0,lhs,3,rhs,"set");	// MATLAB command: set(imageh,'CData',logsng)
			
/*			
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
				*scalarMxD = imHD[channels[i]];
				rhs[0] = scalarMxP;
				mexCallMATLAB(0,lhs,1,rhs,"axes");
				rhs[0] = xcoordP;
				rhs[1] = minmaxbuff;
				rhs[2] = colorstring;
				mexCallMATLAB(1,lhs,3,rhs,"fill");
				patchHD[i] = mxGetScalar(lhs[0]);
				mxDestroyArray(lhs[0]);	// ??
			}
*/
			// Inform user of # of secs remaining
			sprintf(cTemp,"%d seconds remaining",(int)round(secsPerHalfBuff*(numberOfHalfBuffersToRead-loopCount)));
			rhs[0] = (mxArray*) textHP;
			rhs[1] = string2StrMxP;
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
	// Set the size of the sng matrix to the actual size acquired
	mxSetN(sngP,loopCount*nperblock);
	
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

	if (video) {
		serport.ToggleRecord();		// Stop camcorder
//		aih.tc.update(fpout);		// Write the correct time code
		fpout.goTo(filepos-4);
		for (i = 0; i < 4; i++)
			fpout << aih.tc[i];
	}

/*
	mxDestroyArray(minmaxbuff);
	mxDestroyArray(xcoordP);
	mxDestroyArray(patchHP);
	mxDestroyArray(colorstring);
	mxDestroyArray(scalarMxP);
*/
	mxDestroyArray(stringStrMxP);
	mxDestroyArray(string2StrMxP);
	mxDestroyArray(logsngP);
	mxFree(buffer);
	mxFree(halfbuffer);
	mxFree(fourbuf);

	return sngP;
}