#include "mex.h"
#include "SonySerial.h"

void MyTextOut(char *s);
void MyWait();
void ClosePort(void);

void MyTextOut(char *s)
{
	mexPrintf("%s\n",s);
}

void MyWait()
{
}

enum SonyCommands { kStop = 0,
					kPause,
					kPlay,
					kRecord,
					kReadTimeCode,
					kCue,
					kCueReady,
					kScanForward,
					kScanBackward,
					kSlowForward,
					kSlowBackward,
					kRewind,
					kFastForward};

SonySerial ssport(0);

void ClosePort(void)
{
	ssport.SonySerial::~SonySerial();	// Call the destructor
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int mrows,ncols;
	int command;
	VTRTimeCode tc;
	double *dblP;
	static bool firstpass = false;
	bool status;
	int i;

	if (nrhs < 1 || nrhs > 2)
		mexErrMsgTxt("One or two inputs are required");
	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("One or two outputs are required");
	
	mexAtExit(ClosePort);
	if (firstpass) {
		while (!ssport.ReadTimeCode(tc))
			mexPrintf("First time in: warm up by reading time code");
		firstpass = false;
	}
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mrows*ncols != 1)
		mexErrMsgTxt("command must be non-complex scalar double");
	command = (int) mxGetScalar(prhs[0]);
	
	status = false;
	if (command == kCue && nrhs == 2) {
		mrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mrows*ncols != 4)
			mexErrMsgTxt("timecode must be non-complex double with 4 elements");
		dblP = mxGetPr(prhs[1]);
		for (i = 0; i < 4; i++)
			tc[i] = (byte) dblP[i];
		// Execute the command right now
		status = ssport.Cue(tc);
/*
		while (!ssport.CueReady())
			mexEvalString("pause(1)");
*/
	}
	else if (command == kStop)
		status = ssport.Stop();
	else if (command == kPause)
		status = ssport.Pause();
	else if (command == kPlay)
		status = ssport.Play();
	else if (command == kReadTimeCode) {
		if (nlhs < 2)
			mexErrMsgTxt("Must have 2 outputs to read the time code");
		status = ssport.ReadTimeCode(tc);
		plhs[1] = mxCreateDoubleMatrix(1,4,mxREAL);
		dblP = mxGetPr(plhs[1]);
		for (i = 0; i < 4; i++)
			dblP[i] = (double) tc[i];
	}
	else if (command == kCueReady)
		status = ssport.CueReady();
	else if (command == kScanForward)
		status = ssport.ScanForward();
	else if (command == kScanBackward)
		status = ssport.ScanBackward();
	else if (command == kSlowForward)
		status = ssport.SlowForward();
	else if (command == kSlowBackward)
		status = ssport.SlowBackward();
	else if (command == kRewind)
		status = ssport.Rewind();
	else if (command == kFastForward)
		status = ssport.FastForward();
	else
		mexErrMsgTxt("Unrecognized command");

	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	dblP = mxGetPr(plhs[0]);
	*dblP = status;
}