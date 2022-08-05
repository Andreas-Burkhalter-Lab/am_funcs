/*
	Computes sonogram from a disk file. Optionally, graph in real time.
	Tim Holy, 4/27/00
	
	Normalization of output sonogram: the mean across frequencies and time is
	the mean total power. Units are in volts^2.
 */


#include <MacHeaders.h>
#include "mex.h"
#include "iotypes.h"
#include "FileHeaders.cp"
#include <stdio.h>
#include <string>

using namespace std;

extern "C" void spctrmD(float d[], double p[], int m, int k, int ovrlap);

char gtempstr[256];

mxArray* sonogram(const char *filename,uint32 nfreq,int16 navg,int16 plotint,mxArray *imHP);
void MyTextOut(char *s);
void MyWait();

void THMessage(char *s)
{
	MyTextOut(s);
}

const int printing = 1;

// Gateway routine: do the argument parsing
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	if (printing) mexPrintf("Version 1\n");
	const mxArray *imHP;
	mxArray *sngP;

	uint32 nfreq;
	int16 navg;
	int16 plotint;
	double *dblP;
	char *filename;
	int status;
	int mrows,ncols;
	int i;
	int cbufflen;
	
	if (nrhs > 5 )
		mexErrMsgTxt( "Too many input arguments." );
	if (nrhs < 3 )
		mexErrMsgTxt( "Too few input arguments." );
	if (nlhs > 1 )
		mexErrMsgTxt( "Too many output arguments." );

	filename = "";
	
	// Argument order: filename,nfreq,navg,plotint,imH
	// filename: row vector of chars, convert to array of chars
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	if (!mxIsChar(prhs[0]) || mrows != 1)
		mexErrMsgTxt("filename input must be a row vector of characters");
	cbufflen = mrows*ncols*sizeof(mxChar)+1;
	filename = (char*) mxCalloc(cbufflen,sizeof(char));
	status = mxGetString(prhs[0],filename,cbufflen);
	if (status != 0)
		mexWarnMsgTxt("Not enough space, filename string is truncated");
	// nfreq: noncomplex double scalar, turn into uint32
	mrows = mxGetM(prhs[1]);
	ncols = mxGetN(prhs[1]);
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mrows*ncols != 1)
		mexErrMsgTxt("nfreq input must be non-complex scalar double");
	if (mxGetScalar(prhs[1]) <= 0)
		mexErrMsgTxt("nfreq must be positive");
	nfreq = (uint32) mxGetScalar(prhs[1]);
	// navg: noncomplex double scalar, turn into int16
	mrows = mxGetM(prhs[2]);
	ncols = mxGetN(prhs[2]);
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mrows*ncols != 1)
		mexErrMsgTxt("navg input must be non-complex scalar double");
	if (mxGetScalar(prhs[2]) <= 0)
		mexErrMsgTxt("navg must be positive");
	navg = (int16) mxGetScalar(prhs[2]);
	if (navg % 2 != 0)
		mexErrMsgTxt("navg must be an even integer");
	if (nrhs > 3) {
		// plotint: noncomplex double scalar, turn into int16
		mrows = mxGetM(prhs[3]);
		ncols = mxGetN(prhs[3]);
		if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mrows*ncols != 1)
			mexErrMsgTxt("plotint input must be non-complex scalar double");
		plotint = (int16) mxGetScalar(prhs[3]);
	} else
		plotint = 0;
	if (nrhs > 4) {
		// imH: noncomplex double vector
		mrows = mxGetM(prhs[4]);
		ncols = mxGetN(prhs[4]);
		if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
			mexErrMsgTxt("imH input must be non-complex double");
		if (mrows != 1 && ncols != 1)
			mexErrMsgTxt("imH must be a vector");
		imHP = prhs[4];
	} else
		imHP = 0;
	
	sngP = sonogram(filename,nfreq,navg,plotint,(mxArray*) imHP);

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

mxArray* sonogram(const char *filename,uint32 nfreq,int16 navg,int16 plotint,mxArray *imHP)
{
	// Variables for interface with MATLAB
	mxArray *scalarMxP;
	mxArray *sngP,*logsngP;
	mxArray *lhs[1];		// For calling MATLAB interpreter
	mxArray *rhs[3];	// For calling MATLAB interpreter
	mxArray *stringStrMxP,*string2StrMxP;
	double *sngD,*logsngD;
		
	// Miscellaneous
	int npoints,ntimes;
	float *fourbuf,*fbi;
	int16 *diskbuf,*hbi;
	FilePtr fpin;
	AIVidHeader aih;
	int i,j,k;
	int status;
	int16 *dataP;
	char cTemp[256];
	float scalemult,scaleoff;
	int loopCount;
	size_t pointsRead;
	
	// Read file size info
	fpin.open(filename,"rb");
	fpin >> aih;
	
	// Compute number of time points in sonogram
	npoints = 2*nfreq*navg;
	ntimes = aih.nscans/npoints;
	scalemult = aih.scalemult;
	scaleoff = aih.scaleoff;

	// Set up plotting
	if (imHP) {
		stringStrMxP = mxCreateString("CData");
		rhs[0] = imHP;
		rhs[1] = stringStrMxP;
		mexCallMATLAB(1,lhs,2,rhs,"get");	// MATLAB command: logsng = get(imageh,'CData');
		logsngP = lhs[0];
		k = mxGetN(logsngP);
		if (k != ntimes)
			mexErrMsgTxt("Size of image does not match size of sonogram");
		logsngD = mxGetPr(logsngP);
	}

	// Set up stuff for computing sonogram
	diskbuf = (int16*) mxMalloc(npoints*sizeof(int16));
	fourbuf = (float*) mxMalloc(npoints*sizeof(float));
	sngP = mxCreateDoubleMatrix(nfreq,ntimes,mxREAL);
	sngD = mxGetPr(sngP);
	
	// Get going!
	for (loopCount = 0; loopCount < ntimes; loopCount++) {
		// Read in the data
		pointsRead = fread(diskbuf,sizeof(int16),npoints,fpin);
		if (pointsRead != npoints)
			mexErrMsgTxt("File reading error");
		for (fbi = fourbuf,hbi = diskbuf,j = 0; j < npoints; j++,fbi++,hbi++)
			*fbi = *hbi * scalemult + scaleoff;	// converting to volts
		spctrmD(fourbuf-1,sngD-1,nfreq,navg/2,0);		// PSD with no overlap (for speed)
		if (imHP)
			for (j = 0; j < nfreq; j++,sngD++,logsngD++) {
				*sngD *= nfreq;							// Now the mean across freqs is the mean total power
														// ("sum squared amplitude" normalization in NR parlance)
				*logsngD = log10(*sngD);
			}
		else
			for (j = 0; j < nfreq; j++,sngD++,logsngD++)
				*sngD *= nfreq;
	
		// Graphical output of sonogram
		if (imHP && (loopCount % plotint == 0 || loopCount == ntimes-1)) {
			rhs[0] = imHP;
			rhs[1] = stringStrMxP;
			rhs[2] = logsngP;
			mexCallMATLAB(0,lhs,3,rhs,"set");	// MATLAB command: set(imageh,'CData',logsng)
			mexEvalString("drawnow");
		}
	}
	mxDestroyArray(stringStrMxP);
	mxDestroyArray(logsngP);
	mxFree(diskbuf);
	mxFree(fourbuf);

	return sngP;
}