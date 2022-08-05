// Take the rms of waveform over bins (thus decimating waveform)
// Compute the mean value using only the first ~5s of data
// (this should suffice to eliminate zero-offset on amplifier
//  & A/D card)
// Note: The end point is not handled correctly yet when the last
// packet is of substandard size, so I just throw that out now.
//#define CHECKBOUNDS 1
#include <vector>
#include "iotypes.h"
#include "numertypes.h"
#include "FileData.h"
#include "mex.h"

using namespace std;

const long default_nwind = 5*10000; // Number of scans to process at one time

extern void WaveRMSDec(FilePtr &fpin,double *dout,const int32 numCh,const int32 nscans,const vector<int> &chIndex,int decimate);

void WaveRMSDec(FilePtr &fpin,double *dout,const int32 numCh,const int32 nscans,const vector<int> &chIndex,int decimate)
{
  //  Make sure that nwind is a multiple of decimate
  long nwind = (default_nwind/decimate)*decimate;     

  ScanWindow<int16> sw(nscans,numCh,nwind);
  vector<SkipPtr<int16> > channelP(chIndex.size());
  vector<double> means(chIndex.size(),0);
  long i,size;
  double *douttemp;
  int firsttime;

  firsttime = 1;
  douttemp = dout;
  // Loop over disk reads
  do {
    // cout << '.' << flush;
    fpin >> sw;
    // Set up pointers to each channel's data
    for (i = 0; i < chIndex.size(); i++)
      channelP[i] = sw.window(chIndex[i]);
    for (i = 0; i < channelP.size(); i++)
      channelP[i] = channelP[i].begin();

    SkipPtr<int16> &first = channelP[0];
    SkipPtr<int16> tmpp,endp;
    float tmpf;

    // First time through, compute the mean value for each channel
    if (firsttime) {
      for (i = 0; i < channelP.size(); i++) {
		for (; channelP[i] < channelP[i].end(); channelP[i]++)
	  		means[i] += *channelP[i];
		means[i] /= channelP[i].size();
		//mexPrintf("Mean %g:\n",means[i]);
		channelP[i] = channelP[i].begin();         // reset pointer
      }
      firsttime = 0;
    }

    // Loop over packets of size decimate
    long counter = 0;
    while (first + decimate <= first.end()) {
      // Compute size of current packet
      // (it will be decimate, unless it's the very last packet)
      size = min(decimate,int(channelP[0].end() - channelP[0])/int(numCh));
      //if (counter < 20 | size < decimate)
	  //    mexPrintf("size %d, decimate %d, bigsize %d\n",size,decimate,int(channelP[0].end() - channelP[0])/int(numCh));
      // Loop over channels
      for (i = 0; i < channelP.size(); i++) {
		// Compute  variance(this region of data))
		*dout = 0;
		endp = channelP[i]+size;
		for (tmpp = channelP[i]; tmpp < endp; tmpp++) {
	  		tmpf = float(*tmpp) - means[i];
	  		*dout += tmpf*tmpf;
			//if (counter < 20) {
				//mexPrintf("%d %g ",*tmpp,tmpf);
			//	mexPrintf("%d ",*tmpp);
			//	counter++;
			//}
		}
		//if (counter < 20)
		//	mexPrintf("\nsumsq %g\n",*dout);
		*dout = sqrt(*dout/size);
		//if (counter < 20)
		//	mexPrintf("\noutput %g\n",*dout);
		dout++;
		channelP[i] += decimate;
      }
    }
  } while (!sw.atEnd());
  //mexPrintf("Ptr diff: %d\n",dout-douttemp);
}

