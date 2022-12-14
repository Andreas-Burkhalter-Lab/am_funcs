function sng = RecordingWhis(channel,scanrate,nfreq,navg,nperblock,nblocks,gain,npix,axH,textH,usrhdr,filename)
% RecordingWork(channels,scanrate,nscans,buffersize,trigger,gain,npix,axH,textH,usrhdr,filename)
% The real version is a MEX file
% channels = vector of channel #s (0-63), in the scanning order
% scanrate = number of scans/second (1 scan = 1 D/A conversion for each channel)
% nscans = total number of scans to acquire
% buffersize = number of scans to process in a chunk
% trigger: 0 = trigger is put out on PFI0, 1 = acquisition waits for trigger on trigger1/PFI0
% gain: valid choices are -1 (for 0.5), 1, 2, 5, 10, 20, 50, 100 (a scalar)
% npix = number of "bins" to decimate over (for graphical output)
% axH = vector of axis handles for the respective channels
% textH = handle for text messages during plotting
% usrhdr = text string (optional), to save with file
% filename = text string (optional), use only if you want to save the data to a file
