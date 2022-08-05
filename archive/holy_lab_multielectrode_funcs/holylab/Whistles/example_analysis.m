% First compute a sparse sonogram: throw out the "hiss," the broadband
% noise that comes from the detector. First we have to figure out how to
% set the threshold: run the following lines for a few seconds, and then
% interrupt with Ctrl-C. From the graph (the power spectrum in a 1ms
% slice), pick a value for the threshold that is just above the baseline
% noise.
p.threshold=0
p.nfreq=256
p.freqrange=[25000 100000]
p.plot=1
sound2sng('file1.bin',p,'file1.sng')
% OK, now we can pick our threshold, so let's save the sparse sonogram
p.threshold=.28
p.plot=0
sound2sng('file1.bin',p,'file1.sng')
% Time to look at it. The window at the bottom is a browser.
spsngplot('file1.sng')
% Let's calculate the beginning and end of each syllable
wt = whistimes('file1.sng',whistimesdefaults)';
save('file1.whistimes','wt','-ascii','-double');
% Compute some parameters of each syllable
p = classifywhis('file1.sng')
% Categorize the syllables by their pitch jumps
[syl,usyl,n] = whistagsyllables(p);
