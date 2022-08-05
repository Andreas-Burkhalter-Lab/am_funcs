% sngparms_default
function x = sngparms_default(threshold,freqrange,nfft)
x = struct('threshold',threshold,'plot',0,'progressbar',0,'freqrange',freqrange,'nfreq',nfft/2);
end