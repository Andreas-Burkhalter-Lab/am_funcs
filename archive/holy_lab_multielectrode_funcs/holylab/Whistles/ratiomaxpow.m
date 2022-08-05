function [ratio,t] = ratiomaxpow(filename,trange,band)
% RATIOMAXPOW: calculate power ratio in bands, to detect whistles
% [ratio,t] = ratiomaxpow(filename,trange,band)
% where
%    filename is the name of the .bin file (AI type, raw waveform recording)
%    trange is the range of time to use, in seconds (default: whole file)
%    band is the 4-vector defining the signal and control bands, [sig ctl],
%         in kHz (default: [30 95 20 30])
% and
%    ratio is the ratio of max(sig)/max(ctl)
%    t is the vector of time markers for ratio

% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>

  
if (nargin < 3)
        band = [30 95 20 30];
end
[fid,message] = fopen(filename,'r');
if (fid < 1)
        error(message);
end
h = ReadAIHeader(fid);
if (nargin < 2 | isempty(trange))
        trange = [0 h.nscans/h.scanrate];
end

status = fseek(fid,h.headersize,'bof');        % Move to the beginning of data
if (status)
        error(ferror(fid))
end
status = fseek(fid,round(trange(1)*h.scanrate)*2,'cof'); % Move to beginning of chosen range
if (status)
        error(ferror(fid))
end
start = trange(1);
% Process data in chunks
nfft = 2^8;
npts = round(diff(trange)*h.scanrate);
% Set up output
nt = ceil(npts/nfft);
ratio = zeros(1,nt);
t = linspace(trange(1),trange(2),nt);
f = linspace(0,h.scanrate/2000,nfft/2+1);
fsi = find(f >= band(1) & f < band(2));
fci = find(f >= band(3) & f < band(4));
counter = 1;
start = 0;
% Loop through file
while (start < npts)
        z = ReadBinaryData(fid,1,[0 nfft-1]);
        B = fft(z);
        Bs = max(abs(B(fsi)));
        Bc = max(abs(B(fci)));
        ratio(counter) = Bs/Bc;
        start = start+nfft;
        counter = counter+1;
end
