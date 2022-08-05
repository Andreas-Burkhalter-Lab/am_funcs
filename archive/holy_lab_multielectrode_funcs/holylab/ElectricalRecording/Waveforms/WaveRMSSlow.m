function ds = WaveRMSSlow(filename,decfac,trange)
% WAVERMSSlow: take rms value of waveform in bins
% ds = WaveAbs(filename,decfac,trange)
% where
%        decfac is the decimation factor
%        trange (optional) is the time range to use (in seconds)
% The output (ds) is a numch-by-npts matrix
[fid,message] = fopen(filename,'r');
if (fid < 1)
        error(message);
end
h = ReadAIHeader(fid);
if (nargin < 3)
        trange = [1 h.nscans];
else
        trange = IntersectIntervals(round(trange*h.scanrate),[1 h.nscans]);
end
npts = trange(2)-trange(1)+1;
nptsout = ceil(npts/decfac);
ds = zeros(h.numch,nptsout);
chunksize = decfac*ceil(50000/decfac);
start = trange(1);
startout = 1;
firsttime = 1;
while (start < trange(2))
        currchunk = IntersectIntervals([start,start+chunksize-1],trange);
        currchunkout = IntersectIntervals([startout,startout-1+chunksize/decfac],[1 nptsout]);
        %currchunk
        %currchunkout
        d = ReadBinaryData(fid,h.numch,currchunk-currchunk(1));
        if (firsttime)
                mn = mean(d');
        end
        for i = 1:h.numch
                dtemp = reshape(d(i,:)-mn(i),decfac,size(d,2)/decfac);
                %size(dtemp)
                ds(i,currchunkout(1):currchunkout(2)) = sqrt(mean(dtemp.^2));
        end
        start = start+chunksize;
        startout = startout+chunksize/decfac;
end
fclose(fid);
