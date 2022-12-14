function [minv,maxv] = envelopens(fid,n,numch,range)
% envelope(fid,n,range)
% Loads binary data and calculates an "envelope,"
% consisting of min/max pairs over blocks of data.
% It is assumed that there are numch channels in the data
% Very useful for plotting large data sets
% See fillmm for plotting
xacc = 0;
nint = floor(n);
nfrac = n - nint;
npts = ceil((range(2) - range(1)+1)/n);
minv = zeros(numch,npts);
maxv = zeros(numch,npts);
status = fseek(fid,numch*range(1)*2,'cof'); % *2 because int16s
if status
        error(ferror(fid))
end
j = 1;
while (j <= npts & ~feof(fid))
        xacc = xacc + nfrac;
        nload = nint + floor(xacc);
        xacc = xacc - floor(xacc);
        mtemp = fread(fid,[numch,nload],'int16');
        minv(:,j) = min(mtemp');
        maxv(:,j) = max(mtemp');
        j = j+1;
end
