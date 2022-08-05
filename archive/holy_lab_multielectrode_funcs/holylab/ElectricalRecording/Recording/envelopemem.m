function envmm = envelopemem(v,n)
% envelopemem(v,n)
% Calculates an "envelope,"
% consisting of min/max pairs over blocks of data.
% Very useful for plotting large data sets
% v is n_channels-by-T (T = # of samples), n is the blocksize over which
% min/max pairs are calculated.
%
% See fillmm for plotting

xacc = 0;
nint = floor(n);
nfrac = n - nint;
npts = ceil(size(v,2)/n);
numch = size(v,1);
minv = zeros(numch,npts);
maxv = zeros(numch,npts);
j = 1;
left = 1;
while (j <= npts)
    xacc = xacc + nfrac;
    nload = nint + floor(xacc);
    if (left+nload > size(v,2))
      nload = size(v,2) - left;
    end
    if (nload > 0)
      xacc = xacc - floor(xacc);
      mtemp = v(:,left:left+nload-1);
      envmm.min(:,j) = min(mtemp');
      envmm.max(:,j) = max(mtemp');
    end
    j = j+1;
    left = left+nload;
end