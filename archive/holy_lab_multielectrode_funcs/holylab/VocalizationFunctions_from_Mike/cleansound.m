function yc = cleansound(y,blocksize,thresh,cfac)
% CLEANSOUND: removes background noise by threshholding in fourier domain
% yc = cleansound(y,blocksize,thresh,cfac)
% cfac = 0 can also be done with a combination of SPARSESNG and SNG2SOUND
if (nargin < 4)
        cfac = 0;
end
N = length(y);
nblocks = floor(N/blocksize);
yc = zeros(1,nblocks*blocksize);
window = triang(blocksize)';
start = 1;
while (start + blocksize -1 <= N)
        cindx = start:start+blocksize-1;
        ycc = cleanblock(y(cindx),thresh,cfac);
        yc(cindx) = yc(cindx) + window.*ycc;
        start = start+blocksize/2;
end

function yc = cleanblock(y,thresh,cfac)
N = length(y);
yt = fft(y);
%yta = abs(yt(2:N/2));
%indx = find(yta > thresh)+1;
ytc = yt;
%ikeep = [1,indx,N/2+1,N+2-indx];
yta = abs(yt);
ired = find(yta < thresh);
ytc(ired) = ytc(ired)*cfac;
yc = real(ifft(ytc));
