function [dch,sdsng] = specdiscont(sng)
% SPECDISCONT: estimate the spectral discontinuity
% Syntax:
%  d = specdiscont(sng)
% where
%   sng is a sonogram of size nfreq-by-ntimes (may be sparse), 
% and
%   d is the spectral discontinuity at each time slice, of size
%     1-by-(ntimes-1).
% The spectral discontinuity is defined as the sum of the abs of the
% difference between spectra in two adjacent time slices, with two
% modifications:
%    1. The total power is normalized in each time bin
%    2. Each spectrum is allowed to shift by a small integral multiple of
%       the frequency bin size.
% The value of d in a time bin ranges from 0 to 2

% Copyright Timothy E. Holy 2004-08-11

offset = -3:3; % The number of bins in either direction
moff = max(abs(offset));
frange = 1+moff:size(sng,1)-moff;
pow = sum(sng)';
[i,j,s] = find(sng);
% The "normalized" sonogram
% (sum of each column is 0 or 1)
nsng = sparse(i,j,s./pow(j),size(sng,1),size(sng,2));
sdsng = zeros(length(offset),size(sng,2)-1);
for i = 1:length(offset)
  dsng = nsng(frange,2:end) - nsng(frange+offset(i),1:end-1);
  sdsng(i,:) = sum(abs(dsng));
end
% Now find the discontinuity for the best shift at each time slice
dch = [min(sdsng),2];
zindx = find(pow == 0);
dch(zindx) = 2;   % For bins with no data, set to highest discontinuity
