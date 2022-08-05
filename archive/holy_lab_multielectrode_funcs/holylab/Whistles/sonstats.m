function sts = sonstats(filenames,tsplit,frange)
% SONSTATS: compute summary statistics for sonograms
% Computes the total power in a given frequency band
% over the time, using only those bins in which the
% power high/power low ratio is greater than 1.
%
% sts = sonstats(filenames,tsplit,frange)
% where
%   filenames is a cell array of *sng.mat files
%   tsplit (optional or empty) contains the bin boundaries
%     in time, in seconds (default: one bin containing entire file)
%   frange (optional, though must supply tsplit if supplying
%      frange) is a 4-vector of [bandlow bandhigh], in kHz.
%      The ratio of mean powers over bandhigh/bandlow is
%      computed, and only those timebins with ratio>1 are
%      included in the total power. Only the power in bandhigh
%      contributes.
%      (default: [25 35 40 95]). 
  
% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>

if (nargin < 3)
  frange = [25 35 40 95];
end
if (nargin < 2)
  tsplit = [];
end
if ~iscell(filenames)
  filenames = {filenames};
end
bandsig = frange(3:4);
bandctl = frange(1:2);
nfiles = length(filenames);
sts = zeros(nfiles,length(tsplit)+1);
for i = 1:nfiles
  load(filenames{i});
  [ratsng,t] = SngRatio(sng,p,bandsig,bandctl);
  deltat = diff(t(1:2));
  % Compute the total power over regions where the ratio > 1
  ff = linspace(0,p.scanrate/2000,p.nfreq);        % Converts to kHz
  fsigi = find(ff >= bandsig(1) & ff <= bandsig(2));
  ti = find(ratsng > 1);
  % Now split these regions into time bins
  %    First convert times into indices
  [ts,indx] = sort([t tsplit]);
  isplit = find(indx > length(t));
  if ~isempty(isplit)
    isplit = isplit - (0:length(isplit)-1);
  end
  %    Second, parcel out the power
  if (isempty(isplit))
    sts(i,1) = sum(sum(sng(fsigi,ti)));
  else
    ilast = 0;
    for j = 1:length(isplit)
      ij = ti(find(ti >= ilast & ti < isplit(j)));
      sts(i,j) = sum(sum(sng(fsigi,ij)));
      ilast = isplit(j);
    end
    ij = ti(find(ti >= ilast));
    sts(i,end) = sum(sum(sng(fsigi,ij)));
  end
  % Make normalization equal to mean square amplitude
  sts(i,:) = sts(i,:)/p.nfreq;        % (un-does a dumb thing in C source code)
  binwidth = diff([0 tsplit p.tacq]);
  sts(i,:) = deltat*sts(i,:)./binwidth;
end
