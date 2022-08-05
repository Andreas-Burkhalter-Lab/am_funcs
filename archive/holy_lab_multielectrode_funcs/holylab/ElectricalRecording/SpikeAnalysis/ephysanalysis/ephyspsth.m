function [rout,tout,rerrout,npb] = ephyspsth(ephys,binwidth,index,fieldname)
% EPHYSPSTH: compute peri-stimulus time histogram
% Syntax:
%   [rout,tout,rerrout,npb] = ephyspsth(ephys,binwidth,index,fieldname)
% where
%   ephys is the input ephys structure (should only contain elements
%     describing epochs you want included in the PSTH; i.e., all elements
%     should have identical stimulus content & duration);
%   binwidth is the width of the bins, in seconds;
%   index is a scalar giving the index of the channel or cell to be
%     analyzed;
%   fieldname (optional) is a string containing 'sniptimes' or
%     'celltimes', depending on whether you want multi-unit or
%     single-unit data, respectively;
% and
%   rout is a vector containing the mean firing rate as a function of
%     time, in Hz;
%   tout is the time of each bin center, in seconds;
%   rerrout is the standard error in the mean firing rate, in Hz;
%   npb is the number of spikes per bin on each trial.
%
% Note this will divide the interval into evenly-spaced bins, so input
% binwidth is only approximate.

% Copyright 2001 by Timothy E. Holy

% 2007-06-22 (RCH): modified so will expunge any intervals which have "NaN"
%        instead of any celltimes at all (designed to accomodate cells
%        whose sorting is only defined over some of the intervals held in a
%        large structure)
% 2006-03-26 (RCH): altered removal of spike times from bins so that 
%        if there is only one set of spiketimes and it is not within an
%        extra cell structure it will still be read out (useful if eg
%        unsorted data from VNO neurons is being treated like celltimes)
  
if (nargin < 4)
  if isfield(ephys,'celltimes')
    fieldname = 'celltimes';
  else
    fieldname = 'sniptimes';
  end
end
check_celltimes = strcmp(fieldname,'celltimes') && isfield(ephys,'cellscanrange');

% if (ischar(index) & strcmp(index,'all'))
%   if strcmp(fieldname,'sniptimes')
%     index = 1:length(ephys(1).channels);
%   else
%     index = 1:length(ephys(1).cellnumbers);
%   end
% end
nrpts = length(ephys);
[duration,toff] = ephysvalidatetimerange(ephys);
zspike = cell(1,nrpts);
% Convert to secs, shift to correct offset, and identify NaN intervals
isnanint = zeros(1,nrpts);
for i = 1:nrpts
  %zspike{i} = getfield(ephys(i),fieldname,{index});
  if index==1 & ~iscell(ephys(i).(fieldname))
	  zspike{i} = ephys(i).(fieldname);
  else
	  zspike{i} = ephys(i).(fieldname){index};
  end
  zspike{i} = (zspike{i} - ephys(i).scanrange(1)) / ephys(i).scanrate;
  if isnan(zspike{i})
      isnanint(i) = 1;
  end
  if check_celltimes && ~isequal(ephys(i).scanrange, ...
      ephys(i).cellscanrange{index})
    isnanint(i) = 1;
  end
end
% expunge NaN intervals...
zspike(find(isnanint)) = [];
nrpts = length(zspike);
% get back to business!
nbins = round(duration/binwidth);
binwidth = duration/nbins;
tsplit = (1:nbins-1)*binwidth;
npb = zeros(nrpts,nbins);
for i = 1:nrpts
  npb(i,:) = HistSplit(zspike{i},tsplit);
end
if (nrpts > 1)
  rout = mean(npb)/binwidth;
  rerrout = std(npb)/(binwidth*sqrt(nrpts-1));
else
  rout = npb/binwidth;
  rerrout = NaN*npb;
end
tout = [binwidth/2,tsplit + binwidth/2] + toff;
if nrpts == 0 % means we expunged all of them!
    rout = NaNs(1,nbins);
    rerrout = NaNs(1,nbins);
    npb = NaNs(1,nbins);
end