function [dr,drerr,tagsout,rtrial,ntrials] = deltarmax(ephys,tags,trange,fieldname,tmin)
% DELTARMAX: compute rate difference upon stimulation
% [dr,drerr,tagsout,rtrial,ntrials] = deltarate(ephys,tags,trange,fieldname,tmin)
% where
%   ephys is a structure array of type EPHYS (1-by-n)
%   tags is a string or cell array of strings which determines which elements
%     of ephys are studied. If the tags = 'all', then all tags are used.
%   trange specifies the time bin boundaries in seconds, relative to
%     valve opening. For example, [-10 0 5] would compute the difference
%     in rates between the periods [-10 0] and [0 5], where 0 is the time
%     of valve opening (determined by ephys.toffset). It can be good to
%     offset these times by a small factor, say 0.1, so that any stimulus
%     artifacts have a better chance of cancelling.  For example, [-10 0
%     10]+0.1 would be a good choice.  trange may be supplied as a matrix
%     of size (ntags-by-3), allowing different timeranges to be used for
%     each tag.
%   fieldname determines which field (sniptimes or celltimes) is used. If
%     not specified, it will assume celltimes unless celltimes is not
%     present.
%   tmin (default 0): specifies the minimum interval needed to measure rates.
%
%   dr is a ntags-by-ncells (or ntags-by-nchannels) matrix giving the
%     rate changes, where rate is measured as the mean inverse interspike
%     interval;
%   drerr is the sem of dr across trials
%   tagsout is the cell array of tag strings used (useful for tags =
%     'all')
%   rtrial is a cell array of length ntags-by-ncells (or
%     ntags-by-nchannels), with each element containing a
%     nrepeats-by-2 matrix giving the average inverse interspike interval
%     for each trial in the two periods.
%   ntrials is a ntags-by-ncells (or ntags-by-nchannels) matrix giving
%     the number of valid trials for each stimulus/cell combination.
%     This may not be the same for all cells, since each cell can have
%     its own start/stop times.
%
% See also: DELTARATE, DELTARMAX_CALCULATION, EPHYS, EPHYSSUBRANGE, INTERVALSFROMSTIM.
  
% Copyright 2008 by Timothy E. Holy

% dtmin default
if (nargin < 5)
  tmin = 0;
end

% fieldname default behavior
if (nargin < 4)
  if isfield(ephys,'celltimes')
    fieldname = 'celltimes';
  else
    fieldname = 'sniptimes';
  end
end

% tags default behavior
if strcmp(tags,'all')
  tags = unique({ephys.tag});
end
if ischar(tags)
  tags = {tags};
end
ntags = length(tags);
if (ntags > 1 && size(trange,1) == 1)
  trange = repmat(trange,ntags,1);
end
if (size(trange,1) ~= ntags)
  error('Number of tags and number of time ranges do not match');
end
% Get number of channels/cells, and check that it's the same
% throughout ephys
nc = length(ephys(1).(fieldname) );  % Number of channels/cells
for i = 2:length(ephys)
  if (length(ephys(i).(fieldname) ) ~= nc)
    error(['Not all ',fieldname,' are of the same length']);
  end
end
% Decide if we need to check cell ranges
checkcellranges = false;
if isfield(ephys,'cellscanrange')
  checkcellranges = true;
end
% Get started
dr = zeros(ntags,nc);
drerr = zeros(ntags,nc);
ntrials = zeros(ntags,nc);
for i = 1:ntags
  sindx = strmatch(tags{i},{ephys.tag},'exact');
  if isempty(sindx)
    dr(i,:) = NaN;
    drerr(i,:) = NaN;
  else
    nrpts = length(sindx);
    [duration,toff] = ephysvalidatetimerange(ephys(sindx));
    if (trange(1) < toff)
        error('The left boundary of trange was not included in the ephys trange');
    end
    if (trange(end) > duration+toff)
      error('The right boundary of trange was not included in the ephys trange');
    end
    tsplit = trange(i,:);
    dtsplit = diff(tsplit);
    for cindex = 1:length(ephys(1).(fieldname) )   % channum or cellnum
      rtrial_tmp = nan(nrpts,2);
      for j = 1:nrpts
        cur = sindx(j);
        % Check to make sure that the cell is defined in this time interval
        if checkcellranges && ~isequal(ephys(cur).scanrange,ephys(cur).cellscanrange{cindex})
          continue
        end
        % Convert to secs, shift to correct offset
        zspike = ephys(cur).(fieldname){cindex};
        zspike = (zspike - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
        % Compute the rate in the "pre" period
        flag = zspike >= tsplit(1) & zspike < tsplit(2);
        rtrial_tmp(j,1) = sum(flag) / dtsplit(1);
	% Compute the rate in the "post" period
    options.tmin = tmin;
	rtrial_tmp(j,2) = deltarmax_calculation(zspike - tsplit(2),options);
      end
      rtrial{i,cindex} = rtrial_tmp;
      % Average across trials
      dr(i,cindex) = nanmean(diff(rtrial_tmp,1,2));
      ntrials(i,cindex) = sum(~isnan(rtrial_tmp(:,1)));
      drerr(i,cindex) = nanstd(diff(rtrial_tmp,1,2))/sqrt(ntrials(i,cindex));
    end
  end
end
if (nargout > 2)
  tagsout = tags;
end
