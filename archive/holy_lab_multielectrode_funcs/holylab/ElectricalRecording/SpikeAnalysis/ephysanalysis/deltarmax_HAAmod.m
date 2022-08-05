function [deltar,deltarerr,tagsout,rtrial,ntrials] = ...
  deltarmax_HAAmod(ephys,tags,trange,fieldname,options)
%tmin,pre)
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
% HAA -- Added ability to calculate the deltarmax for the time interval
%     prior to stimulus delivery instead of calculating the average firing
%     rate
%   options contains the following fields:
%     tmin: sets the minimum time interval from stim onset to calculate
%      drmax over.  Default = 1
%     pre: is a number that gives the time interval over which to measure 
%      deltarmax prior to stimulus delivery.  If absent, will calculate baseline
%      using average firing rate.  Default = 0
%     dr_over_t: if true, returns ddr (deltaR/time) and ddrerror instead of
%      dr and drerr.  Default = false
%     set_time: if true, calculate deltaR at the same time for each repeate
%      base on what time gives you the max dr.  Default = false
%     dr_slide: if true, calcuate deltaRmax using a sliding time window of
%      tmin second intervals.  Default = false
%
%
%

% See also: DELTARATE, DELTARMAX_CALCULATION, EPHYS, EPHYSSUBRANGE, INTERVALSFROMSTIM.
  
% Copyright 2008 by Timothy E. Holy


if (nargin < 5)
  options = struct;
end
options = default(options,'tmin',1,'pre',0,'dr_over_t',false,'set_time',false,'dr_slide',false); 

% 
% % dtmin default
% if (nargin < 6)
%   pre = 0;
% end
%   
% if (nargin < 5)
%   tmin = 0;
% end

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
ddr = zeros(ntags,nc);
drerr = zeros(ntags,nc);
ddrerr = zeros(ntags,nc);
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
      ttrial_tmp = nan(nrpts,2);  % times of rmax
      for j = 1:nrpts
        cur = sindx(j);
        % Check to make sure that the cell is defined in this time interval
        if checkcellranges && ~isequal(ephys(cur).scanrange,ephys(cur).cellscanrange{cindex})
          continue
        end
        % Convert to secs, shift to correct offset
        zspike = ephys(cur).(fieldname){cindex};
        zspike = (zspike - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
        
        if options.pre ~= 0
          if options.dr_slide == 1
            [rtrial_tmp(j,2),ttrial_tmp(j,2),rtrial_tmp(j,1),ttrial_tmp(j,1)] ...
              = deltarmax_slide(zspike - tsplit(2),tsplit,options);
          else          
          % Compute the rate in the "pre" and "post" period using drmax over 15 sec
          [rtrial_tmp(j,2),ttrial_tmp(j,2),rtrial_tmp(j,1),ttrial_tmp(j,1)]...
            = deltarmax_calculation(zspike-tsplit(2),options);
          end
        else        
          if options.dr_slide == 1
            [rtrial_tmp(j,2),ttrial_tmp(j,2),] = deltarmax_slide(zspike - tsplit(2),tsplit,options);
          else
            % Compute the rate in the "post" period
            [rtrial_tmp(j,2),ttrial_tmp(j,2)] = deltarmax_calculation(zspike - tsplit(2),options);
            %           [rtrial_tmp(j,2),ttrial_tmp(j,2),rtrial_tmp(j,1),ttrial_tmp(j,1)] ...
            %             = deltarmax_calculation(zspike - tsplit(2),tmin);
          end
          % Compute the rate in the "pre" period
          flag = zspike >= tsplit(1) & zspike < tsplit(2);
          rtrial_tmp(j,1) = sum(flag) / dtsplit(1);
        end
      end

      % Now calculate drmax using a set time
      if options.set_time == true
        max_rtrial = max(rtrial_tmp(:,2));
        if max_rtrial == 0 
          rtrial_tmp(j,2) = 0;
        else
          tmax = ttrial_tmp(find(rtrial_tmp(:,2) == max_rtrial),2);
          if length(tmax) > 1
            tmax = tmax(1);
          end
          for j= 1:nrpts
            cur = sindx(j);
            zspike = ephys(cur).(fieldname){cindex};
            zspike = (zspike - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
          
            if isempty(zspike) % There are no spikes
              rtrial_tmp(j,2) = 0;
            else
              spike = zspike(zspike >0);
              spike = spike(spike <= tmax);
              if isempty(spike)
                rtrial_tmp(j,2) = 0;
              else
                n = 1:length(spike);
                rtrial_tmp(j,2) = n(length(n))/tmax;
              end
            end
          end
        end
      end
      
      rtrial{i,cindex} = rtrial_tmp;
      ttrial{i,cindex} = ttrial_tmp;
      % Average across trials
      dr(i,cindex) = nanmean(diff(rtrial_tmp,1,2));
      ddr(i,cindex) = nanmean((diff(rtrial_tmp,1,2)./ttrial_tmp(:,2)));
      ntrials(i,cindex) = sum(~isnan(rtrial_tmp(:,1)));
      drerr(i,cindex) = nanstd(diff(rtrial_tmp,1,2))/sqrt(ntrials(i,cindex));
      ddrerr(i,cindex) = nanstd((diff(rtrial_tmp,1,2)./ttrial_tmp(:,2)))/sqrt(ntrials(i,cindex));
    end
  end
end
if (nargout > 2)
  tagsout = tags;
end
if options.dr_over_t == true
  deltar = ddr;
  deltarerr = ddrerr;
else
  deltar = dr;
  deltarerr = drerr;
end
  
  
