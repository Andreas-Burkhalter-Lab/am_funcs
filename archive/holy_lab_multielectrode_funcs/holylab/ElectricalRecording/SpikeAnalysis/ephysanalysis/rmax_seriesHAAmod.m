function [rtrial,tpeak] = rmax_series(ephys,tags,trange,fieldname,nbins,tmin,offset)
% MAX_SERIES: compute neural response to concentration series
% Syntax:
%   [rtrial,tpeak] = rmax_series(ephys,tags,trange,fieldname,nbins,tmin)
% where
%   ephys is a structure array of type EPHYS (1-by-n)
%   tags is a string or cell array of strings which determines which elements
%     of ephys are studied. THESE MUST BE IN ORDER OF INCREASING CONCENTRATION.
%   trange specifies the time bin boundaries in seconds, relative to
%     valve opening. For example, [-10 0 5] would compute the difference
%     in rates between the periods [-10 0] and [0 5], where 0 is the time
%     of valve opening (determined by ephys.toffset). 
%   fieldname determines which field (sniptimes or celltimes) is used. If
%     not specified, it will assume celltimes unless celltimes is not
%     present.
%   nbins (default 20): the number of time bins used to separate the
%     response period [tmin trange(3)]. This affects computation time and
%     how finely-grained you measure your firing rates.
%   tmin (default 0): specifies the minimum interval needed to measure rates.
%
%   rtrial is a cell array of size ntags-by-ncells (or
%     ntags-by-nchannels), with each element containing a
%     matrix of size nrepeats-by-2 giving the response size on each
%     trial; the first column is the baseline (pre-stimulus) and the
%     second column is the response (peri-stimulus).
%   tpeak is an ntags-by-ncells matrix, for each cell/channel the set
%     of times at which the firing rate is evaluated for each stimulus.
%
% The response size is measured in the following way: all spikes that
% occur within the selected time interval are used to calculate a
% _cumulative_ firing rate from the beginning of the interval.  For each
% concentration, the time at which the maximum firing rate (averaged over
% all trials) occurs is noted, and denoted as tpeak. rtrial is then the
% firing rate on each trial, evaluated from the beginning of the interval
% to the tpeak for that concentration.
%
% The one subtlety is the following: tpeak is constrained to be
% non-increasing with concentration.
%
% The baseline is measured as the rmax with a tmin that is at least as
% large as max(tpeak).  For this reason, you want the baseline period to
% be at least the same length as the response period.
%
% See also: DELTARATE, DELTARMAX, MAXSUM_MONOTONIC, DELTARMAX_CALCULATION, EPHYS, EPHYSSUBRANGE, INTERVALSFROMSTIM.
  
% Copyright 2008 by Timothy E. Holy

  % tmin default
  if (nargin < 6)
    tmin = 0;
  end

  % nbins default
  if (nargin < 5)
    nbins = 20;
  end
  if (nbins < 2)
    nbins = 2;
  end

  % fieldname default behavior
  if (nargin < 4)
    if isfield(ephys,'celltimes')
      fieldname = 'celltimes';
    else
      fieldname = 'sniptimes';
    end
  end

  % parse the tags
  if ischar(tags)
    tags = {tags};
  end
  ntags = length(tags);

  % Check the trange
  trange = trange(1:3);
  dtrange = diff(trange);
  if (diff(dtrange) > 0)
    error('Baseline period must be at least as large as the response period');
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

  % Find the ephys indices that correspond to each tag, and validate the
  % time ranges
  tagIndex = cell(1,ntags);
  for tagI = 1:ntags
    sindx = strmatch(tags{tagI},{ephys.tag},'exact'); 
    tagIndex{tagI} = sindx;
    [duration,toff] = ephysvalidatetimerange(ephys(sindx));
    if (trange(1) < toff)
      error('The left boundary of trange was not included in the ephys trange');
    end
    if (trange(end) > duration+toff)
      error('The right boundary of trange was not included in the ephys trange');
    end
  end

  % Create the bin boundaries for the response period
%   tbin = [0 linspace(tmin,diff(trange([2 end])),nbins+tmin)]';  
%   nbins = length(tbin);
  rtrial = cell(ntags,nc);
  tpeak = nan(ntags,nc);
  for cellIndex = 1:nc
    % For the current cell/channel, collect the spike times on all trials,
    % referencing them to the beginning of the interval.  Then compute the
    % cumulative firing rate
    
    %MOVED FROM ABOVE HAA
    tbin = [0 linspace(tmin,diff([offset(cellIndex) trange(end)]),nbins+tmin)]';  % MOD HAA
    nbins = length(tbin);
    %%%% End edits   
    
    rtrial_binned = cell(1,ntags);
    rbinned = nan(nbins-1,ntags);
    tspike_baseline = cell(1,ntags);
    for tagI = 1:ntags
      % Extract the spikes
      nrpts = size(tagIndex{tagI},1);
      nspikes = nan(nbins,nrpts);
      keepFlag = false(1,nrpts);
      tspike_baseline{tagI} = cell(1,nrpts);
      for rptI = 1:nrpts
        % Check to make sure that the cell is defined in this time interval
        if checkcellranges && ~isequal(ephys(rptI).scanrange,ephys(rptI).cellscanrange{cellIndex})
          continue
        end
        keepFlag(rptI) = true;
        % Convert to secs, shift to correct offset relative to the ephys interval
        tspike = ephys(tagIndex{tagI}(rptI)).(fieldname){cellIndex};
        if ~isempty(tspike)
            tspike = (tspike(:) - ephys(tagIndex{tagI}(rptI)).scanrange(1)) / ...
                ephys(tagIndex{tagI}(rptI)).scanrate + ephys(tagIndex{tagI}(rptI)).toffset;
            % Extract spikes in the "baseline" interval
            flag = tspike >= trange(1) & tspike < offset(cellIndex);  %trange(2); HAA CHANGED
            tspike_baseline{tagI}{rptI} = tspike(flag) - trange(1);
            % Extract spikes in the "response" interval
            flag = tspike >= offset(cellIndex) & tspike < trange(3); %trange(2) & tspike < trange(3); HAA CHANGED
            tspike = tspike(flag) - offset(cellIndex); % trange(2); HAA CHANGED  
            % Count the spikes in the response interval in bins
            if ~isempty(tspike)
                nspikes(:,rptI) = histc(tspike,tbin);
            else
                nspikes(:,rptI) = 0;
            end
        else
            tspike_baseline{tagI}{rptI} = [];
            nspikes(:,rptI) = 0;
        end
      end
      nspikes = nspikes(:,keepFlag);
      tspike_baseline{tagI} = tspike_baseline{tagI}(keepFlag);

      nrpts = size(nspikes,2);
      if (nrpts > 0)
        % Calculate the firing rate (measured from the beginning of the
        % interval) on each trial
        ntmp = cumsum(nspikes,1);
        rtrial_binned{tagI} = ntmp(1:end-1,:) ./ ...
          repmat(tbin(2:end),1,nrpts);
        % Compute the average firing rate across trials
        rbinned(:,tagI) = mean(rtrial_binned{tagI},2);
      end
    end

    % Find the optimum index
    tIndex = maxsum_monotonic(rbinned);
    % Evaluate the peak time, and calculate the firing rate on each trial
    tpeak(:,cellIndex) = tbin(tIndex+1)';
    for tagI = 1:ntags
     % indx = sub2ind(size(rtrial_binned{tagI}),tIndex,1:ntags); 
      rtrial{tagI,cellIndex} = rtrial_binned{tagI}(tIndex(tagI),:)';
    end

    % Evaluate the baseline rate using the longest time
    for tagI = 1:ntags
      nrpts = length(tspike_baseline{tagI});
      if (nrpts > 0)
        rbase = zeros(nrpts,1);
        tbase = sort(cat(1,tspike_baseline{tagI}{:}));
        [rmax,tmax] = deltarmax_calculation(tbase',struct('tmin',max(tpeak(:,cellIndex))));
        for rptI = 1:nrpts
          rbase(rptI) = sum(tspike_baseline{tagI}{rptI} <= tmax)/tmax;
        end
        rtrial{tagI,cellIndex} = [rbase rtrial{tagI,cellIndex}];
      end
    end
  end
end
