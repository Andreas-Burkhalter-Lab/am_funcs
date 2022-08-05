function [dr,drerr,tagsout,rtrial,ntrials,tmax_all] = drmax_acrosstrials(ephys,tags,trange,fieldname,options)
% DRMAX_ACROSSTRIALS: compute rate difference upon stimulation. Finds max dr across all trials
% [dr,drerr,tagsout,rtrial,ntrials,tmax] = drmax_acrosstrials(ephys,tags,trange,fieldname,options);
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
%   options: a structure which may have the following fields:
%     tmin (default 1): specifies the minimum interval needed to measure rates.
%
%   dr is a ntags-by-ncells (or ntags-by-nchannels) matrix giving the
%     rate changes upon stimulation
%   drerr is the sem of dr across trials
%   tagsout is the cell array of tag strings used (useful for tags =
%     'all')
%   rtrial is a cell array of length ntags-by-ncells (or
%     ntags-by-nchannels), with each element containing a
%     nrepeats-by-2 matrix giving the average firing rate
%     for each trial in the two periods.
%   ntrials is a ntags-by-ncells (or ntags-by-nchannels) matrix giving
%     the number of valid trials for each stimulus/cell combination.
%     This may not be the same for all cells, since each cell can have
%     its own start/stop times.
%   tmax is a ntags-by-ncells array giving the time of the maximum
%     cumulative firing rate in the post-stimulus period
%
% See also: DELTARATE, DELTARMAX_CALCULATION, EPHYS, EPHYSSUBRANGE, INTERVALSFROMSTIM.
%
% Oct 2010 HAA


% dtmin default
if (nargin < 5)
    options.tmin = 1;
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
tmax_all = zeros(ntags,nc);
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
        % first get all the spikes for each trial
        for j = 1:nrpts
            cur = sindx(j);
            % Check to make sure that the cell is defined in this time interval
            if checkcellranges && ~isequal(ephys(cur).scanrange,ephys(cur).cellscanrange{cindex})
                continue
            end
            zspike = ephys(cur).(fieldname){cindex};
            % get all the dimensions in the right order
            if isempty(zspike)
                tmpsp{j} = NaN;
            else
                zspike = reshape(zspike,length(zspike),1);
                tmpsp{j} = (zspike - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
            end
        end
        a = cat(1,tmpsp{:});
        sp_trial = sort(a(~isnan(a)),'ascend'); % all spike times across all trials
        % Compute the rate in the "post" period
        ops = options;
        ops.pre = tsplit(1)-tsplit(2);
        [rmax,tmax] = deltarmax_calculation(sp_trial-tsplit(2),ops); % subtract the offset
        tmax = tmax+tsplit(2); % compensate for the offset
        rmax = rmax/nrpts; % the average firing rate across trials at that time point
        tmax_all(i,cindex) = tmax;

        % now use that time information and calculate the firing rates for each trial
        rtrial_tmp = nan(nrpts,2);
        for j = 1:nrpts
            cur = sindx(j);
            zspike = ephys(cur).(fieldname){cindex};
            zspike = (zspike - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
            % Compute the rate in the "pre" period
            flag = zspike >= tsplit(1) & zspike < tsplit(2);
            rtrial_tmp(j,1) = sum(flag) / dtsplit(1);
            % now compute the rate in the "post" period until the max time
            flag = zspike >= tsplit(2) & zspike <= tmax;
            rtrial_tmp(j,2) = sum(flag) / diff([tsplit(2) tmax]);
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
            
           