function [dr,drerr,tagsout,npbo,ntrials] = deltarate(ephys,tags,trange,fieldname)
% DELTARATE: compute rate difference upon stimulation
% [dr,drerr,tagsout,npb,ntrials] = deltarate(ephys,tags,trange,fieldname)
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
%     not specified, it will assume celltimes unless celltimes is not present.
%
%   dr is a ntags-by-ncells (or ntags-by-nchannels) matrix giving the
%     rate changes;
%   drerr is the sem of dr
%   tagsout is the cell array of tag strings used (useful for tags =
%     'all')
%   npb is a cell array of length ntags-by-ncells (or
%     ntags-by-nchannels), with each element containing a
%     nrepeats-by-2 matrix giving the number of spikes in each interval
%     for each trial.
%   ntrials is a ntags-by-ncells (or ntags-by-nchannels) matrix giving
%     the number of valid trials for each stimulus/cell combination.
%     This may not be the same for all cells, since each cell can have
%     its own start/stop times.
%
% See also: EPHYS, EPHYSSUBRANGE, INTERVALSFROMSTIM.

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
    for cindex = 1:length(ephys(1).(fieldname) )   % channum or cellnum
      zspike = cell(1,nrpts);
      % Convert to secs, shift to correct offset
      for j = 1:nrpts
        cur = sindx(j);
        %zspike{j} = getfield(ephys(cur),fieldname,{cindex});
        zspike{j} = ephys(cur).(fieldname){cindex};
        zspike{j} = (zspike{j} - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
      end
      tsplit = trange(i,:);
      npb = zeros(nrpts,4);
      for j = 1:nrpts
        tmp = HistSplit(zspike{j},tsplit);
        npb(j,1:length(tmp)) = tmp;
      end
      % Here we check to see if the given cell was defined during each
      % interval---if not, put NaNs in the npb for those trials
      if strcmp(fieldname,'celltimes') && isfield(ephys,'cellscanrange')
        for j = 1:nrpts
          cur = sindx(j);
          if ~isequal(ephys(cur).scanrange,ephys(cur).cellscanrange{cindex})
            npb(j,:) = nan;
          end
        end
      end
      if (nargout > 3)
        npbo{i,cindex} = npb(:,[2 3]);
      end
      % convert the number per bin into a rate
      drtmp = diff(npb(:,2:3)./repmat(diff(trange(i,:)),nrpts,1),1,2);
      nrpts_tmp = sum(~isnan(drtmp));
      dr(i,cindex) = nanmean(drtmp);
      if (nrpts_tmp > 1)
        drerr(i,cindex) = nanstd(drtmp)/sqrt(nrpts_tmp);
      else
        drerr(i,cindex) = NaN;
      end
      ntrials(i,cindex) = nrpts_tmp;
    end
  end
end
if (nargout > 2)
  tagsout = tags;
end
