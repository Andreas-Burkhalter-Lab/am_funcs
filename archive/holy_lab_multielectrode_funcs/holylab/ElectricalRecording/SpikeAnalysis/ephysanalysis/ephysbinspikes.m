function [npb,tagsout] = ephysbinspikes(ephys,tags,tsplit,fieldname)
% EPHYSBINSPIKES: bin spikes from ephys structure
% [npb,tagsout] = ephysbinspikes(ephys,tags,binbounds,fieldname)
% where
%   ephys is a structure array of type EPHYS (1-by-n)
%   tags is a string or cell array of strings which determines which elements
%     of ephys are studied. If the tags = 'all', then all tags are used.
%   tsplit specifies the time bin boundaries in seconds, relative to
%     valve opening. For example, [-10 0 5] would compute the difference
%     in rates between the periods [-10 0] and [0 5], where 0 is the time
%     of valve opening (determined by ephys.toffset). It can be good to
%     offset these times by a small factor, say 0.1, so that any stimulus
%     artifacts have a better chance of cancelling. For example, [-10 0
%     10]+0.1 would be a good choice.
%   fieldname determines which field (sniptimes or celltimes) is used. If
%     not specified, it will assume celltimes unless celltimes is not present.
%
%   npb is a 1-by-ntags cell array, containing an
%     ncells-by-nrpts-by-ntimebins matrix with the number of spikes in each
%     bin.
%   tagsout is the cell array of tag strings used (useful for tags = 'all')
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
use_celltimes = strcmp(fieldname,'celltimes');
% tags default behavior
if strcmp(tags,'all')
  tags = unique({ephys.tags});
end
if ischar(tags)
  tags = {tags};
end
ntags = length(tags);
% Get number of channels/cells, and check that it's the same
% throughout ephys
nc = length(ephys(1).(fieldname));  % Number of channels/cells
for i = 2:length(ephys)
  if (length(ephys(1).(fieldname)) ~= nc)
    error(['Not all ',fieldname,' are of the same length']);
  end
end
% Check that tsplit is contained in the time intervals of ephys
for i = 1:length(ephys)
  if (tsplit(1) - ephys(i).toffset < 0 | ...
      tsplit(end) - ephys(i).toffset > diff(ephys(i).scanrange)/ ...
               ephys(i).scanrate)
    error('The times in tsplit exceed the bounds of the data');
  end
end
% Get started
npb = cell(1,ntags);
ntimes = length(tsplit)+1;
for i = 1:ntags
  sindx = strmatch(tags{i},{ephys.tag},'exact');
  nrpts = length(sindx);
  if isempty(sindx)
    npb{i} = zeros(nc,nrpts,ntimes);
  else
    [duration,toff] = ephysvalidatetimerange(ephys(sindx));
    for cindex = 1:length(ephys(1).(fieldname))   % channum or cellnum
      % Convert to secs, shift to correct offset
      for j = 1:nrpts
        cur = sindx(j);
	if use_celltimes && ~isequal(ephys(cur).cellscanrange{cindex}, ...
				     ephys(cur).scanrange)
	  % Cells may not be valid over the entire given range, if so
          % enter NaNs for the values in these intervals
	  npb{i}(cindex,j,:) = nan;
	else
        %zspike = getfield(ephys(cur),fieldname,{cindex});
	zspike = ephys(cur).(fieldname){cindex};
	zspike = (zspike{1} - ephys(cur).scanrange(1)) / ephys(cur).scanrate + ephys(cur).toffset;
        npb{i}(cindex,j,:) = HistSplit(zspike,tsplit);
      end
    end
  end
end
if (nargout > 2)
  tagsout = tags;
end
