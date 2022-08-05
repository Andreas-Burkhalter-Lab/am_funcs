function [duration,toff,stimduration] = ephysvalidatetimerange(ephys)
% EPHYSVALIDATETIMERANGE: check for consistent time range across repeats
% [duration,toffset,stimduration] = ephysvalidatetimerange(ephys)
% ephys is a structure array of length nrepeats
% duration is the length of the interval in seconds
% toffset is the time offset (valve opening) in seconds
% stimduration is a vector of stimulus duration, in seconds
%
% See also: EPHYS
scanr = reshape([ephys.scanrange],2,length(ephys));
duration_secs = (diff(scanr)+1)./[ephys.scanrate];
if any(abs(duration_secs - duration_secs(1)) > 0.001)
  error('Not all segments are of the same duration');
end
duration = duration_secs(1);
% Fix case where times are in seconds (this is a hack!)
if all([ephys.scanrate] == 1)
  duration = duration-1;
end

if isfield(ephys,'toffset')
  toff = [ephys.toffset];
  if isempty(toff)
    toff = zeros(1,nrpts);
  end
  if any(toff ~= toff(1))
    error('Do not all have the same toffset');
  end
  toff = toff(1);
else
  toff = 0;
end
if (nargout > 2)
    stimduration = zeros(size(ephys));
    nephys = prod(size(ephys));
    for i = 1:nephys
        stimduration(i) = diff(ephys(i).stimulus(2,[2 3]))/ephys(i).scanrate;
    end
end