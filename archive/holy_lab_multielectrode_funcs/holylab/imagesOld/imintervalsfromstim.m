function [indx,identity,vlvsout] = imintervalsfromstim(ip,trange,vlvs)
% [indx,identity,vlvsout] = imintervalsfromstim(imphysin,trange,vlvs)
% where
%   imphysin is a structure array of type IMPHYS
%   trange is in seconds, and gives the [tbegin tend] relative to valve
%     openings---or, if trange is set to the string 'open', just the time
%     range in which the valve is open is used.
%   vlvs (optional) is the list of good valve numbers (if absent or set
%     to 'all', chooses all non-zero values).
% and
%   indx is a cell array of vectors, each giving indices of all elements
%     in imphysin corresponding to a chosen valve transition;
%   identity is a vector of the same size, each holding the
%     valve number associated with the given transition
%   vlvsout is a vector of unique valve numbers that were used
%
% See also: IMPHYS.
  
  if ~isstruct(ip)
    error('imphysin must be a structure (see IMPHYS)');
  end
  if (nargin < 3)
    vlvs = 'all';
  end
  if ischar(vlvs)
    if strcmp(vlvs,'all')
      vlvs = setdiff([ip.stimulus],0);
    else
      error(['Don''t recognize vlvs string ' vlvs]);
    end
  end
  vlvs = intersect(vlvs,[ip.stimulus]);
  vlvsout = vlvs;
  % Parse trange
  justopen = 0;
  if ischar(trange)
    if strcmp(trange,'open')
      justopen = 1;
    else
      error(['Don''t recognize trange string ' trange]);
    end
  end
  % Find all transitions
  newindx = find(diff([ip.stimulus])) + 1;
  % Loop over transitions; if it's of the selected type, then fill in the
  % range
  indx = {};
  identity = [];
  for i = 1:length(newindx)
    cindx = newindx(i);
    if ~isempty(find(vlvs == ip(cindx).stimulus))
      if justopen
        if (i < length(newindx))
          indx{end+1} = newindx(i):newindx(i+1)-1;
        else
          error('Selected valve runs beyond end of file');
        end
      else
        indx{end+1} = find(...
            [ip.stacktime] >= ip(cindx).stacktime+trange(1) & ...
            [ip.stacktime] <= ip(cindx).stacktime+trange(2));
      end
      identity(end+1) = ip(cindx).stimulus;
    end
  end
