function t = sortheader_absolute_sniptime(sh,index)
% SORTHEADER_ABSOLUTE_SNIPTIME: Get the absolute spike time in seconds
% Syntax:
%   t = sortheader_absolute_sniptime(shc,index)
% where
%   shc is a structure array of channel-specific sortheaders (see
%     SORTHEADER_IMPORTCHAN);
%   index is a cell array of snippet indices, one for each entry in shc
%     (can instead be a vector if only one snippet file is used) or the
%     string 'all' if you want them all;
% and
%   t is a cell array of vectors of spike times, in absolute seconds
%     (i.e., not relative to the start of the file, but relative to some
%     fixed time, e.g. Jan 1 in 0A.D.).  t will be a vector, rather than a
%     cell array of vectors, if index was supplied as a vector.
%
% See also: SORTHEADER_ABSOLUTE_STARTTIME, SORTHEADER_IMPORTCHAN.
  
% Copyright 2005 by Timothy E. Holy
  
  nfiles = length(sh);
  outcell = 1;
  useall = 0;
  if ischar(index)
    useall = 1;
  elseif (nfiles == 1 && ~iscell(index))
    index = {index};
    outcell = 0;
  end
  tstart = sortheader_absolute_starttime(sh);
  for i = 1:nfiles
    if useall
      t{i} = double(sh(i).sniptimes')/sh(i).scanrate + tstart(i);
    else
      t{i} = double(sh(i).sniptimes(index{i})')/sh(i).scanrate + ...
             tstart(i);
    end
  end
  if (~outcell & length(t) == 1)
    t = t{1};
  end