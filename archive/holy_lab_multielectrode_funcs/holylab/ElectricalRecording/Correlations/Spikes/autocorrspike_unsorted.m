function [tac,indx] = autocorrspike_unsorted(t,tmax,nbins)
% AUTOCORRSPIKE_UNSORTED: auto-correlations for unsorted "event" data
%
% The syntax is identicial to autocorrspike, but it first sorts the data
% and handles the index gymnastics for going back to the original
% representation.
%
% See also: AUTOCORRSPIKE.

% Copyright 2009 by Timothy E. Holy

  binning = nargin > 2;
  
  if ~issorted(t)
    [ts,tso] = sort(t);
  else
    ts = t;
    if ~binning
      tso = 1:length(t);
    end
  end
  if ~binning
    [tac,indx] = autocorrspike(ts,tmax);
    indx = tso(indx);
  else
    tac = autocorrspike(ts,tmax,nbins);
  end
end