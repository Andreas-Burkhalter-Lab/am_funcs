function [tcc,indx] = crosscorrspike_unsorted(t1,t2,tmax,nbins)
% CROSSCORRSPIKE_UNSORTED: cross-correlations for unsorted "event" data
%
% The syntax is identicial to crosscorrspike, but it first sorts the data
% and handles the index gymnastics for going back to the original
% representation.
%
% See also: CROSSCORRSPIKE.

% Copyright 2009 by Timothy E. Holy
  
  binning = nargin > 3;

  if ~issorted(t1)
    [t1s,t1so] = sort(t1);
  else
    t1s = t1;
    if ~binning
      t1so = 1:length(t1);
    end
  end
  if ~issorted(t2)
    [t2s,t2so] = sort(t2);
  else
    t2s = t2;
    if ~binning
      t2so = 1:length(t2);
    end
  end
  if ~binning
    [tcc,indx] = crosscorrspike(t1s,t2s,tmax);
    indx = [t1so(indx(1,:)); t2so(indx(2,:))];
  else
    tcc = crosscorrspike(t1s,t2s,tmax,nbins);
  end
end
