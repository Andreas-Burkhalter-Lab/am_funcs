function tIndex = maxsum_monotonic(R)
% MAXSUM_MONOTONIC: find max of sum of firing rates evaluated at decreasing times
%
% Syntax:
%   tIndex = maxsum_monotonic(R)
% where
%   R is K-by-N, where K is the number of time bins and N is the number
%   of different concentrations; R(i,n) contains the cumulative firing
%   rate from the start of the trial until the ith timeslice, for the nth
%   concentration.
% and
%   tIndex is 1-by-N, containing the time index for each concentration
%   that solves the optimization problem described below.
%
% Given a matrix of firing rates R(timeindex,stimindex), compute a vector
% tIndex(stimindex) such that the sum of R(tIndex(stimindex),stimindex)
% over all stimindex is maximal. The tIndex are constrained to be
% monotonic in stimindex, tIndex(stimindex+1) <= tIndex(stimindex).
% While the algorithm is generic, this is intended to cover the situation
% in which stimulus concentration increases with stimindex, so the peak
% in firing occurs earlier with increasing stimindex.
%
% This can be solved efficiently with a dynamic programming algorithm:
% suppose with n of the concentrations, we know the optimum settings for
% any set of intervals up to a timeindex of T1 (T1 ranges from 1 to
% nbins, the number of timebins). Then by considering all pairs T2 >= T1,
% we can find optimum settings for n+1 of the concentrations as a
% function of T2.  This algorithm is quadratic in nbins, but is linear in
% the number of stimuli.
%
% See also: RMAX_SERIES.
  
% Copyright 2008 by Timothy E. Holy
  
  [nbins,nstim] = size(R);
  % Initialize for n = 1
  tIndex = ones(nbins,nstim);  % we have to maintain T-intermediate results
  rmax = R(1,end);  % Start at the highest concentration
  rsum = nan(nbins,1);
  rsum(1) = rmax;
  for T = 2:nbins
    if (R(T,end) > rmax)
      rmax = R(T,end);
      tIndex(T,end) = T;
    else
      tIndex(T,end) = tIndex(T-1,end);
    end
    rsum(T) = rmax;
  end
  % Add lower concentrations iteratively
  for tagI = nstim-1:-1:1
    rsum_old = rsum;
    tIndex_old = tIndex;
    rmax = rsum(1) + R(1,tagI);
    rsum(1) = rmax;
    for T2 = 2:nbins  % an index on the "new" column
      maxfound = false;
      for T1 = 1:T2   % an index on the "old" columns
        rtmp = rsum_old(T1) + R(T2,tagI);
        if (rtmp > rmax)
          rmax = rtmp;
          tIndex(T2,tagI:end) = [T2 tIndex_old(T1,tagI+1:end)];
          maxfound = true;
        end
      end
      if ~maxfound
        % Copy over the previous setting
        tIndex(T2,tagI:end) = tIndex(T2-1,tagI:end);
      end
      rsum(T2) = rmax;
    end
  end
  % The optimum solution is the final row of tIndex
  tIndex = tIndex(end,:);
end
