function [t,Mc] = resample_evenly_preserve_integral(tc,N)
% RESAMPLE_EVENLY_PRESERVE_INTEGRAL: evenly-resample a set of 1d signals, invariantly
%
% Given a set of signals sampled at different non-uniform time intervals,
% put all of the signals on a common temporal axis by finding a set of
% common times and projection matrices that preserve to total integral.
%
% Syntax:
%   [t,Mc] = resample_evenly_preserve_integral(tc,N)
% where
%   tc is a cell array, each element corresponding to a single "signal",
%     where the element is a list of monotonic-increasing times at which a
%     sample was taken (conceptually, each value is the right edge of the
%     sampling interval)
%   N is the number of intervals over which you'd like the function
%     resampled
% and
%   t is the vector of sampling times, evenly spaced in time
%     (like tc, marking the right edge of an interval)
%   Mc is a cell array of sparse matrices; when a column vector sampled
%     times tc{i} is pre-multiplied by Mc{i}, one obtains the resampled
%     values at regularly-spaced times t
%
% See also: resample_preserve_integral.
  
% Copyright 2009 by Timothy E. Holy
  
  if ~iscell(tc)
    tc = {tc};
  end
  n_samples = length(tc);

  tmin = cellfun(@(p) p(1),tc);
  tmax = cellfun(@(p) p(end),tc);
  dt = cellfun(@(p) median(diff(p)),tc);
  trange = [min(tmin-dt),max(tmax)];
  trange(2) = trange(2)+10*eps(trange(2));
  deltat = diff(trange)/N;
  t = (1:N)*deltat + trange(1);

  Mc = cell(1,n_samples);
  for sampIndex = 1:n_samples
    n = length(tc{sampIndex});
    T = (tc{sampIndex} - trange(1))/deltat;
    Ti = floor(T);
    indexSplit = find(diff(Ti))+1;
    frac = (Ti(indexSplit) - T(indexSplit-1)) ./ (T(indexSplit) - T(indexSplit-1));
    v = ones(size(Ti));
    v(indexSplit) = frac;
    Tis = Ti;
    Tis(indexSplit) = Ti(indexSplit)-1;
    Mc{sampIndex} = sparse([Tis Ti(indexSplit)]+1,[1:n indexSplit],[v 1-frac]);
  end
