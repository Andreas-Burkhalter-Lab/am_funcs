function index = timesininterval(times,intervals)
% TIMESININTERVAL: select those indices within particular intervals
% index = timesininterval(times,intervals)
% where
%   times is a vector of times. MUST be sorted on input.
%   intervals is an n-by-2 matrix of time ranges, [tstart tend].
%     No particular order required.
% and
%   index is a 1-by-n cell array of indices, times(index{i}) are the
%     times that fall in the range intervals(i,:).
%
% The times are specified in half-open fashion, i.e.
% times falling within the range [tstart tend) are retained.
  [nintervals,sizeint] = size(intervals);
  if (sizeint ~= 2)
    error('intervals must be a n-by-2 matrix');
  end
  if (nintervals == 0)
    index = {zeros(0,1)};
    return
  end
  index = cell(1,nintervals);
  times = times(:);
  % Define the included region to be half-open, [tstart tend)
  % This is the reason for swapping orders in the 2 sorts
  [tbeg,tbegi] = sort([intervals(:,1);times]);
  [tend,tendi] = sort([times;intervals(:,2)]);
  nt = length(times);
  % Find the marks to the interval boundaries
  ibegi = find(tbegi <= nintervals);
  iendi = find(tendi > nt);
  % Map these marks back into the original time vector
  % (don't be confused by inclusion of interval boundaries)
  ibeg = ibegi - (1:nintervals)';
  iend = iendi - (1:nintervals)';
  % Permute these back to the original input order
  [dummy,pbeg] = sort(tbegi(ibegi));
  [dummy,pend] = sort(tendi(iendi));
  ibeg = ibeg(pbeg);
  iend = iend(pend);
  % Select out the indices
  for i = 1:nintervals
    index{i} = ibeg(i)+1:iend(i);
  end
