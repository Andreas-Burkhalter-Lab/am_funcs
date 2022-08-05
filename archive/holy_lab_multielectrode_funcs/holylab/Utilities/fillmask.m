function xf = fillmask(x,isgood)
% FILLMASK: replace masked values by interpolated ones
% Syntax:
%   xf = fillmask(x,isgood)
% where
%   x is the signal (one-dimensional) to be "filtered"
%   isgood is a logical vector, 1 where the corresponding value x is to be
%     preserved, and 0 where the value of x is to be interpolated using
%     the nearest good points on either side
% and
%   xf is the output, of the same size as x.
%
% Note: "bad" values (isgood = 0) that are on either edge of the vector
% cannot be interpolated, and will simply be copies of the original values.
% A warning will be issued.
  
% Copyright 2006 by Timothy E. Holy
  
  forwardindex = nan(1,numel(x));
  backwardindex = forwardindex;
  % Filter forward, keeping the most recent good value
  lastgoodindex = nan;
  for i = 1:length(x)
    if isgood(i)
      lastgoodindex = i;
    end
    forwardindex(i) = lastgoodindex;
  end
  % Filter backward
  lastgoodindex = NaN;
  for i = length(x):-1:1
    if isgood(i)
      lastgoodindex = i;
    end
    backwardindex(i) = lastgoodindex;
  end
  % Check to see if we have problems on the edges
  if any(isnan(isgood([1 end])))
    warning('fillmask: Can''t interpolate on edges');
  end
  % Now combine them, interpolating if necessary
  xf = x;
  interpindex = find(forwardindex ~= backwardindex & ...
     ~isnan(forwardindex) & ~isnan(backwardindex));
  alpha = (backwardindex(interpindex) - interpindex)./...
      (backwardindex(interpindex) - forwardindex(interpindex));
  xf(interpindex) = alpha .* x(forwardindex(interpindex)) + ...
      (1-alpha) .* x(backwardindex(interpindex));
