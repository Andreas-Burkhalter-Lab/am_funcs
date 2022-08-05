function [thresh,fc] = thresh_discrim(x,y)
% THRESH_DISCRIM: compute discriminability as a function of threshold
% Given two vectors of data, x & y, determine the percent-correct
% classification as a function of a threshold.
%
% Syntax:
%   [thresh,fc] = thresh_discrim(x,y)
% where
%   x and y are vectors of observations (not necessarily of the same
%     length);
% and
%   thresh is a vector of thresholds
%   fc is a vector of fraction correct, fc(i) is the fraction correct when
%     choosing a threshold thresh(i).  This fraction correct is the average
%     of the fraction correct for each group, so it is not influenced by
%     differences in the sizes of the two sets of observations.

% Copyright 2006 by Timothy E. Holy
  
  sx = sort(x);
  sx(end+1) = sx(end)+eps;
  sy = sort(y);
  sy(end+1) = sy(end)+eps;
  [ux,xr] = unique(sx);
  [uy,yr] = unique(sy);
  % Compute the fraction of observations that are less than each unique value
  xr = (xr-1)/length(x);
  yr = (yr-1)/length(y);
  thresh = unique([ux(:);uy(:)]');
  xindx = 1;
  yindx = 1;
  fc = zeros(size(thresh));
  for i = 1:length(thresh)
    % Keep xindx & yindx pointing to the corresponding x & y ranks
    while (xindx < length(ux) && ux(xindx) < thresh(i))
      xindx = xindx+1;
    end
    while (yindx < length(uy) && uy(yindx) < thresh(i))
      yindx = yindx+1;
    end
    fc(i) = (xr(xindx) + 1 - yr(yindx))/2;
    fc(i) = max(fc(i),1-fc(i));
  end
