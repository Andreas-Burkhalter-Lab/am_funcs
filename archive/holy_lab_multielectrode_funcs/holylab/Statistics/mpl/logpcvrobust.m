function [lPcv,isOK] = logpcvrobust(Pcv,outlier_frac)
% LOGPCVROBUST: robustly compute log(Pcv), given outliers
%
% This function also allows you to "throw out" outliers by setting the
% "outlier_frac" parameter described below.
%
% Syntax:
%   lPcv = logpcvrobust(Pcv)
%   lPcv = logpcvrobust(Pcv,outlier_frac)
% where
%   Pcv is the cross-validated density (a vector)
%   outlier_frac: to make the calculation robust, you can specify that a
%     certain fraction of the points are to be treated as outliers, and
%     discarded in terms of their contribution to the log-likelihood. For
%     example, if outlier_frac is 0.05, then the 5% of points with the
%     lowest logPcv will be discarded.
% and
%   lPcv is the estimate of log(Pcv)
%   isOK is 1 for the points that are retained, and 0 for the outliers.

% Copyright 2006 by Timothy E. Holy

  M = length(Pcv);
  lPcv = -inf(size(Pcv));
  is_robust = nargin > 1;
  isz = (Pcv == 0);
  lPcv(~isz) = log(Pcv(~isz));
  z_index = find(isz);
  if is_robust
    % Throw out a number of the most-outlying points
    [slPcv,sort_order] = sort(lPcv);
    threshIndex = round(outlier_frac*M);
    slPcv(1:threshIndex) = 0;
    isOK = logical([zeros(1,threshIndex) ones(1,M-threshIndex)]);
    [sso,sort_inv] = sort(sort_order);
    lPcv = slPcv(sort_inv);
    isOK = isOK(sort_inv);
  end
  if is_robust
    isOK = reshape(isOK,size(Pcv));
  end
  