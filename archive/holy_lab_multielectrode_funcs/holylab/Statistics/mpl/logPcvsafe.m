function [lPcv,isOK] = logPcvsafe(Pcv,X,outlier_frac,is_sorted)
% LOGPCVSAFE: safely compute log(Pcv), even when Pcv == 0
% The idea here is that when Pcv == 0, we can approximate it as exp(-|dX|)
% in one dimension, or exp(-|dX|^2/2) in multidimensions, where |dX| is the
% distance to the closest point.  Thus, if we also have the distances on
% hand, then we can estimate log(Pcv) even when Pcv is zero.
%
% This function also allows you to "throw out" outliers by setting the
% "thresh" parameter described below.
%
% Syntax:
%   lPcv = logPcvsafe(Pcv,X)   (one dimension)
%   lPcv = logPcvsafe(Pcv,sd)  (multidimensions)
%   lPcv = logPcvsafe(...,outlier_frac)
%   lPcv = logPcvsafe(...,outlier_frac,is_sorted)   (one dimension)
% where
%   Pcv is the cross-validated density (a vector)
%   X (one dimensions) is the position of each point, or its projection
%     (1-by-M vector)
%   sd (multidimensions) is the square distance between each pair of
%     points (a M-by-M matrix)
%   outlier_frac: to make the calculation robust, you can specify that a
%     certain fraction of the points are to be treated as outliers, and
%     discarded in terms of their contribution to the log-likelihood. For
%     example, if outlier_frac is 0.05, then the 5% of points with the
%     lowest logPcv will be discarded.
%   is_sorted: in one dimension, you can save time by setting this to true
%     if your data are already sorted.
% and
%   lPcv is the (robust) estimate of log(Pcv)
%   isOK is 1 for the points that are retained, and 0 for the outliers.

% Copyright 2006 by Timothy E. Holy

  is_multidim = ~isvector(X);
  M = size(X,2);
  is_robust = nargin > 2;
  isz = (Pcv == 0);
  lPcv(~isz) = log(Pcv(~isz));
  z_index = find(isz);
  if ~isempty(z_index)
    % First check to see whether we're discarding outliers, and if so,
    % whether the number of points with Pcv=0 is less than or equal to
    % the number of points to be discarded.  When that's true, we can
    % save time by not worrying about these points
    if (is_robust && length(z_index) <= outlier_frac*M)
      lPcv(z_index) = -Inf;
    else
      % OK, we have to do the estimation.
      % For the ones where Pcv == 0, estimate Pcv = exp(-dX/2), where
      % dX is the distance to the closest point (get rid of the /2
      % in one dimension)
      if is_multidim
        dX = zeros(1,length(z_index));
        for i = 1:length(z_index)
          dX(i) = min(X(z_index(i),:)); % X is matrix of square distances
        end
        lPcv(z_index) = -dX/2;
      else
        if (nargin > 3 && is_sorted)
          Xs = [-Inf X Inf];
          z_index_s = z_index+1;
        else
          [Xs,sort_order] = sort(X);
          [sso,inv_sort_order] = sort(sort_order);
          Xs = [-Inf Xs Inf];  % protect the edges
          z_index_s = inv_sort_order(z_index)+1;
        end
        dX = min([Xs(z_index_s+1) - Xs(z_index_s); ...
                  Xs(z_index_s) - Xs(z_index_s-1)]);
        lPcv(z_index) = -dX;
      end
    end
  end
  if is_robust
    % OK, throw out a number of the most-outlying points
    slPcv = sort(lPcv);
    threshIndex = round(outlier_frac*M);
    if (threshIndex > 0)
      thresh = slPcv(threshIndex);
      isOK = lPcv > thresh;
      lPcv(~isOK) = 0;
    else
      % outlier_frac was too small to throw out even a single point
      isOK = logical(ones(size(Pcv)));
    end
  end
  lPcv = reshape(lPcv,size(Pcv));
  if is_robust
    isOK = reshape(isOK,size(Pcv));
  end