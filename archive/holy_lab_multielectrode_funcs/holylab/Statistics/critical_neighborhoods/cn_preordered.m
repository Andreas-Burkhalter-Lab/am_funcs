function T2 = cn_preordered(dx,varargin)
% cn_preordered: neighborhood statistics for a pre-defined inflation sequence
% 
% This function computes the T^2 statistic for collections of points that
% grow by "inflation" (accretion).  One value is computed at each step of
% inflation.
%
% Syntax: T2 = cn_preordered(dx,sortOrder,covarianceModel)
% Syntax: T2 = cn_preordered(dx,w,sortOrder,covarianceModel)
%   where
%     dx is a d-by-N matrix of data points, each data point on a
%       column. These measure the displacement of the data point from a
%       basepoint.
%     w (optional) is a 1-by-N vector of weights, one per data point
%       (default is w=1)
%     sortOrder is a 1-by-N vector, containing the sequence in which the
%       points should be added to the neighborhood (unit-offset as is
%       typical for Matlab)
%     covarianceModel (optional) is 'isotropic', 'diagonal', or 'full'.
%        (Only the first character is used.) Default is 'isotropic'. See
%        dist_mahalanobis.
%   and
%     T2 is a vector of the length of sortOrder.  T2(i) contains the T^2
%       statistic for the displacement of the mean for a neighborhood that
%       includes the points sortOrder(1:i).  Note that if you are using
%       weights, then the significance of T2(i) should be tested using
%       n = round(cw(i)), where cw = cumsum(w), rather than n = i.
%
% See also: cn_neighborhoodstatistics, dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy

% Implemented as a MEX file
