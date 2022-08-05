function [xf,n,n_iter] = cn_filter(x,pvalue,isacausal,covarianceModel)
% cn_filter: filter a signal in time using critical neighborhoods
% Syntax:
%   xf = cn_filter(x,pvalue)
% The simplest syntax, performs causal (looking into the past only)
% filtering on the signal using the specified pvalue (a scalar). x is a
% 1-by-N vector of values, or (if you want to filter a multi-valued signal)
% a d-by-N matrix.  On output, xf has the same size as x.  For multi-valued
% signals, the default for the accumulator is to use an isotropic
% covarianceModel (see below).
%
%   xf = cn_filter(x,pvalue,isacausal)
% Set isacausal true if you want to use neighborhoods that are symmetric with
% respect to past and future, false if you just want to do causal
% filtering.
%
%   xf = cn_filter(x,pvalue,isacausal,covarianceModel)
% For multi-valued signals, this lets you specify the covarianceModel
% explicitly (see dist_mahalanobis).
%
%   [xf,n,n_iter] = cn_filter(x,...)
% This returns the number of points in the neighborhood (n), and the number
% of iterations required to get to the peak (n_iter) for each time point.
%
% See also: dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy

% Implemented as a MEX file
