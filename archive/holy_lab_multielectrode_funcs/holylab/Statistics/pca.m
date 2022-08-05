function [projectDirections,sqrt_lambda,centroid] = pca(X,options)
% PCA: do Principal Component Analysis and return all components & eigenvalues
% 
% This version of PCA makes no attempt to determine the significance of
% the individual components.
%
% Syntax:
%   [projectDirections,sqrt_lambda,centroid] = pca(X,options)
% where
%   X is a nX-by-d data matrix (each row is an observation in d dimensions);
%   options is a structure array with the following possible fields:
%     maxpoints (default Inf): if nX > maxpoints, a subset of points will
%       be used for pca.
%     nu (default 1): if nu<1, does power-law scaled PCA (PLS-PCA), a
%       variant of PCA more robust to outliers.  This is somewhat slower
%       because the computation of the centroid is more involved than
%       with the mean.
%     nocenter: if present & true, skips the step of subtracting the mean
%       (or other centroid, see above);
% and
%   projectDirections is a d-by-ndims matrix, where each column is a
%     principal component;
%   sqrt_lambda is a vector containing the square root of the
%     eigenvalues of the covariance matrix (note if you're doing PLS-PCA,
%     and you want the eigenvalues in the units of the unscaled data, you
%     may want to use eigenvalues = sqrt_lambda^(2/nu));
%   centroid contains the mean or other centroid of the data set, or NaN
%     if nocenter was true.
%
% See also: CENTROID_PLS, PCA_LANDMARK, PCA_PARALLEL_ANALYSIS.

% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  if ~isfield(options,'maxpoints')
    options.maxpoints = Inf;
  end
  if ~isfield(options,'nu')
    options.nu = 1;
  end
  if ~isfield(options,'nocenter')
    options.nocenter = 0;
  end
  [nX,d] = size(X);
  if (nX > options.maxpoints)
    skip = ceil(nX / options.maxpoints);
    X = X(1:skip:end,:);
    nX = size(X,1);
  end
  % Shift & scale the data
  [dX,centroid] = plsdata(X,options.nu,options);
  % Now do the PCA
  [U,S,projectDirections] = svd(dX,0);
  ndims = min(nX-1,d);
  projectDirections = projectDirections(:,1:ndims);
  sqrt_lambda = diag(S)/sqrt(nX-1);
