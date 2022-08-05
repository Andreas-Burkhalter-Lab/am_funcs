function [Y,centroid] = plsdata(X,nu,options)
% PLSDATA: power-law scale data
% 
% Given a set of data vectors x_i and a scaling parameter nu, return the
% rescaled vectors
%     y_i = (x_i - mu)/|x_i - mu|^(1-nu),
% where mu (the centroid) is set so that the y_i have zero mean.
%
% Syntax:
%   Y = plsdata(X,nu)
%   Y = plsdata(X,nu,options)
%   [Y,centroid] = plsdata(...)
% where
%   X is a nX-by-d data matrix (each row is an observation in d
%     dimensions);
%   nu is the scaling power, see above (default nu = 1, i.e., no scaling);
%   options is a structure array with the following possible fields:
%     nocenter: if present & true, skips the step of subtracting the
%       centroid (mu). 
% and
%   Y contains the shifted & scaled data;
%   centroid contains the mean or other centroid of the data set, or NaN
%     if nocenter was true.
%
% See also: CENTROID_PLS, PCA.

% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 3)
    options = struct;
  end
  if (nargin < 2)
    nu = 1;
  end
  if ~isfield(options,'nocenter')
    options.nocenter = 0;
  end
  [nX,d] = size(X);
  % Compute centroid
  if options.nocenter
    centroid = NaN;
    Y = X;
  else
    centroid = centroid_pls(X,nu);
    Y = X - repmat(centroid,nX,1);
  end
  if (nu ~= 1)
    % Scale by a power of the magnitude
    Ydenom = sum(Y.^2,2).^((1-nu)/2);
    Y = Y ./ repmat(Ydenom,1,d);
  end
