function [y,nbrList,S,nbrList0] = meanshift_errorbars(x,xerr,y,options)
% meanshift_errorbars: perform meanshift on a set of points with known uncertainties
%
% In cases where observations come with errorbars, one has an intrinsic
% notion of statisical goodness-of-fit. Consequently, defining a
% neighborhood of points for any given model is straightforward. This
% algorithm flows any candidate model ("probe point") to its local-density
% maximum.
%
% Syntax:
%   [yout,nbrList,S,nbrList0] = meanshift_errorbars(x,xerr,y,options)
% where
%   x is a d-by-N matrix containing a set of observations (i.e., data
%     points, N points in d dimensions)
%   xerr is a d-by-N matrix containing the uncertainty in each coordinate
%     of these observations
%   y is a d-by-1 vector containing a starting "model"
% options is a structure which is either of the form
%     options.mode = 'manual';
%     options.thresh = thresh;  % where threshold is a scalar for chi^2
% or
%     options.mode = 'chi2';    % use chi^2 statistics
%     options.pvalue = pvalue;  % pvalue will be converted into a threshold
% or
%     options.mode = 'T2';   % use Student-T statistics rather than chi^2
%     options.n = n;  % # of trials for each coordiante (can be a scalar,
%       vector, or matrix)
%     options.pvalue = pvalue;
% pvalue is the threshold for goodness of fit (e.g., 0.01)
%
% On exit,
%   yout is the model consistent with the largest number of points (not
%     necessarily the starting model, however)
%   nbrList is the list of neighbors "visible" from yout
%   S is the statistic computed for each point in x with the model yout
%   nbrList0 is the list of neighborhs "visible" from the starting position
%     y
%
% See also; cn_flow_by_distance.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end
  options = default(options,'mode','chi2','pvalue',0.05);
  switch options.mode
    case 'chi2'
      options.thresh = chi2inv(1-options.pvalue,size(x,1));
    case 'T2'
      if any(options.n < 5)
        warning('T2:invalidn','Adjusting n to a minimum of 5, this model cannot handle n<5');
        options.n(options.n < 5) = 5;
      end
      thresh = cn_neighborhoodstatistics(options.n,size(x,1),options.pvalue,'diagonal',0);
      options.thresh = thresh(end);
    case 'manual'
      % Do nothing
    otherwise
      error('mode not recognized')
  end
    
  w = 1./xerr.^2;  % for computing the weighted mean
  
  history = cn_neighborhoodhistory;
  first = true;
  while true
    % Compute the statistic
    dx = bsxfun(@minus,x,y);
    S = sum((dx./xerr).^2,1);
    % Determine the points that fit within the threshold of significance
    % (and their ordering)
    Stmp = S; Stmp(S >= options.thresh) = inf;
    [Stmps,sortOrder] = sort(Stmp);
    nbrList = sortOrder(~isinf(Stmps));
    % Add to history so we can check to see if we are done
    if first
      nbrList0 = nbrList;
      first = false;
    end
    history = history.add(nbrList);
    if history.isAtMax
      break
    end
    % Find the new model
    y = sum(w(:,nbrList).*x(:,nbrList),2)./sum(w(:,nbrList),2);
  end
end
