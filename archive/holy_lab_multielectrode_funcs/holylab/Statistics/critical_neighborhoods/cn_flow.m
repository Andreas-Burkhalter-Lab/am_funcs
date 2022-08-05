function [nbrList,triggerSource,x0,C,n] = cn_flow(x,w,params,T2thresh,x0,C)
% cn_flow: move a point to its peak by critical neighborhood mean shift
%
% Syntax:
%   [nbrList,triggerSource,x0,C,n] = cn_flow(x,w,params,T2thresh,x0,C0)
% where
%   x is the d-by-N matrix of data points
%   w (optional) is a 1-by-N vector of weights (supply [] if all 1)
%   params is a structure with the following fields:
%     .nMin (optional): the minimum number of points per neighborhood
%       (default 1)
%     .covarianceModel (optional): 'isotropic', 'diagonal', or 'full'
%       (default 'isotropic', or inferred from C0). See dist_mahalanobis.
%     .updateCovariance (optional): if true, the distance to points is
%       computed using the most recent neighborhood sample covariance
%       matrix (default true)
%   T2thresh is 1-by-N vector containing the threshold for significance for
%     the T2 statistic for neighborhoods of any number of points
%   x0 is the starting location (a d-by-1 vector)
%   C0 (optional) is the initial covariance matrix to use for calculating
%     distances. Default is the identity. See dist_mahalanobis for the format.
%
% See also: dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy
  
% Note: this version does not yet use pointwise statistics
  
  %% Parse arguments
  [d,N] = size(x);
  useWeight = ~isempty(w);
  params = default(params,'nMin',1,'updateCovariance',true);
  covarianceModelIndex = strmatch(params.covarianceModel,{'isotropic','diagonal','full'});
  if isempty(covarianceModelIndex)
    error('covarianceModel not recognized')
  else
    covarianceModelIndex = covarianceModelIndex-1;
  end  
  % Initialize covariance
  if ~isfield(params,'covarianceModel') || isempty(params.covarianceModel)
    if isempty(C)
      params.covarianceModel = 'isotropic';
    else
      if (numel(C) == 1)
        params.covarianceModel = 'isotropic';
      elseif (numel(C) == d)
        params.covarianceModel = 'diagonal';
      elseif (numel(C) == d*d)
        params.covarianceModel = 'full';
      else
        error('Can''t determine which covariance model is being used');
      end
    end
  end
  if (nargin < 6 || isempty(C))
    switch params.covarianceModel
      case 'isotropic'
        C = 1;
      case 'diagonal'
        C = ones(d,1);
      case 'full'
        C = eye(d,d);
      otherwise
        error('Covariance model not recognized');
    end
  end
        
  %% Do the flow
  history = NeighborhoodHistory;
  nOld = 0;
  while ~history.isAtMax
    % Order the points by increasing distance
    dx = bsxfun(@minus,x,x0);
    z2 = dist_mahalanobis(dx,C,params.covarianceModel);
    [~,sortOrder] = sort(z2);
    % Compute the statistic
    if useWeight
      T2N = chisqNP(dx,w,sortOrder,covarianceModelIndex);
    else
      T2N = chisqNP(dx,sortOrder,covarianceModelIndex);
    end
    % Find the first point that exceeds significance
    triggerSource = 1;
    T2Check = T2N; T2Check(1:params.nMin-1) = 0;
    firstIndex = find(T2Check > T2thresh,1,'first');
    if isempty(firstIndex)
      firstIndex = N+1;
      triggerSource = 0;
    end
    n = firstIndex-1;
    % Update the neighborhood parameters and covariance used in the metric
    nbrList = sortOrder(1:n);
    xsel = x(:,nbrList);
    if useWeight
      wsel = w(:,nbrList(1:n));
      [n,x0,Cnew] = moments_from_points(xsel,params.covarianceModel,wsel);
    else
      [n,x0,Cnew] = moments_from_points(xsel,params.covarianceModel);
    end
    if triggerSource == 0
      C = Cnew;
      return
    end
    if params.updateCovariance
      % Do it so that radical decreases in number do not unnecesarily
      % change C
      if (n > nOld/2)
        C = Cnew;
      end
%       C = (nOld*C + n*Cnew)/(nOld + n);
      nOld = n;
    end
    % Add the neighborhood to history, to test whether we are done
    history = history.add(nbrList);
  end
  n = history.n;
  