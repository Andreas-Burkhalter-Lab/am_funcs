function [nbrList,triggerSource,x0,C,n] = cn_flow_gmm(x,psum,params,T2thresh,x0,C)
% cn_flow_gmm: move a point to its peak by critical neighborhood mean shift
% using an algorithm appropriate for Gaussian mixture models
%
% Syntax:
%   [nbrList,triggerSource,x0,C,n] = cn_flow_gmm(x,psum,params,T2thresh,x0,C0)
% where
%   x is the d-by-N matrix of data points
%   psum is a 1-by-N vector of the sum-of-likelihoods of the _other_ gaussian
%     components (supply all zeros if this is the first component)
%   params is a structure with the following fields:
%     .nMin (optional): the minimum number of points per neighborhood
%        (default 1)
%     .covarianceModel (optional): 'isotropic', 'diagonal', or 'full'
%        (default 'isotropic', or inferred from C0). See dist_mahalanobis.
%     .updateCovariance (optional): if true, the distance to points is
%        computed using the most recent neighborhood sample covariance
%        matrix (default true)
%     .initialW1 (default true): if true, on the first iteration the GMM
%        weights are set to 1; after that, the covariance is used to help
%        set the GMM weights. If false, the supplied covariance C0 (see
%        below) is used on the first iteration (set to false only if the
%        initial covariance is meaningful).
%   T2thresh is 1-by-N vector containing the threshold for significance for
%     the T2 statistic for neighborhoods of any number of points
%   x0 is the starting location (a d-by-1 vector)
%   C0 (optional) is the initial covariance matrix to use for calculating
%     distances. Default is the identity. See dist_mahalanobis for the format.
%
% See also: cn_gmm, dist_mahalanobis.

% Copyright 2011 by Timothy E. Holy
  
% Note: this version does not yet use pointwise statistics
  
  %% Parse arguments
  [d,N] = size(x);
  params = default(params,'nMin',1,'updateCovariance',true,'initialW1',true);
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
  useC = ~params.initialW1;
        
  %% Do the flow
  history = NeighborhoodHistory;
  nOld = 0;  % used as a sanity check in updating the covariance for the distance computation
  triggerSourceOld = -1;
  while ~history.isAtMax
    % Order the points by increasing distance
    dx = bsxfun(@minus,x,x0);
    z2 = dist_mahalanobis(dx,C,params.covarianceModel);
    [~,sortOrder] = sort(z2);
    % Compute the weight
    if useC
      p = exp(-z2/2)/sqrt((2*pi)^d * det_cov(C,d,params.covarianceModel));
      w = p./(p+psum);
      w(p == 0) = 0;  % handle NaNs
    else
      w = ones(1,N);
    end
    % Compute the statistic
    T2N = chisqNP(dx,w,sortOrder,covarianceModelIndex);
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
    wsel = w(:,nbrList);
    [n,x0,Cnew] = moments_from_points(xsel,params.covarianceModel,wsel);
    % If we've taken the whole data set twice in a row, quit
    if triggerSource == 0 && triggerSourceOld == 0
      C = Cnew;
      return
    end
    triggerSourceOld = triggerSource;
    useC = true;
    if params.updateCovariance
      % Do it so that radical decreases in number do not unnecesarily
      % change C
      if (n > nOld/2)
        C = Cnew;
      end
%       C = (nOld*C + n*Cnew)/(nOld + n);
      nOld = n;
    else
      % Just update the magnitude of the covariance---because we need this
      % for setting w---but don't change any other aspects
      C = C * (det_cov(Cnew,d,params.covarianceModel)/det_cov(C,d,params.covarianceModel))^(1/d);
    end
    % Add the neighborhood to history, to test whether we are done
    history = history.add(nbrList,round(n));
  end
  n = history.n;
  