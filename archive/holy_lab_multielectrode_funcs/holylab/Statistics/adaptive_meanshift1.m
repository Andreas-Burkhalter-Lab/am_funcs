function [xc,R2,indexContrib,msd] = adaptive_meanshift1(x,x0,options)
% ADAPTIVE_MEANSHIFT1: shift a single landmark once by adaptive meanshift
% Syntax:
%   xc = adaptive_meanshift1(x,x0)
%   [xc,R2] = adaptive_meanshift1(x,x0)
%   [xc,R2,indexContrib,msd] = ...
%   ... = adaptive_meanshift1(x,x0,options)
% where
%   x is a d-by-N matrix of data points in d-dimensional space;
%   x0 is a column vector of length d, the input position of the
%     landmark;
%   options is a structure with the following fields:
%     factor (default 3): the number of times that the step size needs to
%       exceed the standard error by in order to be considered
%       statistically significant;
%     minN (default 10): the minimum number of neighbors to use when
%       initially considering whether the step is statistically
%       significant. (It's possible for the final number of neighbors to
%       be smaller than this, but only if statistical significance is
%       first triggered beyond this minimum number.)
% and
%   xc is the output position of the landmark;
%   R2 is the squared distance to the farthest point contributing to the
%     new mean position;
%   indexContrib is a list of points that contribute to the new
%     mean position;
%   msd is the vector of mean squared displacements (one for each
%     coordinate) of points that contribute to the new mean position.
%
% Note density can be estimated in terms of n_contrib/(V_d R^d), where d
% is the dimensionality and V_d is the volume coefficient for a ball in d
% dimensions.
%
% References:
%   See the paper by T. E. Holy, ???.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'minN')
    options.minN = 10;   % Don't check the first points initially (noisy)
  end
  if ~isfield(options,'factor')
    options.factor = 3;  % The number of standard errors needed to be sure
  end
  if ~isfield(options,'backtrack')
    options.backtrack = true;
  end
  if ~isfield(options,'protect_n')
    options.protect_n = false;
  end
  if ~isfield(options,'covariance')
    options.covariance = false;
  end

  [d,N] = size(x);
  dx = x - repmat(x0,1,N);
  sd = sum(dx.^2,1);  % Square distance between x and x0
  [ssd,sort_order] = sort(sd);   % Sort by increasing distance
  dxs = dx(:,sort_order);
  % Calculate both the step and its uncertainty. For the following we'd
  % have to divide by n to get the true step, but there's no need until
  % we've picked the correct R.
  dxsum = cumsum(dxs,2);   % Sum over points; except for /n, this is the step
  if ~options.covariance
    stepsqr = sum(dxsum.^2,1);  % Sum over dimensions
    %mse_component = cumsum(ssd);
    cummsd = cumsum(dxs.^2,2);
    mse_component = sum(cummsd,1);
  else
    cov_elements = repmat(reshape(dxs,[d 1 N]),[1 d 1]) .* ...
      repmat(reshape(dxs,[1 d N]),[d 1 1]);
    cumcov = cumsum(cov_elements,3);
    isdone = false;
    stepsqr = [];
    while ~isdone
      curindex = length(stepsqr)+1;
      Cinv = ams1inv(cumcov(:,:,curindex));
      stepsqr(curindex) = dxsum(:,curindex)' * Cinv * dxsum(:,curindex);
      mse_component(curindex) = 1;
      isdone = stepsqr(curindex) > options.factor^2;
    end
  end
  % Now compare the square step against the mean square error in the
  % step. Find an R where it's clear that the step first exceeds the
  % expected error by a substantial margin.
  chosen_sure = find(stepsqr(options.minN:end) > ...
    options.factor^2*mse_component(options.minN:end),...
    1,'first');
  if isempty(chosen_sure)
    chosen = length(sd);
  else
    chosen_sure = chosen_sure + options.minN;  % correct for excluding first pts
    if options.backtrack
      % Now that we know it's well in excess of our threshold, back up a
      % bit to the point at which it just crosses the standard error.  This
      % way we pick up the point at which the mean first starts to come out
      % of the noise
      chosen = find(stepsqr(1:chosen_sure-1) <= ...
        mse_component(1:chosen_sure-1),...
        1,'last');
      if isempty(chosen)
        %chosen = options.minN;
        chosen = 2;
      else
        chosen = chosen+1;
      end
      %chosen = max(chosen,options.minN);
    else
      chosen = chosen_sure;
    end
  end
  if options.protect_n
    if isnumeric(options.protect_n)
      chosen_min = round(options.nOld*options.protect_n);
    elseif ischar(options.protect_n)
      if strcmp(lower(options.protect_n),'sqrt')
        chosen_min = options.nOld - round(sqrt(options.nOld));
      else
        error('options.protect_n not recognized');
      end
    end
    chosen = max(chosen,chosen_min);
  end
  
  xc = x0+dxsum(:,chosen)/chosen;  % Calculate the new center
  R2 = ssd(chosen);
  msd = cummsd(:,chosen)/(chosen-1);
  indexContrib = sort_order(1:chosen);
  
  
  
function Cinv = ams1inv(C)
  [U,S,V] = svd(C);
  sdiag = diag(S);
  sdiag(sdiag == 0) = Inf;
  Cinv = V*diag(1./sdiag)*U';
  