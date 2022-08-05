function [info,params] = cn_flow_by_distance(x,params,varargin)
%   [info,params] = cn_flow_by_distance(x,params,b)
%   [info,params] = cn_flow_by_distance(x,params,b,C0)
%   [info,params] = cn_flow_by_distance(x,params,info_in)
  
  [d,N] = size(x);
  params = default(params,'updateCovariance',true,'alphaeps',sqrt(eps(class(x))));
  [params,b,C0] = cn_startposition_defaults(params,varargin{:});
  % Initialize T2thresh
  if ~isfield(params,'T2thresh')
    try
      [params.T2thresh,params.T2lsthresh] = cn_T2thresh_frommc(d,N,params);
    catch ME
      params.T2thresh = cn_neighborhoodstatistics(N,d,params.pvalue,params.covarianceModel,params.nMin);
    end
  end
  if ~isfield(params,'T2lsthresh')
    warning('cn:linesearch','The parameters for the line search are not set, expect poor performance!');
    params.T2lsthresh = params.T2thresh;
  end
  haveFunc = isfield(params,'func');

  info = struct('nbrList',[],'T2',[],'xpre',b,'xpost',[],'Cpre',C0,'interiorIndex',[],'n_iter',0);

  alpha = inf;
  nOld = 0;
  history = cn_neighborhoodhistory;
  while (abs(alpha) > params.alphaeps)
    %% Check, and update, the history
    if (info.n_iter > 1)
      if history.isAtMax
        break
      end
    end
    if (info.n_iter > 0)
      history = history.add(info.nbrList);
    end
    
    info.n_iter = info.n_iter+1;
    info.Cpre = C0;

    %% Shift the base point
    [info.nbrList,alpha,b,Cnew,alphamax] = cn_baseshift_by_distance(x,params,b,C0);
    info.xpost = b;
    [~,info.interiorIndex] = min(alphamax);
    
    %% User interaction
    if haveFunc
      params = params.func(b,info.nbrList,params);
    end

    %% Update the covariance matrix used for ranking distance
    if params.updateCovariance
      n = length(info.nbrList);
      % Don't let radical decreases in number change C
      if (n > nOld/2 && ~(strcmp(params.covarianceModel,'full') && n <= d))
        C0 = Cnew;
        nOld = n;
      end
    end
  end
  return
  
  
  
  alpha_history = zeros(0,2);
  nOld = 0;
  info = struct('nbrList',[],'T2',[],'xpre',b,'xpost',[],'Cpre',C0,'interiorIndex',[],'n_iter',0);
  n_iter = 0;

  while true
    %% Shift the base point to include more data points (if possible)
    n_iter = n_iter+1;
    [sortOrder,nvec,mu,C] = cn_moments_by_distance(x,params,b,C0);
    [deltax,info] = cn_baseshift(nvec,mu,C,params);
    info.xpre = b;
    info.Cpre = C0;
    if ~isempty(deltax)
      b = b + deltax;
      info.nbrList = sortOrder(1:info.alphaindex);
    else
      b = b + mu(:,end);
      info.nbrList = sortOrder;
    end
    info.xpost = b;
    info.n_iter = n_iter;
    info

    %% User interaction
    if haveFunc
      params = params.func(b,info.nbrList,params);
    end
    
    %% Check for convergence
    if (isempty(deltax) || info.alpha(1) < params.alphaeps)
      return
    end
    % Rarely, one can get into a cycle---check for this too
%     if (info.n == info.alphaindex)  % are we reaching steady-state?
      if any(all(abs(bsxfun(@minus,alpha_history,info.alpha)) < params.alphaeps,2),1)
        return
      end
      alpha_history(end+1,:) = info.alpha;  %#ok<AGROW>   % add to history
%     end
    
    %% Update the covariance matrix used for ranking distance
    if params.updateCovariance
      n = info.alphaindex;
      % Don't let radical decreases in number change C
      if (n > nOld/2 && ~(strcmp(params.covarianceModel,'full') && n <= d))
        switch params.covarianceModel
          case {'isotropic','diagonal'}
            C0 = C(:,n);
          case 'full'
            C0 = C(:,:,n);
        end
      end
      nOld = n;
    end
  end
end
