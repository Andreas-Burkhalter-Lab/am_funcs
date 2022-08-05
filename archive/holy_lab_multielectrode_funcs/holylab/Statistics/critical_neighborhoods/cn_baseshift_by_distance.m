function [sortOrder,alpha,b,C,alphamax] = cn_baseshift_by_distance(x,params,b,C0)
%   b = basestruct.b + basestruct.alpha*basestruct.deltab;
%   [sortOrder,n_all,mu_all,C_all] = cn_moments_by_distance(x,params,b,basestruct.C);
  [sortOrder,n_all,mu_all,C_all] = cn_moments_by_distance(x,params,b,C0);
  T2 = cn_T2(n_all,mu_all,C_all,params.covarianceModel,params.T2thresh);
  n = find(T2 > params.T2thresh,1,'first');
  if isempty(n)
    % Nothing violated significance, let's quit
    alpha = 1;
    alphamax = inf(1,n);
    [~,b,C] = moments_from_points(x,params.covarianceModel);
    return
  end
  
  % Do the line search in the direction of the first violation 
  deltab = mu_all(:,n);
  if ~isinf(params.T2thresh(n-1))
    n = n-1;
  end
  sortOrder = sortOrder(1:n);
  params.positive_only = true;
  [alpha,alphamax,~,info.alphaindex] = cn_linesearch(mu_all,C_all,deltab,params);
  if (alpha > 1)
    alpha = 1;
  end
  
  % Flow to the 1-d peak, or alpha = 1
  xsel = x(:,sortOrder(1:n)); % we can look at just the points that contributed
%   xsel = x(:,sortOrder); % look at all points
  params.positive_only = false;
  alpha_history = [];
  n_history = [];
  while true
    % Check for cycling
    indx = find(abs(alpha_history - alpha)<params.alphaeps,1,'last');
    if ~isempty(indx)
      nmax = max(n_history(indx:end));
      if (n == nmax)
        break
      end
    end
    % Add to history
    alpha_history(end+1) = alpha; %#ok<AGROW>
    n_history(end+1) = n; %#ok<AGROW>
    % Pick new point
    [sortNew,~,mu_tmp,C_tmp] = cn_moments_by_distance(xsel,params,b+alpha*deltab,C0);
    if all(C_tmp(:) == 0)
      alpha1 = 0;
      break
    end
    [alpha1,alphamax,~,n] = cn_linesearch(mu_tmp,C_tmp,deltab,params);
    if (abs(alpha1) < params.alphaeps || alpha+alpha1 < 0 || (alpha == 1 && alpha1>0))
      break
    else
      alpha = alpha+alpha1;
      if (alpha > 1)
        alpha = 1;
      end
    end
  end
  if abs(alpha1) < params.alphaeps
    alpha = alpha+alpha1;  % This avoids some problems with roundoff errors?
  end
  if (alpha < 0 && abs(alpha) > params.alphaeps) % || alpha+alpha1 < 0)
    error('Alpha is negative')
  end
  if (alpha > 1)
    alpha = 1;
  end
  
%   fprintf('n = %d, alpha = %g, nsteps = %d\n',n,alpha,length(alpha_history)+1)
  
  b = b+alpha*deltab;
  sortOrder = sortOrder(sortNew(1:n));
  xsel = x(:,sortOrder);
  if (nargout > 3)
    [~,~,C] = moments_from_points(xsel,params.covarianceModel);
  end
  