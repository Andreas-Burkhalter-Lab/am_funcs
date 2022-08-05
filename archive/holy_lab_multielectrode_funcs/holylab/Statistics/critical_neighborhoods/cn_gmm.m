function [n_all,mu_all,C_all] = cn_gmm(x,params)
  % cn_gmm: create a Gaussian mixture model using critical neighborhoods
  % Syntax:
  %   [n,mu,C] = cn_gmm(x,params)
  % where
  %   x is the d-by-N data matrix
  %   params is a structure whose fields are described in cn_flow
  % and
  %   n is a 1-by-k vector of amplitudes
  %   mu is a d-by-k matrix of centers
  %   C is a an array of covariances, whose size depends upon the
  %     covariance model chosen in params (isotropic: 1-by-1-by-k;
  %     diagonal: d-by-1-by-k; full: d-by-d-by-k).
  %
  % See also: gmm_optimize.

  % Copyright 2011 by Timothy E. Holy
  
  %% Parse arguments
  params = default(params,'nMin',ceil(-log(params.pvalue)),'covarianceModel','isotropic','updateCovariance',true);
  [d,N] = size(x);
  covarianceModelIndex = strmatch(params.covarianceModel,{'isotropic','diagonal','full'});
  if isempty(covarianceModelIndex)
    error('covarianceModel not recognized')
  else
    covarianceModelIndex = covarianceModelIndex-1;
  end
  
  %% Determine the threshold for statistical significance
  Nlist = 1:N;
  T2thresh = inf(size(Nlist));
  switch params.covarianceModel
    case 'isotropic'
      keepFlag = Nlist > 1;
      Nlist = Nlist(keepFlag);
      T2thresh(keepFlag) = ((d*Nlist-1)./(Nlist-1)).*finv(1-params.pvalue,d,d*(Nlist-1));
    case 'diagonal'
      % Use the 2-moment approximation
      n2 = ((d+2)*Nlist-5*d+2)/3;
      s = ((Nlist-3).*n2)./(d*(Nlist-1).*(n2-2));
      keepFlag = Nlist >= 5;
      T2thresh(keepFlag) = finv(1-params.pvalue,d,n2(keepFlag))./s(keepFlag);
    case 'full'
      keepFlag = Nlist > d;
      Nlist = Nlist(keepFlag);
      T2thresh(keepFlag) = (d*(Nlist-1)./(Nlist-d)).*finv(1-params.pvalue,d,Nlist-d);
    otherwise
      error('Covariance model not recognized');
  end
  
  %% Iteratively add Gaussians
  ez2 = zeros(1,N);  % will hold the summed e^(-z^2) from all Gaussians
  psum = zeros(1,N); % will hold the summed likelihood from all Gaussians
  n_all = [];
  mu_all = {};
  C_all = {};
%   while (min(ez2) < params.pvalue)
  while (min(ez2) < 0.1)
    % Choose a random point in proportion to 1/p
    psel = ez2;
    psel(psel < params.pvalue^2) = params.pvalue^2;
    psel = 1./psel;
    indx = pick_randomly_weighted(psel);
    fprintf('starting index %d\n',indx)
    x0 = x(:,indx);
    ez2(indx) = inf;
    % Flow to peak
    [nbrList,triggerSource,x0flow,Cflow] = cn_flow_gmm(x,psum,params,T2thresh,x0);
    % Update the weights
    dx = bsxfun(@minus,x,x0flow);
    z2 = dist_mahalanobis(dx,Cflow,params.covarianceModel);
    [~,sortOrder] = sort(z2);
    pnew = exp(-z2/2)/sqrt((2*pi)^d * det_cov(Cflow,d,params.covarianceModel));
    w = pnew./(pnew+psum);
    w(pnew == 0) = 0;  % avoid NaNs
    % Determine the right number to be in the neighborhood, by determining
    % the peak of the statistic for just the first n points on the next
    % iteration
%     if (triggerSource > 0)
      T2N = chisqNP(dx,w,sortOrder(1:min(length(nbrList),ceil(sum(w)))),covarianceModelIndex);
      T2N(1:params.nMin-1) = 0;
      T2N(T2N < T2thresh(1:length(T2N))) = 0;
      [~,n] = max(T2N);
      if (n > 1)
        % If n == 1, then all T2 were 0; the data set contains a single peak
        nbrList = nbrList(1:n);
      end
%     end
    % Calculate the moments of this neighborhood
    xsel = x(:,nbrList);
    wsel = w(:,nbrList);
    [n,mu,C] = moments_from_points(xsel,params.covarianceModel,wsel);
    % Correct the parameters
    dx = bsxfun(@minus,xsel,mu);
    z2 = dist_mahalanobis(dx,C,params.covarianceModel);
    if (triggerSource ~= 0)
      [n,covfac] = correct_gaussian_parameters(d,n,max(z2)/(sum(z2.*wsel)/n));
      C = C*covfac;
    end
    % Update the probabilities of all points
    dx = bsxfun(@minus,x,mu);
    z2 = dist_mahalanobis(dx,C,params.covarianceModel);
%     this_ez2 = exp(-z2/2);
    this_ez2 = 1 - gammainc(z2/2,d/2).^(n+1);
    ez2 = ez2 + this_ez2;
    pnew = this_ez2/sqrt((2*pi)^d * det_cov(C,d,params.covarianceModel));
    psum = psum+pnew;
    % Store the parameters
    n_all(end+1) = n;
    mu_all{end+1} = mu;
    C_all{end+1} = C;
  end
  
  mu_all = cat(2,mu_all{:});
  C_all = cat(3,C_all{:});
end


    