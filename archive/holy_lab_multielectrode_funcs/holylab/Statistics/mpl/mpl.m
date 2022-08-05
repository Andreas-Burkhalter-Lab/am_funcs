function [p,action,im] = mpl(x,n,W,options)
% MPL: calculate the maximum penalized likelihood density.
%
% The discrete action is minimized, and then the density at the supplied
% points is calculated.
%
% Syntax:
%   p = mpl(x,n,W,options)
%   [p,action] = mpl(x,n,W,options)
%   [p,action,chi2mtrx] = mpl(x,n,W,options)
% where
%   x contains the data points, and is of size d-by-Nunique;
%   n contains the multiplicity of each point, and is of size
%     1-by-Nunique;
%   W is the smoothing factor, where the projected data points W*x are
%     used with standard-width kernels;
% and
%   p is a vector containing the density at each point in x, or in
%     options.xeval as described below;
%   action is a 2-vector, [S SCV], where the 2nd is the cross-validated
%     action. This output is provided only when options.kernel_only (see
%     below) is false.
%   im is a structure of "intermediates," containing useful fields like
%     lambda, L (or G), psi, etc.
%
% You can adjust the results with the structure "options". options may
% contain the following fields:
%   xeval (default x): a vector containing the set of points at which the
%     density is to be evaluated (p(i) = density at xeval(i));
%     In d=1, the nu=1 MPL estimate is described in terms of its value at
%     the data points and exponentials between them.  To evaluate the
%     density at arbitrary locations, we simply find the two points on
%     either side (which might include +/- infinity) and then calculate
%     the exponentials. This gives us a fast algorithm.
%   kernel_only (default false): if true, calculates the kernel density
%     estimate rather than MPL.
%   outlier_frac (default 0): The fraction of points to be considered as
%     outliers, and eliminated from contributing to the cross-validated
%     action (see LOGPCVROBUST).
%
% See also: MPL_OPTW, MPL_OPT, LOGPCVROBUST.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'kernel_only')
    options.kernel_only = 0;
  end
  if (nargout > 1 && options.kernel_only)
    error('Multiple outputs provided only when options.kernel_only is false');
  end
  if ~isfield(options,'outlier_frac')
    options.outlier_frac = 0;
  end

  [d,M] = size(x);
  if (d > 1)
    % We have to use the slow algorithm. Particularly for xeval-type
    % inputs, this is actually easier to program, so let's do it first.
    X = W*x;
    N = sum(n);
    if options.kernel_only
      if isfield(options,'xeval')
        sd = sqrdist(X,W*options.xeval);
      else
        sd = sqrdist(X,X);
      end
      G = exp(-sd/2)/(2*pi)^(d/2);
      p = (n*G)*(abs(det(W))/N);
    else
      % MPL
      sd = sqrdist(X,X);
      G = exp(-sd/2)/(2*pi)^(d/2);
      V = exp(-sd/4)/(4*pi)^(d/2);
      ops.use_a = 1;
      [a,lambda] = mpl_initial_guess(G,V,n,ops);
      [a,lambda] = mpl_opt(a,lambda,G,V,n,ops);
      if (nargout > 1)
        psi = G*a;
        action0 = lambda*(a'*psi) - lambda;
        action(1) = action0 - 2*(n*log(psi)) - N*log(abs(det(W)));
        Gcv = G;
        for i = 1:size(G,1)
          Gcv(i,i) = 0;
        end
        psicv = Gcv*a;
        [lPcv,isOK] = logpcvrobust(psicv.^2,options.outlier_frac);
        action(2) = action0 - n*lPcv  - sum(isOK)*log(abs(det(W)));;
        if (nargout > 2)
          im.lambda = lambda;
          im.G = G;
          im.V = V;
          im.psi = psi;
          im.a = a;
        end
      end
      if isfield(options,'xeval')
        sd = sqrdist(X,W*options.xeval);
        G = exp(-sd/2)/(2*pi)^(d/2);
      end
      psi = a'*G;
      p = psi.^2*abs(det(W));
    end
    return
  end
  % OK, this is all in d=1
  [X,sort_order] = sort(W*x);
  dX = diff(X);
  L1 = calcL1d1(dX);
  N = sum(n);
  if options.kernel_only
    %u = (L1\n')/N;  % Here u is P
    u = tridiag_inv(L1{1},L1{2},L1{3},n'/N);
  else
    V = calcVd1(dX);
    ops.use_a = 0;
    [u,lambda] = mpl_initial_guess(L1,V,n,ops);  % here u is psi
    [u,lambda] = mpl_opt(u,lambda,L1,V,n,ops);
    if (nargout > 1)
      Lpsi = (n'/lambda)./u;
      action0 = lambda*(u'*Lpsi) - lambda;
      action(1) = action0 - 2*(n*log(u)) - N*log(abs(W));
      psicv = crossval1d(u',dX)';%u - Lpsi/2;
      [lPcv,isOK] = logpcvrobust(psicv.^2,options.outlier_frac);
      action(2) = action0 - (n*lPcv) - sum(isOK)*log(abs(W));
      if (nargout > 2)
        im.lambda = lambda;
        im.L = L1;
        im.V = V;
        im.psi = u;
      end
    end
  end
  if ~isfield(options,'xeval')
    [sos,sort_invert] = sort(sort_order);
    % Just evaluate at x
    if options.kernel_only
      p = abs(W)*u(sort_invert);
    else
      p = abs(W)*u(sort_invert).^2;
    end
    return
  end
  % OK, we have to evaluate at xeval
  p = zeros(1,length(options.xeval));  % Initially, this will be either P or psi
                                       % depending on kernel_only
  [Xeval,sort_order1] = sort(W*options.xeval);
  sdX = sinh(dX);
  expdX = exp(dX);
  % Evaluate for the points less than the first data point
  ie = 0;
  while (Xeval(ie+1) <= X(1))
    ie = ie+1;
  end
  if (ie > 0)
    p(1:ie) = u(1)*exp(Xeval(1:ie)-X(1));
  end
  % Evaluate for all the middle points
  ix = 2;
  while (ie <= length(Xeval) && ix <= length(X))
    ie_start = ie+1;
    while (ie+1 < length(Xeval) && Xeval(ie+1) <= X(ix))
      ie = ie+1;
    end
    if (ie >= ie_start)
      p(ie_start:ie) = ...
          (expdX(ix-1)*u(ix-1)-u(ix))/(2*sdX(ix-1))*exp(-(Xeval(ie_start:ie)-X(ix-1))) + ...
          (expdX(ix-1)*u(ix)-u(ix-1))/(2*sdX(ix-1))*exp((Xeval(ie_start:ie)-X(ix)));
    end
    ix = ix+1;
  end
  % Now get the points that are beyond the last X
  ie_start = ie+1;
  if (ie_start <= length(Xeval))
    p(ie_start:length(Xeval)) = u(end)*exp(X(end)-Xeval(ie_start:end));
  end
  % Now convert to real probability
  if options.kernel_only
    p = p*abs(W);     % p was P
  else
    p = p.^2*abs(W);  % p was psi
  end
  [so1s,sort_invert] = sort(sort_order1);
  p = p(sort_invert);
