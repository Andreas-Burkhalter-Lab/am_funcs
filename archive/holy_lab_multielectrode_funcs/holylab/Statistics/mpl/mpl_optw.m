function [p,W] = mpl_optw(x,n,options)
% MPL_OPTW: optimize the smoothing length scale matrix W in MPL
%    (MPL = Maximum Penalized Likelihood density estimation)
% Syntax:
%   [p,W] = mpl_optw(x,n,options)
% where
%   x contains the data points, and is of size d-by-Nunique
%     (d is the dimensionality);
%   n contains the multiplicity of each point, and is of size 1-by-Nunique;
% and
%   p is the vector of estimated densities, p(i) is the density at x(:,i);
%   W is the optimized smoothing length scale matrix;
% You can control aspects of the computation with options, a structure
% which may have the following fields:
%     outlier_frac (default 0): helps with handling outliers, which can
%       distort the optimal W. Specifies the fraction of points to be
%       treated as outliers, and thus not contribute to the
%       cross-validated action. 0 retains all points; 0.05 would discard
%       the 5% least-likely points.
%     mode: 'full' allows W to be a full matrix (d-by-d),'diag' constrains
%       W to be a diagonal matrix, and 'isotropic' constrains W to be a
%       scalar times the identity matrix. This parameter obviously has no
%       effect in d=1.
%     start: starting guess for the optimum W. In d=1, this should actually
%       be a 2-vector [Wmin Wmax], with the minimum bracketed. Of course,
%       you can omit this field if you don't want to provide any starting
%       information.
%     kernel_only (default false): if true, just do kernel density
%       estimation, not true MPL.  For d > 1 and W full, it is recommended
%       that you first find W this way, because the gradient of the action
%       with respect to W can be easily computed (making minimization
%       faster) for kernel density estimation. Then you can polish up the
%       magnitudes using a call in 'isotropic' mode and full MPL.
%    optimization_fcn: for d>1, this determines which
%      optimization function to use.  The default, @fminunc, is in the
%      Optimization Toolbox, an add-on product. If it doesn't exist, then
%      one of two things will happen:
%        1. If 'mode' = 'full', this function will bail and instead call
%           PROJPURSCG. This is a custom-written conjugate-gradient solver
%           and thus uses derivative information. It is set up to handle
%           tricky issues relating to dimensionality reduction (which here
%           you don't care about, or you'd be calling PROJPURSCG
%           directly), so it makes it more complicated than a
%           run-of-the-mill CG solver. But it gets the job done, without
%           wandering around practically forever along the degenerate
%           directions (arising from the fact that an overall
%           rotation of the data has no impact), which is what would have
%           happened if you used a solver which couldn't use gradient
%           information.
%        2. Otherwise, it defaults to @fminsearch. Note that for
%           'isotropic' mode this is the default anyway.
%
%    display (default 'iter'): controls the amount of display in the matlab
%      minimization routines. See help for fminunc (d>1), fminbnd (d=1).
%    sorted (default false): if true and d=1, you're promising that the
%      data are already in sorted order. Bad Things will happen if this
%      isn't true.
%
% See also: PROJPURSCG, PPEK, LOGPCVROBUST.

% Copyright 2006 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  [d,Np] = size(x);
  if (d > Np)
    error(['x has more dimensions than points. Did you mean to supply the ' ...
           'transpose?']);
  end
  if (Np ~= length(n))
    error('The number of points does not match the number of multiplicities');
  end
  N = sum(n);
  if ~isfield(options,'kernel_only')
    options.kernel_only = 0;
  end
  if ~isfield(options,'mode')
    options.mode = 'full';
  end
  if isfield(options,'Display')
    options.display = options.Display;  % make case insensitive
  end
  if ~isfield(options,'display')
    options.display = 'iter';
  end
  if ~isfield(options,'outlier_frac')
    options.outlier_frac = 0;
  end
  if options.kernel_only
    options.gradobj = 'on';
  else
    options.gradobj = 'off';
    if (d > 1 && strcmp(options.mode,'full'))
      error('mpl_optw:mode_not_supported',...
            ['Optimizing a full multidimensional W with MPL is not supported.\nFirst call with options.kernel_only = 1, then project the data with W, and finally\ncall a second time with options.kernel_only = 0, options.mode = ''isotropic''.\nThe final W is the product of Wdiag*W.  See paper for details.']);
    end
  end
  % gradobj, test for full & mpl, etc
  if (d > 1 & ~isfield(options,'optimization_func'))
    if strcmp(options.mode,'isotropic')
      options.optimization_func = @fminsearch;
    elseif (exist('fminunc') == 2)
      options.optimization_func = @fminunc;
    elseif strcmp(options.mode,'full')
      warning('mpl_optw:optimizer','Need an optimizer which can use derivative information.\nIf you want to run this function, you need to supply one.\nI''m defaulting to projpurscg instead.')
      W = projpurscg(x,n,d,struct('sequential',0,'start','pca'));
      S = modwork('SCV',options,W,x,n);
      [p,lambda] = modwork('report',options);
      p = abs(det(W))*p;
      return
    else
      options.optimization_func = @fminsearch;
      gradobj = 'off';
    end
  end
  
  if (d == 1)
    % One dimension
    if ~isfield(options,'sorted')
      options.sorted = 0;
    end
    if ~options.sorted
      [x,sort_index] = sort(x);
      [ssi,sort_invert] = sort(sort_index);
    end
    dx = diff(x);
    if ~isfield(options,'start')
      options.start = 1./[diff(x([1 end])), median(dx)];
    end
    % Note SCV is more quadratic with respect to log(W) than it is to
    % W itself; hence, minimization proceeds somewhat better if
    % minimize with respect to log(W)
    ops = optimset('Display',options.display);
    lgw = fminbnd(@(lgw) mo1work('SCV',options,lgw,dx,n),...
      log(options.start(1)),log(options.start(2)),ops);
    W = exp(lgw);
    [psi,lambda] = mo1work('report',options);
    if ~options.sorted
      psi = psi(sort_invert);
    end
  else
    % Multidimensions
    if ~isfield(options,'start')
      % Starting guess by PCA. Since we have multiplicity information, we
      % need to do the PCA in a way that uses it.
      xm = mean(x,2);
      dx = x - repmat(xm,1,Np);
      nd = spdiags(n',0,Np,Np);  % multiplicity info
      Cx = (dx*nd*dx')/(N-1);
      % Note that if this is too high-dimensional, you're in trouble
      % anyway.
      % See ppek or projpurscg in such cases.
      [U,S,V] = svd(Cx);  % do PCA
      dS = diag(S);
      switch options.mode
        case 'full'
          W0 = diag(1./sqrt(dS))*V';
        case 'diag'
          W0 = dS;
        case 'isotropic'
          W0 = median(dS);
        otherwise
          error(['Mode ' options.mode ' not recognized.']);
      end
    else
      W0 = options.start;
    end
    minops = optimset('Display','iter','GradObj',options.gradobj);
    W = fminunc(@(W) modwork('SCV',options,W,x,n),W0(:),minops);
    W = reshape(W,size(W0));
    [psi,lambda] = modwork('report',options);
  end
  % Now convert to probability
  if options.kernel_only
    p = sqrt(det(W'*W))*psi;  % psi was actually P
  else
    p = sqrt(det(W'*W))*psi.^2;
  end
  
function [out,lambda_out] = mo1work(task,options,lgw,dx,n)
  persistent psi P lambda
  switch task
   case 'SCV'
    W = exp(lgw);
    dX = W*dx;
    L1 = calcL1d1(dX);
    N = sum(n);
    if options.kernel_only
      lambda = N;
      P = tridiag_inv(L1{1},L1{2},L1{3},n'/N);
      %Pcv = P - n'/(2*N);
      Pcv = crossval1d(P',dX);
      [lPcv,isOK] = logpcvrobust(Pcv,options.outlier_frac);
      out = -n*lPcv - sum(isOK)*lgw;  % The cross-validated log likelihood
    else
      V = calcVd1(dX);
      ops.use_a = 0;
      [psi,lambda] = mpl_initial_guess(L1,V,n,ops);
      [psi,lambda] = mpl_opt(psi,lambda,L1,V,n,ops);
      %Lpsi = L1*psi;
      Lpsi = (n'/lambda)./psi;
      %psicv = psi - Lpsi/2;
      psicv = crossval1d(psi',dX)';
      [lPcv,isOK] = logpcvrobust(psicv.^2,options.outlier_frac);
      out = lambda*(psi'*Lpsi) - (n*lPcv) - lambda - sum(isOK)*lgw;
    end
   case 'report'
    if options.kernel_only
      out = P;
    else
      out = psi;
    end
    lambda_out = lambda;
  end

function [out,outgrad] = modwork(task,options,W,x,n)
  persistent psi P lambda
  switch task
   case 'SCV'
    [d,M] = size(x);
    switch options.mode
     case 'full'
      W = reshape(W,length(W)/d,d);
      X = W*x;
     case 'diag'
      X = spdiags(W,0,d,d)*x;
     case 'isotropic'
      X = W*x;
     otherwise
      error(['Mode ' options.mode ' not recognized']);
    end
    sd = sqrdist(X,X);
    G = exp(-sd/2)/(2*pi)^(d/2);
    N = sum(n);
    if options.kernel_only
      lambda = N;
      P = (G*n')/N;
      %Pcv = P - G(1,1)*n'/N;
      Gcv = G;
      for i = 1:size(G,1)
        Gcv(i,i) = 0;
      end
      Pcv = (Gcv*n')/N;
      [lPcv,isOK] = logpcvrobust(Pcv,options.outlier_frac);
      Nok = sum(isOK);
      % The cross-validated log likelihood
      out = -n*lPcv; % will add -N*log(abs(det(W))) shortly
      % The gradient
      if ~strcmp(options.mode,'isotropic')
        outgrad = gradwllcv(x,X,n,G,Pcv,options.mode);
      else
        outgrad = (W/Nok) * (n./Pcv')*(sd.*G)*n';
      end
      switch options.mode
        case 'full'
          out = out - Nok*log(abs(det(W)));
          outgrad = outgrad - Nok*inv(W*W')*W;
        case 'diag'
          out = out - Nok*sum(log(abs(W)));
          outgrad = outgrad - Nok./W;
        case 'isotropic'
          out = out - Nok*d*log(abs(W));
          outgrad = outgrad - Nok*d/W;
      end
    else
      V = exp(-sd/4)/(4*pi)^(d/2);
      ops.use_a = 1;
      [a,lambda] = mpl_initial_guess(G,V,n,ops);
      [a,lambda] = mpl_opt(a,lambda,G,V,n,ops);
      psi = G*a;
      Gcv = G;
      for i = 1:size(G,1)
        Gcv(i,i) = 0;
      end
      psicv = Gcv*a;
      [lPcv,isOK] = logpcvrobust(psicv.^2,options.outlier_frac);
      Nok = sum(isOK);
      out = lambda*(psi'*a) -(n*lPcv) - lambda;
      switch options.mode
        case 'full'
          out = out - Nok*log(abs(det(W)));
        case 'diag'
          out = out - Nok*sum(log(abs(W)));
        case 'isotropic'
          out = out - Nok*log(abs(W))*d;
      end
    end
   case 'report'
    if options.kernel_only
      out = P;
    else
      out = psi;
    end
    outgrad = lambda;
  end
  