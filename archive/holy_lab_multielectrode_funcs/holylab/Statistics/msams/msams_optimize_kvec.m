function [kvec,step,w,converged,func_iters] = msams_optimize_kvec(x,x0,options)
% This may be recommended only for points that have converged as judged by
% regular MSAMS. Best way to call this is with R2 input from clustinfo.
  [d,N] = size(x);
  if ~isequal(size(x0),[d 1])
    error('x0 must be d-by-1');
  end
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'linked_coords',mat2cell(1:d,1,ones(1,d)),...
		    'thresh',3,'max_iter_bracket',20,'backtrace',false);
  options = default(options,'lb',1:length(options.linked_coords));
      
  ws = warning;
  warning('off','optim:fmincon:SwitchingToMediumScale');
      
  thresh = options.thresh^2; % Penalty is the squared ratio
  scalar_linked_coords = {1:d};
  dX = x - repmat(x0,1,N);
  
  func_iters = 0;
  if isfield(options,'R2')
    kscalar = 1./sqrt(options.R2);
    E = msams_vg_wrapper(kscalar,scalar_linked_coords);
    thresh = E;  % just keep whatever E we inherited
    fprintf('E0: %g\n',E);
  else
    %% Get a "scalar" value for k (i.e., homogenous scaling of coordinates)
    % Starting point: computed in terms of closest neighbors
    sd = sqrdist(x,x0);
    ssd = unique(sd);
    indx = 3;   % find 3rd closest neighbor
    if (ssd(1) == 0)
      indx = indx+1;
    end
    % Use this to construct an initial guess
    kscalar = 1/sqrt(ssd(indx));
    kscalarLast = Inf;
    % Bracket the root of E-thresh
    E = msams_vg_wrapper(kscalar,scalar_linked_coords);
    iter = 0;
    while (E > thresh)
      kscalarLast = kscalar;
      kscalar = 2*kscalar;
      E = msams_vg_wrapper(kscalar,scalar_linked_coords);
    end
    while (E < thresh && iter < options.max_iter_bracket)
      kscalarLast = kscalar;
      kscalar = kscalar/2;
      E = msams_vg_wrapper(kscalar,scalar_linked_coords);
      iter = iter+1;
    end
    if (iter >= options.max_iter_bracket)
      error('Failed to bracket kscalar');
    end
    % Find the root
    rootfun = @(ks) msams_vg_wrapper(ks,scalar_linked_coords) - thresh;
    kscalar = fzero(rootfun,[kscalar kscalarLast]);
  end
  
  %% Minimize prod(kvec) or sum(w) subject to E <= thresh
  %objfunc = @(kv) prod_gp(kv);
  objfunc = @(kv) sumw(kv,options.linked_coords);
  confunc = @(kv) mvg_fmc_wrapper(kv,options.linked_coords);
  fmcops = optimset('GradObj','on','GradConstr','on','TolX',1e-3*kscalar,'RelLineSrchBnd',0.1,'RelLineSrchBndDuration',Inf,'Display','iter');
  %fmcops = optimset('GradObj','on','GradConstr','on','TolX',1e-3*kscalar);
  kvec0 = kscalar * ones(1,length(options.linked_coords));
  [kvec,fval,exitstatus] = fmincon(objfunc,kvec0,[],[],[],[],options.lb,[],confunc,fmcops);
  converged = true;
  if (exitstatus < 1)
    [tmp,discrep] = confunc(kvec);
    warning('msams:failedConvergence','k-vector optimization failed to converge (achieved a value of %g, with constraint %g)',fval, discrep);
    converged = false;
  end
  
  %% Minimize prod(kvec) subject to E = 1.
  if options.backtrace
    % We only backtrace on the magnitude of kvec, we don't adjust the ratio
    % of the components.
    E = thresh;
    thresh = 1;
    alpha = 1;
    iter = 0;
    while (E > thresh && iter < options.max_iter_bracket)
      alphaLast = alpha;
      alpha = 2*alpha;
      E = mvg_scale_magnitude(alpha,options.linked_coords);
      iter = iter+1;
    end
    if (iter >= options.max_iter_bracket)
      error('Failed to bracket alpha during backtrace');
    end
    % Find the root
    rootfun = @(al) mvg_scale_magnitude(al,options.linked_coords) - thresh;
    alpha = fzero(rootfun,[alpha alphaLast]);
    kvec = alpha*kvec;
    %confunc = @(kv) mvg_fmc_wrapper(kv,options.linked_coords);
    %kvec = fmincon(objfunc,kvec,[],[],[],[],[],[],confunc,fmcops);
  end
  
  kvec = abs(kvec);
  warning(ws);
  
  %% Nested functions for computations
  function [E,gradE] = msams_vg_wrapper(kv,linked_coords)
    nlc = length(linked_coords);
    for i = 1:nlc
      kvfull(linked_coords{i}) = kv(i);
    end
    if (nargout == 1)
      [step,stepvar,w] = msams_val_gradient(x,x0,kvfull);
    else
      [step,stepvar,w,gradEfull] = msams_val_gradient(x,x0,kvfull);
    end
    E = (kvfull.^2 * step.^2) / (kvfull.^2 * stepvar);
    if (nargout > 1)
      gradE = zeros(nlc,1);
      for i = 1:nlc
        gradE(i) = sum(gradEfull(linked_coords{i}));
      end
    end
    func_iters = func_iters+1;
  end
  
  function [cineq,ceq,Gcineq,Gceq] = mvg_fmc_wrapper(kv,linked_coords)
    ceq = [];
    if (nargout > 2)
      Gceq = [];
      [cineq,Gcineq] = msams_vg_wrapper(kv,linked_coords);
      cineq = cineq-thresh;
    else
      cineq = msams_vg_wrapper(kv,linked_coords) - thresh;
    end
  end

  function E = mvg_scale_magnitude(alpha,linked_coords)
    kv = alpha*kvec;
    E = msams_vg_wrapper(kv,linked_coords);
  end

  function [val,grad] = prod_gp(kv)
    val = abs(prod(kv/kscalar));
    if (nargout > 1)
      grad = val ./ kv;
    end
  end
  
  function [val,grad] = sumw(kv,linked_coords)
    nlc = length(linked_coords);
    for i = 1:nlc
      kvfull(linked_coords{i}) = kv(i);
    end
    %kvfull

    [step,stepvar,w] = msams_val_gradient(x,x0,kvfull);
    val = -sum(w)/N;
    if (nargout > 1)
      gradfull = 2*kvfull(:) .* sum(repmat(w,d,1) .* dX.^2,2);
      grad = zeros(nlc,1);
      for i = 1:nlc
        grad(i) = sum(gradfull(linked_coords{i}))/N;
      end
    end
  end
end
