function [g,psig,sqrtdetJ,err,params] = register_multigrid_relax(g,psiM_0,lambda,covariantflag,params)
  if covariantflag
    warp_options = struct('covariant',true,'sqrt',true);
  else
    warp_options = struct('covariant',false,'sqrt',false);
  end
  warp_options.sqrtdetJ0 = params.sqrtdetJ0;

  [psig,w,sqrtdetJ] = register_warp(psiM_0,g,warp_options);
  dimKeep = params.pixel_spacing > 0;
  if isempty(g{1})
    sz = size(params.psiF);
    g = register_g0(sz(dimKeep));
  end
  gradE = register_gradEreg(g,params.psiF,psiM_0,psig,w,params.pixel_spacing(dimKeep),lambda,sqrtdetJ,covariantflag);
  % Define the penalty function (use a nested function, so intermediates
  % can be read out again)
  fmu = @(mu) register_E_wrapper(mu);
  % Step in the direction of the negative gradient, but only as far as
  % the function value decreases
  %if isempty(params.mu)
  %  [mu,err] = linmin(fmu,params);
  %  fprintf('linmin: mu = %g\n',mu);
  %else
  [mu,err] = lindec(fmu,params);
  %end
  % Recover the variables that we want to output
  params.mu = 2*mu;  % double it so the step size grows after successful steps
  params.startval = err; % save for the next iteration
  g = g_cur;

  function err = register_E_wrapper(mu)
    g_cur = g;
    for dimIndex = 1:length(g)
      g_cur{dimIndex} = g_cur{dimIndex} - mu*gradE{dimIndex};
    end
    fprintf('mu = %g\n',mu);
    try
      [psig,w,sqrtdetJ] = register_warp(psiM_0,g_cur,warp_options);
    catch
      err = Inf;
      return
    end
    err = register_E(params.psiF,psig,w,sqrtdetJ,params.pixel_spacing,lambda);
  end
  
end
