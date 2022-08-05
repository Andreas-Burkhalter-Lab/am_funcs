function [g,imM,params] = register_multigrid_minimize(g,imM_0, ...
						  lambda,covariantflag,params,maxIter)
  % This is a dumb minimizer! Make this smarter!
  if (maxIter > 0)
    warp_ops = struct('covariant',covariantflag);
    psiM_0 = register_im2psi(imM_0,warp_ops);
    for i = 1:maxIter
      [g,psig,sqrtdetJ,err,params] = register_multigrid_relax(g,psiM_0, ...
        lambda,covariantflag,params);
      fprintf('%g (%g) ',err,params.mu);
      if (err == 0)
        break;
      end
    end
    fprintf('\n');
    imM = register_psi2im(psig,warp_ops);
  else
    imM = imM_0;
    if isempty(g{1})
      dimKeep = params.pixel_spacing > 0;
      sz = size(params.psiF);
      g = register_g0(sz(dimKeep));
    end
  end