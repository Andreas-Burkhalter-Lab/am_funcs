function [u,err,psig] = register_slow_multilevel(psi1,psi2,uOld,alpha)
% REGISTER_SLOW_MULTILEVEL: incremental warping for registration
% Syntax:
%   [u,err,psig] = register_slow_multilevel(psi1,psi2,u_guess,alpha)
% where
%   psi1 and psi2 are the input images (see REGISTER_SLOW);
%   u_guess is the hierarchy of warps (g - identity) that specifies your
%     initial guess;
%   alpha (default 1) controls the amount that you rely on your guess for
%     how warping should increase with finer levels of the hierarchy (1 =
%     use all the info, 0 = use none);
%  and
%    u is the output hierarchy of warps;
%    err is a vector of fitting errors, one per level of the hierarchy;
%    psig is the warped image, using the finest level of the hierarchy.
%
% See also: REGISTER_SLOW.

% Copyright 2006 by Timothy E. Holy
 
  n_g = length(uOld);
  n_dims = ndims(psi1);
  for gIndex = 1:n_g
    g_sz{gIndex} = size(uOld{gIndex}{1});
  end
  u0 = uOld{1};
  gIndex = 1;
  while (gIndex <= n_g)
    [utmp,err(gIndex),psig] = register_slow(psi1,psi2,u0);
    u{gIndex} = utmp;
    if (gIndex+1 <= n_g)
      % Update u0 for the next round
      gIdent_prev = register_g0(g_sz{gIndex});
      gIdent_next = register_g0(g_sz{gIndex+1});
      for dimIndex = 1:n_dims
        g0{dimIndex} = utmp{dimIndex} + gIdent_prev{dimIndex};
        g0Old{dimIndex} = uOld{gIndex}{dimIndex} + gIdent_prev{dimIndex};
        gOld{dimIndex} = uOld{gIndex+1}{dimIndex} + gIdent_next{dimIndex};
      end
      g0 = register_expandg(g0,g_sz{gIndex+1});
      g0Old = register_expandg(g0Old,g_sz{gIndex+1});
      % Enhance the initial guess by adding in the difference between the
      % optimized u and expanded u on the previous iteration---this
      % basically re-introduces the previous "extra" warp associated with
      % going to the next finer level.
      % alpha < 1 means we only add in a portion of this initial guess;
      % this might help "relax out" spurious values
      for dimIndex = 1:n_dims
        u0{dimIndex} = g0{dimIndex} + alpha * (gOld{dimIndex} - g0Old{dimIndex}) ...
          - gIdent_next{dimIndex};
        u0{dimIndex} = double(u0{dimIndex});
      end
    end
    gIndex = gIndex+1;
  end
