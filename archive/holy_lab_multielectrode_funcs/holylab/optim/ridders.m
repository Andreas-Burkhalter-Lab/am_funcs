function t = ridders(func,tl,tr)
% RIDDERS: bracketed one-dimensional root finding
% Syntax:
%   t = ridders(func,tl,tr)
% where
%   func is a functional handle of a one-dimensional function
%   tl is to the left of a root
%   tr is to the right of a root (in particular, func(tl)*func(tr) must be
%     negative)
% and
%   t is the root.
%
% This is an alternative to the fzero function,one that allows you to
% control the interval in which the root is found.
%
% tl and tr can be vectors, in which case a vector of roots is also found.
% Your function must support elementwise evaluation when supplied with a
% vector.
%
% See also: FZERO.

% Copyright 2008 by Timothy E. Holy

  fl = func(tl);
  fr = func(tr);
  if any(fl .* fr > 0)
    error('Root must be bracketed');
  end
  t = Inf;
  dt = Inf;
  iter = 0;
  maxiter = 50;
  dtthresh = 1e-8;
  still_fixing = true(size(tl));
  while(any(still_fixing) && iter < maxiter)
    tm = (tl(still_fixing) + tr(still_fixing))/2; % bisection step
    fm = func(tm);
    fltmp = fl(still_fixing);
    frtmp = fr(still_fixing);
    sq = sqrt(fm.^2 - fltmp.*frtmp);
    isnz = (sq ~= 0);
    sc = zeros(size(fm));
    sc(isnz) = sign(fltmp(isnz) - frtmp(isnz)).*fm(isnz)./sq(isnz);
    tnew = tm + (tm-tl(still_fixing)) .* sc;
    index_fixing = find(still_fixing);
    fnew = func(tnew);
    index_done = index_fixing(fnew == 0);
    still_fixing(index_done) = false;
    tl(index_done) = tnew(index_done);
    tr(index_done) = tnew(index_done);
    
    % Move the bracketing points
    flag1 = fm .* fnew < 0;
    flag2 = fl .* fnew < 0;
    flag3 = fr .* fnew < 0;
    flag2 = flag2 | flag1;
    flag3 = flag3 & ~flag2;
    
    tl(index_fixing(flag1)) = tm(flag1);
    fl(index_fixing(flag1)) = fm(flag1);
    tr(index_fixing(flag2)) = tnew(flag2);
    fr(index_fixing(flag2)) = fnew(flag2);
    tl(index_fixing(flag3)) = tnew(flag3);
    fl(index_fixing(flag3)) = fnew(flag3);
    
    % Convergence criteria
    still_fixing(index_fixing(~isnz)) = false;
    still_fixing(abs(tr-tl) < dtthresh) = false;
    iter = iter+1;
  end
  t = (tl+tr)/2;
  