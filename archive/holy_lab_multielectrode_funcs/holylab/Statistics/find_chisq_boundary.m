function limit = find_chisq_boundary(func,x,coord,dx,deltachi2)
% FIND_CHISQ_BOUNDARY: vary fitting parameter until chisq increases by 1
%
% This function is useful for finding the confidence intervals in
% parameters under circumstances where the linearity of the model cannot be
% guaranteed.
%
% Syntax:
%   limit = find_chisq_boundary(func,x,coordIndex,dx)
%   limit = find_chisq_boundary(func,x,coordIndex,dx,deltachi2)
% where
%   func is a handle to your chisq function
%   x is location of the minimum
%   coordIndex is the index of the coordinate that you want to vary
%   dx is the (signed) amount by which to start varying the coordinate
%   deltachi2 (default 1) is the amount by which chisq is allowed to
%     increase (this defines the confidence interval; choose 4 for 95% conf
%     intervals).
%
% As an example,
%    limit = find_chisq_boundary(myfunc,x,1,0.1)
% would find the rightmost boundary for coordinate #1, where the initial
% value tested will be x(1)+0.1. Correspondingly,
%    limit = find_chisq_boundary(myfunc,x,1,-0.1)
% would find the leftmost boundary in this same coordinate.

% Copyright 2007 by Timothy E. Holy

  if (nargin < 5)
    deltachi2 = 1;
  end

  chisq = func(x);
  dxvec = 0*x; dxvec(coord) = dx;
  chitest = chisq_optimize_var(func,x + dxvec,coord);
  thresh = deltachi2;
  iter_max = 10;
  
  % Try to bracket the solution
  iter = 0;
  while (chisq + thresh > chitest && iter < iter_max)
    dxvec = 2*dxvec;
    chitest = chisq_optimize_var(func,x + dxvec,coord);
    iter = iter+1;
  end
  
  if (iter == iter_max)
    limit = Inf*dx;
    return
  end
  
  testfunc = @(alpha) chisq_optimize_var(func,x + alpha*dxvec,coord) - chisq - thresh;
  alpha = fsolve(testfunc,0.75);
  
  dxvec = alpha*dxvec;
  limit = dxvec(coord);

  
function ret = chisq_var_fixed(func,zvar,zfixed,isfixed)
% Calculate a function when "variable" and "fixed" coordinates have been
% specified separately
  y = zeros(1,length(isfixed));
  y(isfixed) = zfixed;
  y(~isfixed) = zvar;
  ret = func(y);

  
function [chisq,xout] = chisq_optimize_var(chisqfunc,x,coord)
% Hold one coordinate fixed and optimize the rest
  isfixed = 0*x; isfixed(coord) = 1; isfixed = logical(isfixed);
  optfunc = @(xvar) chisq_var_fixed(chisqfunc,xvar,x(coord),isfixed);
  minopts = optimset('fminsearch'); minopts.MaxFunEvals = 1000;
  [xvar,chisq] = fminsearch(optfunc,x(~isfixed),minopts);
  xout = x;
  xout(~isfixed) = xvar;

  
