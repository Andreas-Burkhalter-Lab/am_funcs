function p = fit1storderss(c,dr,drerr,options)
% FIT1STORDERSS: fit steady-state responses to Michaelis-Menten function
% Syntax:
%  p = fit1storderss(c,dr,drerr)
%  p = fit1storderss(c,dr,drerr,options)
% where
%   c are the concentrations
%   dr and drerr are the mean and s.e.m. of deltar at each concentration
% and
%   p is a parameters structure. The fields labeled "range" are the 68%
%     confidence interval (delta chi^2 = 1) for each fitting parameter
%     individually. The field labeled chisq_flat is the value of chisq for
%     a model in which deltar is independent of concentration (and thus
%     provides a good way to test whether there was a response).
%
% The model used is a first-order process,
%
%    r = rmax/(1 + (K/c))
%
% Optionally, you can also fit a Hill exponent n,
%
%    r = rmax/(1 + (K/c)^n)
%
% (by default, n is fixed as 1; but options.fit_n = true will fit n).
%
% Instead of returning rmax, this function fits an alternative
% parametrization,
%
%    r = s/(1/K^n + 1/c^n)
%
% where s = rmax/K^n.  s is fit reliably even if K cannot be precisely
% determined (i.e., in cases where little or no saturation is observed).
% Note that for n=1, s has units rate/concentration, i.e., is the initial
% slope of the response curve.
%
% Options:
%    fit_n (default false): if true, allows n to be different from 1
%    find_ranges (default true): if true, find the 68% confidence interval for
%      each fitting parameter.
%
% See also: fit_hill_equation, RATE1STORDER, CONCCOMPARE.

% Copyright 2007 by Timothy E. Holy

  default_options('fit_n',false,'find_ranges',true);
  
  drmax = max(dr);
  if any(drerr == 0)
    % At least one point has zero error bars, so numerical fit is
    % meaningless (and probably the data are weird anyway)
    p = struct('s',NaN,'log10K',NaN','K',NaN,'n',NaN,'chisq',Inf,'chisq_flat',Inf,...
        'drbar',NaN,'drmax',drmax,'dof',length(c)-(2+options.fit_n));
    if options.find_ranges
      p.srange = [NaN NaN];
      p.lKrange = [NaN NaN];
      if options.fit_n
        p.nrange = [NaN NaN];
      end
    end
    return
  end

  K = c(find(dr < drmax/2,1,'last'));
  if isempty(K)
    K = max(c);
  end
  n = 1;
  s = drmax/K^n;
  
  x(1) = s;
  x(2) = log10(K);
  if options.fit_n
    x(3) = n;
  end
  
  chisqfunc = @(y) f1oss_err(y,c,dr,drerr);
  minopts = optimset('fminsearch'); minopts.MaxFunEvals = 1000;
  [x,chisq] = fminsearch(chisqfunc,x,minopts);
  
  p.s = x(1);
  p.log10K = x(2);
  p.K = 10^(x(2));
  if options.fit_n
    p.n = x(3);
  else
    p.n = 1;
  end
  
  p.chisq = chisq;
  % Compare to a model with no concentration-dependence
  drbar = sum(dr ./ drerr.^2) / sum(1./drerr.^2);
  p.chisq_flat = sum((dr - drbar).^2./drerr.^2);
  p.drbar = drbar;
  p.drmax = drmax;
  p.dof = length(c) - 2 - options.fit_n;

  % Now individually vary the parameters until chisq increases by 1
  % (68% confidence interval). We have to re-fit the remaining parameters
  % to do this properly---but all this work is done by find_chisq_boundary.
  if options.find_ranges
    dxs = abs(0.1*x(1));
    if (dxs == 0)
      dxs = 1;
    end
    dxk = abs(0.1*x(2));
    if (dxk == 0)
      dxk = 1;
    end
    p.srange(1) = find_chisq_boundary(chisqfunc,x,1,-dxs);
    p.srange(2) = find_chisq_boundary(chisqfunc,x,1,dxs);
    p.lKrange(1) = find_chisq_boundary(chisqfunc,x,2,-dxk);
    p.lKrange(2) = find_chisq_boundary(chisqfunc,x,2,dxk);
    if options.fit_n
      p.nrange(1) = find_chisq_boundary(chisqfunc,x,3,-0.1*p.n);
      p.nrange(2) = find_chisq_boundary(chisqfunc,x,3,0.1*p.n);
    end
  end
  
  
function chisq = f1oss_err(y,c,dr,drerr)
% Calculate the fitting error as a function of parameters, data
  s = y(1);
  K = 10^(y(2));
  if (K == 0)
    chisq = Inf;
    return
  end
  if (length(y) > 2)
    n = y(3);
  else
    n = 1;
  end
  
  drth = s ./ (1/K^n + 1./c.^n);
  %semilogx(c,[dr(:) drth(:)]); pause
  
  chisq = sum((dr(:) - drth(:)).^2./drerr(:).^2);
  
  
  
  