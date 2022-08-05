function [pout,chisq] = fit_gaussian(x,y,varargin)
% FIT_GAUSSIAN: fit data to a gaussian w/ baseline
% Syntax:
%   p = fit_gaussian(x,y)
%   p = fit_gaussian(x,y,yerr)
%   p = fit_gaussian(...,p0)
%   [p,sse] = fit_gaussian(...)
% where
%   x is the set of locations at which observations are taken;
%   y is a set of observations which you desire to fit to a gaussian;
%   yerr (optional) is the set of error bars associated with each
%     observation;
%   p0 (optional) is an initial guess for the parameters, with the fields:
%     baseline: the baseline value for the data
%     center: the position of the center of the gaussian
%     peak: the height of the peak of the gaussian
%     sigma: the width of the gaussian ( exp(-x^2/(2*sigma^2)) )
% and
%   p is the optimal set of parameters (with the structure given above)
%   sse is the sum of square errors (or chisq, if you specify yerr) for the
%     fit

% Copyright 2007 by Timothy E. Holy

  y = y(:);
  x = x(:);
  pin = [];
  yerr = [];
  for i = 1:length(varargin)
    if isnumeric(varargin{i})
      yerr = varargin{i}(:);
    elseif isstruct(varargin{i})
      pin = varargin{i};
    end
  end
  if (length(x) ~= length(y))
    error('The length of x and y do not match');
  end
  if ~isstruct(pin)
    pin.baseline = min(y);
    pin.peak = max(y) - pin.baseline;
    dy = y - pin.baseline;
    pin.center = sum(dy.*x)/sum(dy);
    pin.sigma = sqrt(sum(dy.*x.^2)/sum(dy));
  end
  if isempty(yerr)
    yerr = ones(size(y));
  end
  if (length(yerr) ~= length(y))
    error('The length of y and yerr do not match');
  end
  
  p = [pin.baseline pin.center pin.peak pin.sigma]; % Convert to a vector
  fitfunc = @(t) fgchisq(t,x,y,yerr.^2);
  [p,chisq] = fminsearch(fitfunc,p);   % Do the fit
  pout = struct('baseline',p(1),'center',p(2),'peak',p(3),'sigma',p(4)); % Convert back to structure

function chisq = fgchisq(p,x,y,yerr2)
  %x = (1:length(y))';
  xsc = (x - p(2))/p(4);
  yf = p(1) + p(3)*exp(-xsc.^2/2);
  chisq = sum((y-yf).^2./yerr2);
  
