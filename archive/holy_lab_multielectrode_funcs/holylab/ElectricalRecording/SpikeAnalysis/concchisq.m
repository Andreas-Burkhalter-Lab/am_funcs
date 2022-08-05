function [chisq,gradchisq,invcov] = concchisq(X,params)
% CONCCHISQ: compute the fitting error for steady-state firing rates
% This is an "inner loop" function for CONCCOMPARE for use in
% minimization by FMINUNC.
%
% [chisq,gradchisq,invcov] = concchisq(X,params)
% where
%   X(1) = rmax, X(2) = log(K) for best stimulus, X(3 ...) =
%     log(relconc) for the other stimuli, and optionally X(end) = r0
%     (baseline).
%   params is a structure with the following fields:
%     stimlabel (see help for CONCCOMPARE)
%     stimconc (")
%     rates (")
%     rateerrs(")
%     Xlabel: a cell array with unique stimulus labels, Xlabel{1} is the
%       name of the best stimulus, Xlabel{2:end} gives the names of the
%       other stimuli in X(3:end)
%     variablebaseline: if true, allows a variable baseline for the fit;
% and
%   chisq is the measure of goodness-of-fit;
%   gradchisq is its derivative with respect to the parameters X;
%   invcov is the inverse of the covariance matrix.
%
% See also: CONCCOMPARE.
  rmax = X(1);
  K = exp(X(2));
  if (K == 0)
      chisq = Inf;
      'K = 0 trouble...'
      return;
  end
  r0 = 0;
  if (isfield(params,'variablebaseline') & params.variablebaseline)
    r0 = X(end);
  end
  % Now build cOverK for each stimulus
  cOverK = params.stimconc/K;
  for i = 2:length(params.Xlabel)
    indx = strmatch(params.Xlabel{i},params.stimlabel,'exact');
    cOverK(indx) = cOverK(indx)*exp(X(1+i));
  end
  % Compute theoretical firing rates
  saturated = (cOverK./(1+cOverK))';
  r = rmax * saturated + r0;
  % Compute chisq
  err = r - params.rates;
  chisq = sum(err.^2./params.rateerrs.^2);
  if isinf(chisq) || isnan(chisq)
    warning('chisq is infinite or nan'); % put a breakpoint here
  end
  % Compute gradient & invcov, if desired
  if (nargout > 1)
    dr = zeros(length(r),length(X));
    dr(:,1) = saturated;
    dfac = r ./(1+cOverK');
    dr(:,2) = -dfac;
    for i = 2:length(params.Xlabel)
      indx = strmatch(params.Xlabel{i},params.stimlabel,'exact');
      tmp = zeros(size(dfac));
      tmp(indx) = dfac(indx);
      dr(:,1+i) = tmp;
    end
    if (isfield(params,'variablebaseline') & params.variablebaseline)
      dr(:,end) = 1;
    end
    gradchisq = 2*(dr'*(err ./ params.rateerrs.^2));
    if (nargout > 2)
      diagerrs = diag(1./params.rateerrs.^2);
      invcov = dr'*diagerrs*dr;
    end
  end
