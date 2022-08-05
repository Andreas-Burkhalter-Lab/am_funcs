function [w,f1,f2,kappa] = bn_preordered_gaussian(dx,d2,options)
% BN_PREORDERED_GAUSSIAN: balanced neighborhood using Gaussian weights
% Syntax:
%   nbrInfo = bn_preordered_gaussian(dx,d2)
% The conventional choice for d2 is
%     d2 = sum(dx.^2,1)
% but it doesn't have to be that way (e.g., in filtering, it could be the
% square temporal displacement).

% Previously there was the ability to supply a user-defined guess, but this
% has been deemed too dangerous.
  
% Copyright 2009 by Timothy E. Holy
  
  [d,N] = size(dx);
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'plot',false,'min_to_check',10);
  
  id2 = 1./d2;
  id2(isinf(id2)) = 0;
  kappa = sqrt(mean(id2));
  %kappa = [];  % ignore the user-provided value, it's too risky
  %if (isempty(kappa) || isinf(kappa) || isnan(kappa))
  %  sd2 = sort(d2);
  %  n = find(sd2 > 0,options.min_to_check,'first');  % smallest possible neighborhood: mostly includes the closest point
  %  kappa = 1/sqrt(mean(sd2(n))); % should be an upper bound for kappa
  %end
  %kappa = abs(kappa);
  % If the neighborhood isn't big enough (criterion < 0), grow it
  kappalast = [];
  while (1)
    fval = bng_criterion(log(kappa));
    wsum = sum(w);
    if (fval < 0 && wsum < N-1)
      kappalast = kappa;
      kappa = kappa/1.4;
    else
      break
    end
  end
  if options.plot
    kapparange = logspace(-1,1,100)*kappa;
    f1plot = zeros(length(kapparange),1);
    f2plot = zeros(length(kapparange),1);
    ws = zeros(1,length(kapparange));
    for i = 1:length(kapparange)
      bng_criterion(log(kapparange(i)));
      f1plot(i) = f1;
      f2plot(i) = f2;
      ws(i) = sum(w);
    end
    hfig = figure;
    %hax = plotyy(kapparange,fval,kapparange,ws);
    %set(hax,'XScale','log')
    loglog(ws,[f1plot.^2 f2plot]);
  end
  if isempty(kappalast)
    % The neighborhood was too big to begin with, but it should be easy to
    % shrink
    lkappa = fzero(@(k) bng_criterion(k),log(kappa),optimset('Display','none'));
    kappa = exp(lkappa);
  else
    if (wsum > N-1)
      kappa = 0;  % the neighborhood includes all points
    else
      lkappa = fzero(@(k) bng_criterion(k),log([kappalast kappa]),optimset('Display','none'));
      kappa = exp(lkappa);
    end
  end
%   if options.plot
%     line([kappa kappa],[-1 1])
%   end


  % The following (nested) function returns zero when the balanced
  % neighborhood criterion is satisfied. It is parametrized in terms of
  % log(kappa) so that the entire interval [-inf inf] gets converted to [0
  % inf]; moreover, log-scaling seems to be more appropriate (i.e., yields
  % more quadratic-like functions)
  function z = bng_criterion(lk)
    factor = 1;
    k = exp(lk);
    w = exp((-k^2/2)*d2);
    wrep = repmat(w,d,1);
    f1 = sum(wrep.*dx,2);
    f2 = sum(wrep.*dx.^2,2);
    z = sum(f1.^2)/sum(f2) - factor;
  end
end
