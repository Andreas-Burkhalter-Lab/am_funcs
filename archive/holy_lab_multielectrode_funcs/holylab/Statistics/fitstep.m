function [xstep,yfit,chi2, deadloop] = fitstep(y,options)
% FITSTEP: fit a sequence of points to a step transition
% The data points y are fit to a model
%     ym = yfit(1) if x < xstep; ym = yfit(2) if x >= xstep
% where x is the index coordinate along y.  The fit minimizes chi-squared.
% Syntax:
%   [xstep,yfit,chi2] = fitstep(y)
% or
%   [xstep,yfit,chi2] = fitstep(y,options)
% where
%   y is a vector containing the observations;
%   options is a structure with the following fields:
%     algorithm, a string with values 'linear' or 'log', see
%       below (default 'linear');
%     xstart, if supplied gives a starting value for xstep (used only for
%       'linear');
%     guess, a string with values 'largest' or 'bestfit' which determnes
%       the starting position for the 'linear' method, see below (default
%       'bestfit');
%     report, if true causes it to print out the number of iterations
%       required (default 0);
% and
%   xstep is the transition point;
%   yfit is a 2 vector containing the mean values [before after] the
%     transition, respectively;
%   chi2 is the returned value of chisq.
% 
% The algorithm is basically a search: For 'log', it starts with xstep in
%   the middle (ignoring xstart), determines whether xstep should increase
%   or decrease, and then changes xstep by logarithmically-shrinking sizes.
%   For 'linear', it starts with xstep either at
%     1. xstart, if supplied by the user;
%     2. if guess = 'bestfit', at the point at which the transition in y
%        (meaning y(n+1)-y(n-1)) is closest in value to y(end)-y(1);
%     3. if guess = 'largest', at the point at which the largest transition
%        in y occurs.
%   From this point on, the algorithm iteratively changes xstep by 1 unit
%   in the direction which decreases the error.  
% Note that if the noise is high, you can fall into a local minimum. You
%   might want to median-filter your signal before calling this function.
% In general, if it can make a good starting guess, then 'linear' performs
%   best.
%
% Copyright Timothy E. Holy, 2004-07-22
  
  deadloop=0;

  n = length(y);
  if (n < 2)
    error('There must be at least 2 points to fit to a step');
  end
  if (nargin < 2)
    options = struct;
  end
  if ~isfield(options,'algorithm')
    options.algorithm = 'linear';
  end
  if ~isfield(options,'report')
    options.report = 0;
  end
  if ~isfield(options,'guess')
    options.guess = 'bestfit';
  end
  niter = 0;
    
  if strcmp(options.algorithm,'linear')
    % Pick the starting point for xstep
    if isfield(options,'xstart')
      xstep = options.xstart;
    elseif strcmp(options.guess,'largest')
      [mdy,xstep] = max(abs(diff(y))); % guess xstep at largest y transition
      xstep = xstep+1;
    else   % options.guess = 'bestfit'
      [mdy,xstep] = min(abs(y(3:end)-y(1:end-2) - (y(end)-y(1))));
      xstep = xstep+1;
    end
    % Slide in the better-fit direction
    [chi2,dchi2] = fitstep_error(y,xstep);
    bestchi2 = chi2;
    best_y=y;
    best_xstep=xstep;
    
    s = sign(dchi2); % This will help us keep track of whether we've
                     % changed movement direction
    tnIterations=0;                 
    while (s*dchi2 > 0 & xstep > 1 & xstep < n)
      niter = niter+1;
      xstep = xstep - s;
      [chi2,dchi2] = fitstep_error(y,xstep);
      if (chi2 < bestchi2)
        bestchi2 = chi2;
        best_y=y;
        best_xstep=xstep;
      else
        xstep = xstep + s;  % Oops, we went too far, go back a step (will
                            % exit on the next check)
      end
      tnIterations=tnIterations+1;
      if(tnIterations>10000) 
         % dbstop in fitstep at 99
         % chi2=bestchi2;
         % y=best_y;
         % xstep=best_xstep;
         % break;
         deadloop=1;
         yfit=[]; 
         return;
      end
    end
    if (xstep == 1 | xstep == n)
      error('Got to the edge of the interval');
    end
    yfit(1) = mean(y(1:xstep-1));
    yfit(2) = mean(y(xstep:end));
    chi2 = bestchi2;
    if options.report
      fprintf('fitstep niter %d, length(y) = %d\n',niter,n);
    end
    return;
  else % 'log'
    xstep = ceil(n/2);
    xjump = xstep/2;  % This doesn't have to be an integer
    [chi2,dchi2] = fitstep_error(y,xstep);
    bestchi2 = chi2;
    niter = 1;
    xold = xstep;
    xstep = round(xstep - sign(dchi2)*xjump);
    while (xjump >= 1 & xstep > 1 & xstep < n)
      [chi2,dchi2] = fitstep_error(y,xstep);
      if (chi2 < bestchi2)
        bestchi2 = chi2;
        xold = xstep;
      else
        xstep = xold;  % This didn't improve things, go back
      end
      niter = niter+1;
      xjump = xjump/2;
      xstep = round(xstep - sign(dchi2)*xjump);
    end
    if (xstep == 1 | xstep == n)
      error('Got to the edge of the interval');
    end
    xstep = xold;
    yfit(1) = mean(y(1:xstep-1));
    yfit(2) = mean(y(xstep:end));
    chi2 = bestchi2;
    if options.report
      fprintf('fitstep niter %d, length(y) = %d\n',niter,n);
    end
    return;
  end
  
  
function [chi2,dchi2,yfit] = fitstep_error(y,xstep)
  % Note this assumes that xstep is an integer!
  yfit(1) = mean(y(1:xstep-1));
  yfit(2) = mean(y(xstep:end));
  err = zeros(size(y));
  err(1:xstep-1) = y(1:xstep-1) - yfit(1);
  err(xstep:end) = y(xstep:end) - yfit(2);
  chi2 = sum(err.^2);
  % Calculate the change in chi2 if xstep->xstep+1
  dchi2 = (y(xstep)-yfit(1))^2 - (y(xstep)-yfit(2))^2;
  