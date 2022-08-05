function [mu,fmu,d] = lindip(f,mu,options)
% LINDIP: find a point in a trough, along a ray starting from 0
%
% This function is useful in descent minimization techniques, where you
% want to find an improvement in an objective function in the descent
% direction.  This function works along a ray starting from 0.  It seeks
% a point lower than the starting point, but also verifies that another
% step outward would take it back uphill, thus insuring that the chosen
% step is of reasonable size (inside the "dip").
%
% Syntax:
%   mu = lindip(f,mu)
%   [mu,fmu,d] = lindip(f,mu,options)
% where
%   f is a function, mu a scalar > 0, and f(mu) produces data.  In typical
%     examples, f(mu) is the value of the objective function at the
%     position mu along the ray. However, in cases where the objective
%     function cannot be readily represented by a single scalar, you can
%     have f(mu) return relevant data, and then provide a "compare_func" in
%     the options (see below).
%   mu is the starting stepsize along the 1-d ray
%   options can take the following fields:
%     muinc (default 10): the factor by which to increase mu in stepping
%       outwards
%     mudec (default 11): the factor by which to decrease mu in stepping
%       inwards
%     itermax (default 20): the number of times to try to bracket a
%       minimum before giving up and returning mu=0
%     f0 (optional): if supplied, it skips calculating f(0) and uses the
%       supplied value instead (this can save time if function evaluation
%       is costly, and you've already had to compute it anyway)
%     compare_func (optional): if this function is supplied, then the
%       data returned by f(mu) are evaluated using the
%       compare_func. compare_func must take a single cell array as input,
%       and return a scalar value for each cell: specifically,
%              compare_func({f1,f2})
%       must return a 2-vector and
%              compare_func({f1,f2,f3})
%       must return a 3-vector.  See IMCOMPARE for an example of a
%       compare_func (use compare_func = @(imc) imcompare(imref,imc)).
% and
%   mu is the value in the dip
%   fmu is f(mu).
%   d is the value corresponding to fmu returned from the compare_func.
%
% See also: LINMIN, LINDEC, IMCOMPARE, GRADDESCENT.

% Copyright 2007 by Timothy E. Holy
  
if (nargin < 3)
  options = struct;
end
if isstruct(mu)
  options = mu;
  if isfield(options,'mu')
    mu = options.mu;
  else
    mu = 1;
  end
end
options = default(options,'muinc',10);
options = default(options,'mudec',11);
options = default(options,'itermax',20);

if isfield(options,'f0')
  f0 = options.f0;
else
  f0 = f(0);
end

if isfield(options,'compare_func')
  compare_func = options.compare_func;
else
  compare_func = @lb_compare;
end

fmu = f(mu);
d = compare_func({f0,fmu});
% Initialize with sentinels so we can do this with a single bit of logic
muvec = [0 mu mu];
fmuvec = {f0,fmu,fmu};
d(3) = d(2);

iter = 0;
% Keep trying until the middle point is smaller than either edge point
while (~(d(2) < d(1) && d(2) < d(3)) && iter < options.itermax)
  % OK, the middle point is not smaller. Parsing what to do has to consider
  % the possibility that one or more elements of d may be nan, or that some
  % may be identical.
  if (d(2) < d(1) && isnan(d(3)))
    % Middle mu yields improved results, with nan this is good enough
    break
  end
  if (d(2) == d(1) && d(3) > d(2))
    % Flat between first and 2nd, then rising---this would come up when we
    % pretty much start out at the minimum and shouldn't move. So we can
    % quit. Note, however, that if d(3) == d(2) as well, then instead it's
    % likely that we have not moved far enough for f to start changing, so
    % we definitely don't want to quit in that case (that's a "move
    % outward" case).
    break
  end
  movein = true; % by default, move closer---this handles many nan situations
  if (d(3) <= d(2) && d(2) <= d(1))
    % There's a downhill or flat slope from 0 out to the farthest point, so
    % keep moving outward to find where the curve turns around
    movein = false;
  end
  if movein
    % Make the next point closer to 0
    muvec = [0 [1/options.mudec 1]*muvec(2)];
    fmuvec = {f0,f(muvec(2)),fmuvec{2}};
  else
    % Make the next point farther than the extreme point
    muvec = [0 [1 options.muinc]*muvec(3)];
    fmuvec = {f0,fmuvec{3},f(muvec(3))};
  end
  d = compare_func(fmuvec);
  iter = iter+1;
end
if (iter >= options.itermax)
  warning('optim:noConvergence','Failed to bracket a minimum');
  mu = 0;
  fmu = f0;
  d = d(1);
else
  mu = muvec(2);
  fmu = fmuvec{2};
  d = d(2);
end

function d = lb_compare(v)
d = [v{:}];
