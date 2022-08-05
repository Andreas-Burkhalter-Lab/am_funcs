function [f,finv] = mz2indx(mz,mzmin,options)
% MZ2INDX: convert m/z values to an index
% Syntax:
%   f = mz2indx(mz,mzmin)
% where
%   mz is a set of m/z values acquired during a single scan (sorted in
%     increasing order)
%   mzmin is the minimum m/z value
% and
%   f is a function handle such that f(mzprime) returns an index value
%     for each element of mzprime. Note that the index value is not an
%     integer; this allows you to do linear interpolation.
%
%  [f,finv] = mz2indx(mz,mzmin)
% also returns an inverse function, for converting index values back to
% m/z values.
%
% This function is useful for generating full-matrix versions of scans at
% maximum resolution.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'mode','matching');

  dmz = medfilt(diff(mz),3);
  switch options.mode
   case 'linear'
    % Resolution is constant, so the index is linearly related to mz
    dmz = dmz(end);
    f = @(mz) (mz - mzmin)/dmz;
    finv = @(mzi) dmz*mzi + mzmin;
   case 'matching'
    % In practice, the resolution is approximately logarithmic (meaning
    % deltamz is approximately linear in mz), but it's not perfect, so we
    % fit deltamz vs. mz with a quadratic to account for the slight
    % curvature.
    [p,S,mu] = polyfit(mz(1:end-1),dmz,2);
    xmin = (mzmin - mu(1))/mu(2);
    f = @(mz) mz2indx_fun(mz,p,mu,xmin);
    if (nargout > 1)
      finv = @(i) indx2mz_fun(i,p,mu,xmin);
    end
   otherwise
    error(['Mode ' options.mode ' not recognized']);
  end
end

function indx = mz2indx_fun(mz,p,mu,xmin)
  x = (mz - mu(1))/mu(2);
  z = p(2) + 2*p(1)*x;
  zmin = p(2) + 2*p(1)*xmin;
  delta = 4*p(1)*p(3) - p(2)^2;
  if (delta > 0)
    sdelta = sqrt(delta);
    indx = 2*(atan(z/sdelta) - atan(zmin/sdelta))/sdelta;
  elseif (delta == 0)
    indx = -2./z + 2/zmin;
  else
    sdelta = sqrt(-delta);
    indx = -2*(atanh(z/sdelta) - atanh(zmin/sdelta))/sdelta;
  end
  indx = mu(2)*indx+1;  % +1 is unit offset
end

function mz = indx2mz_fun(indx,p,mu,xmin)
  indx = (indx-1)/mu(2);  % undo the unit offset
  zmin = p(2) + 2*p(1)*xmin;
  delta = 4*p(1)*p(3) - p(2)^2;
  if (delta > 0)
    sdelta = sqrt(delta);
    z = sdelta*indx/2 + atan(zmin/sdelta);
    z = sdelta*tan(z);
  elseif (delta == 0)
    z = -2./(indx-2/zmin);
  else
    sdelta = sqrt(-delta);
    z = -sdelta*indx/2 + atanh(zmin/sdelta);
    z = sdelta*tanh(z);
  end
  x = (real(z)-p(2))/(2*p(1));  % imaginary part, if any, is roundoff
  mz = mu(2)*x+mu(1);
end