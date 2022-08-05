function [xo,eo] = trace_mirror_zernike(x,e,p)
% TRACE_MIRROR_ZERNIKE: trace rays reflecting from a mirror whose shape is specified by zernike coefficients
% Syntax:
%   [xo,eo] = trace_mirror_zernike(x,e,p)
% where
%   x is a 3-by-n_rays matrix containing a point on each ray. The plane of
%     the flat mirror is defined by x(3,:) = 0.
%   e is a 3-by-n_rays matrix, where each column is a unit vector giving
%     the propagation direction.
%   p is a parameter vector. p(1) is the "aperture size" for the mirror,
%     with the Zernike polynomials evaluated in terms of x(1:2,:)/p(1).
%     The remaining parameters are Zernike coefficients, starting at
%     second order. This function requires all coefficients of a given
%     order to be specified. So, if you wanted to use Zernike polynomials
%     up through 3rd order, you'd need to specify 3 coefficients for the 3
%     2nd order terms (Z(2,-2), Z(2,0), and Z(2,2)) and 4 coefficients for
%     the 4 3rd order terms (Z(3,-3), Z(3,-1), Z(3,1), Z(3,3)).
% and
%   xo is the set of strike points on the mirror surface for each ray (a
%     3-by-n_rays matrix)
%   eo is the matrix containing the output propagation directions.

% Copyright 2008 by Timothy E. Holy

  % Repack the parameters. The first parameter is aperture size, the
  % remaining are zernike coefficients of 2nd order and higher
  aperture = p(1);   % aperture size
  p(1) = [];
  % convert parameter vector into zernike orders
  zp = cell(1,0);
  zm = cell(1,0);
  zn = [];
  zf = cell(1,0);
  zfdx = cell(1,0);
  zfdy = cell(1,0);
  n = 2;
  while ~isempty(p)
    mmin = -n;
    m = mmin:2:-mmin;
    if (length(p) < length(m))
      error('All Zernike coefficients of a given order must be specified');
    end
    zp{end+1} = p(1:length(m));
    p(1:length(m)) = [];
    zm{end+1} = m;
    zn = n;
    cindx = length(zp);
    for i = 1:length(m)
      % Calculate the zernike polynomial and its derivative with respect to
      % x and y
      [zf{cindx}{i},zfdx{cindx}{i},zfdy{cindx}{i}] = zernikestringr(n,m(i));
    end
    n = n+1;
  end
  % Find the intersection of each ray and the mirror surface
  val0 = zern_val(zp,zf,x,aperture) - x(3,:);
  if any(val0 .* e(3,:) < 0)
    % The mirror shape is crashing into the previous component!
    xo = [];
    eo = [];
    return;
  end
  % Bracket the root between 0 and something big enough to pass through the
  % mirror surface
  t = -x(3,:)./e(3,:);  % Distance to flat mirror plane
  xc = x + repmat(t,3,1) .* e;
  val1 = zern_val(zp,zf,xc,aperture) - xc(3,:);
  istiny = abs(val1) < 1e-12;
  val1(istiny) = 0;
  is_not_bracketed = val1 .* val0 > 0;
  iter = 0;
  maxiter = 10;
  while any(is_not_bracketed) && iter < maxiter
    t(is_not_bracketed) = 2*(t(is_not_bracketed) + mean(t));
    xc = x + repmat(t,3,1) .* e;
    val1 = zern_val(zp,zf,xc,aperture) - xc(3,:);
    is_not_bracketed = val1 .* val0 > 0;
    iter = iter+1;
  end
  if (iter >= maxiter)
    error('Didn''t bracket');
  end
  % Find the root by Ridder's method
  tr = t;
  tl = zeros(size(t));
  t = tl;
  s = sign(val0-val1);
  if any(s ~= s(1))
    error('Unexpected sign mismatch')
  end
  s = s(1);
  dt = Inf;
  dtthresh = 1e-8;
  iter = 0;
  maxiter = 20;
  while (any(abs(dt) > dtthresh) && iter < maxiter)
    tm = (tl + tr)/2;
    xc = x + repmat(tm,3,1) .* e;
    valm = zern_val(zp,zf,xc,aperture) - xc(3,:);
    sq = sqrt(valm.^2 - val1.*val0);
    isz = sq == 0;
    sc = zeros(size(valm));
    sc(~isz) = valm(~isz)./sq(~isz);
    tnew = tl + s*(tm - tl).*sc;
    dt = tnew - tm;
    
  t2 = abs((val1 + val0)./(val1-val0)) .* t;
  %t2(istiny) = t(istiny);
  while (any(abs(dt) > dtthresh) && iter < maxiter)
    xc = x + repmat(t2,3,1) .* e;
    val2 = zern_val(zp,zf,xc,aperture) - xc(3,:);
    sq = sqrt(val2.^2 - val0.*val1);
    sc = val2 ./ sq;
    sc(val2 == 0 | val1 == 0 | sq == 0) = 0;
    t3 = t2 .* (1+s*sc);
    dt = t3-t2;
    t2 = t3;
    istoobig = t2 > tbracket;  % to deal with roundoff error problems
    if any(istoobig)
      warning('too big');
    end
    t2(istoobig) = tbracket(istoobig);
    iter = iter+1;
  end
  if (iter == maxiter)
    error('Ridder''s didn''t converge');
  end
%   
%     % Newton method for finding t
%     zs = zern_slope(zp,zfdx,zfdy,xc,aperture);
%     slope = sum(zs .* e);
%     slope = slope + sign(slope) * 1e-3;  % keep it from being too small
%     if any(slope == 0)
%       error('Oops!')
%     end
%     dt = -val./slope;
%     % Test and see whether this improves the fit
%     isworse = true;
%     while any(isworse) && iter < 100
%       tnew = t+dt;
%       xc = x + repmat(tnew,3,1) .* e;
%       valnew = zern_val(zp,zf,xc,aperture) - xc(3,:);
%       isworse = abs(valnew) > abs(val) + 100*eps;
%       if any(isworse)
%         dt(isworse) = dt(isworse)/(10 * iter);
%         iter = iter+1;
%       end
%     end
%     dt(isworse) = 0;
%     t = t + dt;
%     val = valnew;
%     iter = iter+1;
%   end
  % Calculate the local normal at each strike position
  zs = zern_slope(zp,zfdx,zfdy,xc,aperture);
  nlocal = zs ./ repmat(sum(zs.^2,1),3,1);
  % Use this to calculate the propagation direction
  edotn = sum(e.*nlocal,1);
  eo = e - 2*repmat(edotn,3,1).*nlocal;
  xo = xc;
end

function val = zern_val(zp,zf,xc,aperture)
  n = length(zp);
  val = zeros(1,size(xc,2));
  x = xc(1,:)/aperture;
  y = xc(2,:)/aperture;
  for i = 1:n
    for j = 1:length(zp{i})
      val = val + zp{i}(j) * eval(zf{i}{j});
    end
  end
end

function slope = zern_slope(zp,zfdx,zfdy,xc,aperture)
  n = length(zp);
  slope = zeros(3,size(xc,2));
  x = xc(1,:)/aperture;
  y = xc(2,:)/aperture;
  for i = 1:n
    for j = 1:length(zp{i})
      slope(1,:) = slope(1,:) + zp{i}(j) * eval(zfdx{i}{j});
      slope(2,:) = slope(2,:) + zp{i}(j) * eval(zfdy{i}{j});
    end
  end
  slope(3,:) = -1;
end
  