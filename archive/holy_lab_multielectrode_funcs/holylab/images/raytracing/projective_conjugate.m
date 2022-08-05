function xo = projective_conjugate(params,x,entry_index)
% Calculate the conjugate point for a projective transformation
% Syntax:
%    xconj = projective_conjugate(params,x)
% 
% params is the parameters structure describing the projective
% transformation, of the format described in OPT2DPROJECTIVE (only
% focalpoints and f need to be specified).
%
% x is a 2-vector, or a 2-by-npts matrix, of input points.
%
% xo is a 2-vector, or a 2-by-npts matrix, containing the coordinates of
% the point(s) conjugate to x. 
%
% See also: OPT2DPROJECTIVE, PROJECTIVE_COMBINE.

  f = params.f;
  if (length(f) == 1)
    f(2) = -f(1);
  end
  ff = f(1)*f(2);
  
  fpvec = diff(params.focalpoints,1,2);
  fpsep = norm(fpvec);
  %normal = fpvec/fpsep; % unit vector along optic axis
  normal = params.normal;
  N = diag(normal);
  R = opt2p_rotmtrx(normal);
  
  if (size(x,1) == 1)
    x = x(:);
  end
  n_pts = size(x,2);

  % Determine which side the points are on, by comparing their separation
  % from the unit points with the sign of the focal lengths
  if (nargin < 3)
    entry_index = ones(1,n_pts);
  end
%   unit_points = params.focalpoints + normal*f;
%   entry_index = nan(1,n_pts);
%   for pointIndex = 1:n_pts
%     dx = unit_points - x(:,[pointIndex pointIndex]);
%     side_indicator = f.*sum(N*dx);
%     entry_index(pointIndex) = find(side_indicator > 0);
%   end

  % Calculate the consequences of the projective transform using the
  % Newton formulation 
  dx = x - params.focalpoints(:,entry_index);
  dx = R*dx;
  z = dx(1);
  zp = ff/z;
  yp = f(1)*dx(2)/z;
  dxo = [zp;yp];
  dxo = R'*dxo; % rotate back to original coordinates
  xo = dxo + params.focalpoints(:,3-entry_index); % offset from opposite fp

function R = opt2p_rotmtrx(normal)
  R = [normal(1) normal(2); -normal(2) normal(1)]; % Rotation matrix
  