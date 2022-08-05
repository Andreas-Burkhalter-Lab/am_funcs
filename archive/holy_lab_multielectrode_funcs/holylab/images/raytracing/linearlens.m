function [x1,s1] = linearlens(x0,s0,xn0,xn1,z0,z1,n0,n1)
% LINEARLENS: an ideal (sine-condition satisfying) lens model
% This implements the linear model in the Holy manuscript.  For the
% linear model, this is more flexible than RAYTRACE_IDEAL_LENS, which
% assumes the optic axis is along the z-axis.
%
% Syntax:
%   [x1,s1] = linearlens(x0,s0,xn0,xn1,z0,z1,n0,n1)
% where 0 means the input size, 1 means the output side and
%   x: ray positions, d-by-n_rays matrix (d = dimensionality)
%   s: propagation directions (d-by-n_rays matrix, each column is normalized)
%   xn: positions of the axial nodal points, d-by-1 vectors.  The vector
%     between these points defines the optic axis (so there must be some
%     gap, even if "eps", between these two even if you want a "thin
%     lens").
%   z: distances to focal surfaces (can be inf for a infinite-conjugate
%     lens).  Be sure to follow the appropriate sign convention;
%     distances are measured from nodal points.
%   n: refractive indices of the media in which the object & image are
%     immersed.
%
% See also: RAYTRACE_IDEAL_LENS.

  if isinf(z0) && isinf(z1)
    error('Telescopes not implemented')
  end
  
  [d,n_rays] = size(x0);  % d = dimensionality
  
  % Calculate the unit vector defining the optic axis
  optic_axis = xn1-xn0;
  optic_axis = optic_axis / sqrt(sum(optic_axis.^2));
  % Calculate the projection matrix that eliminates the component along
  % the optic axis
  P = eye(d,d) - optic_axis * optic_axis';
    
  if ~isinf(z0)
    % The ray(s) enter from a finite-conjugate side. Trace to focal
    % plane.
    xf0 = xn0 + z0*optic_axis;  % the focal point
    dx = x0 - repmat(xf0,1,n_rays);
    t = (optic_axis' * dx) ./ (optic_axis' * s0);
    x0f = x0 - repmat(t,d,1) .* s0;  % position in focal plane
    % Calculate the projected position and the nodal ray
    x0fp = P*x0f;
    sNp = x0fp / z0;  % the nodal ray for each object point
    if isinf(z1)
      % The other side is infinite conjugate trace through lens
      yp = (n0*z0/n1)*(P*s0 - sNp);  % position of output relative to nodal pt
      x1 = repmat(xn1,1,n_rays) + yp;
      t1 = sqrt(1-sum(sNp.^2));
      s1 = sNp + repmat(t1,d,1) .* optic_axis;
    else
      % The other side is also finite-conjugate
      M = z1/z0;  % magnification
      Mprime = M * n1/n0;  % reduced magnification
      x1fp = M*x0fp;
      x1 = x1fp + repmat(xn1 + z1*optic_axis,1,n_rays);
      s1p = (P*s0 - sNp)/Mprime + sNp;
      t1 = sqrt(1-sum(s1p.^2));
      s1 = s1p + optic_axis * t1;
      % Back trace to (rougly) nodal plane
      x1 = x1 - z1*s1;
    end
  else
    % The ray(s) enter from the infinite-conjugate side. Trace to nodal
    % plane.
    dx = x0 - repmat(xn0,1,n_rays);
    t = (optic_axis' * dx) ./ (optic_axis' * s0);
    y = x0 + repmat(t,d,1) .* s0;  % position in nodal plane
    % Trace through lens
    sNp = P*s0;
    x1 = z1 * (repmat(optic_axis,1,n_rays) + sNp) + repmat(xn1,1,n_rays);
    s1p = (n0/(n1*z1))*(P*y) + sNp;
    t1 = sqrt(1-sum(s1p.^2));
    s1 = s1p + repmat(t1,d,1) .* optic_axis;
    % translate back to exit nodal plane, a good proxy for the "exit from
    % the lens"
    t2 = (optic_axis'*(x1 - x1n))./(optic_axis'*s1);
    x1 = x1 + repmat(t2,d,1).*s1;
  end
end
  
    
    