classdef opticalQuadratic < opticalAperture
% opticalQuadratic: spherical, cylindrical, and planar optical surfaces
%
% Constructor syntax:
%   s = opticalQuadratic(c,...)
%   s = opticalQuadratic([c1 c2],...)
% c is the curvature (1/R). When c is a scalar, you get rotationally
% symmetric surfaces, i.e., spherical or planar surfaces. For a cylindrical
% lens, supply the curvature as a pair [c1 c2], where c1 is the curvature
% along the first intrinsic coordinate, and c2 is the curvature along the
% second instrinsic coordinate.
% Note: the sign convention is that curvature is positive when the center
% of curvature is to the "right" (farther down the optic axis) of the
% surface vertex.
% The remaining arguments of the constructor are passed on to
% opticalAperture; this class inherits all the methods of that class.
%
% Methods:
%   rbPost = s.trace(rb,npre,npost)
% where rb is ray bundle (see opticalRayBundle) describing the
% pre-refraction rays; rbPost contains the output rays, which have been
% propagated to the surface and refracted in a new direction.
% npre and npost are vectors containing the refractive index for each ray
% in the input and output side, respectively. See opticalRayBundle for
% information about allowed syntaxes.
%
%   [X,Y,Z] = s.surf3d
%   [X,Y,Z] = s.surf3d(ngrid)
% Provides the coordinates of the surface, relative to the origin at
% the vertex.
% The second syntax lets you control the number of points along each axis of
% the grid used to represent the surface (default: 31).
%
%   [y,z] = s.surf2d(perp)
%   [y,z] = s.surf2d(perp,ngrid)
% Coordinates of the intersection of the surface in the plane perpendicular
% to the vector "perp".
% Here the default ngrid is 101.
%
%   w = s.span_along_optic_axis
% Gives the surfaces maximum displacement on either side of the vertex.
%
%   sflip = s.flip;
% Flips the orientation of the surface along its optic axis.
%
% See also: opticalAperture, opticalRayBundle.

% Copyright 2010 by Timothy E. Holy

  properties (GetAccess = public,SetAccess = private)
    c = [0 0]';
  end
  properties (Access = private)
    cz = [];  % will hold the component of c with max absolute value
  end
  methods
    %% Constructor
    function s = opticalQuadratic(cIn,varargin)
      % aperture
      s = s@opticalAperture(varargin{:});
      % c
      if isscalar(cIn)
        cIn = [1 1]*cIn;
      elseif (length(cIn) ~= 2)
        error('c must be a scalar or 2-vector, [c1 c2]');
      end
      s.c = cIn(:);
      if (abs(s.c(1)) > abs(s.c(2)))
        s.cz = s.c(1);
      else
        s.cz = s.c(2);
      end
      if any(abs(s.c) .* s.sz > 2)
        error(['The aperture exceeds the physical size of the cylinder,' ...
          ' did you forget to supply c as 1/R?']);
      end
    end
    
    %% Ray trace
    function rb = trace(s,rb,npre,npost)
      % Trace rays to aperture, and transform to relative coordinates
      [rb,p] = s.trace@opticalAperture(rb,npre);
      e = s.transform_rotate(rb.e);
      % Propagate from the aperture to the actual surface
      ctot = [s.c; s.cz];
      if (s.cz ~= 0)  % if planar, we are already at the surface
        % This corresponds to solving a quadratic
        coef2 = sum(bsxfun(@times,ctot,e.^2),1);
        coef1 = sum(bsxfun(@times,ctot,e.*p),1) - e(3,:);  % really half of coef1
        coef0 = sum(bsxfun(@times,ctot,p.^2),1) - 2*p(3,:); % last terms should be 0
        sqrtarg = coef1.^2 - coef0.*coef2;
        % Rays can pass through aperture, yet (if at steep angle) fail to
        % strike the surface. So we need to check for this.
        validFlag = sqrtarg >= 0;
        rb.valid = validFlag;
        coef2 = coef2(validFlag);
        coef1 = coef1(validFlag);
        coef0 = coef0(validFlag);
        % Finally, solve the quadratic
        t = solve_quadratic(coef2,2*coef1,coef0,-1);
        % Propagate the rays (it's OK if t is negative)
        rb = rb.propagate(t,npre);
        % Do it manually too on the relative coordinates
        p(:,validFlag) = p(:,validFlag) + bsxfun(@times,t,e(:,validFlag));
        p = p(:,validFlag);
      end
      % Refract the "new" ray bundle; first calculate the normal
      normal = bsxfun(@times,ctot,p);
      normal(3,:) = normal(3,:) - 1;
      normal = s.transform_rotate_inverse(normal);  % put in absolute coords
      rb = rb.refract(normal,npre./npost);
    end
    
    %% Surface geometry in 3d
    function [X,Y,Z] = surf3d(s,ngrid)
      if (nargin < 2)
        ngrid = 31;
      end
      unitgrid = linspace(-0.5,0.5,ngrid);
      if isscalar(s.sz)
        x = unitgrid*s.sz;
        [X,Y] = ndgrid(x,x);
        R2 = X.^2 + Y.^2;
        killFlag = R2 > s.sz^2/4;
        X(killFlag) = NaN;
        Y(killFlag) = NaN;
      else
        x = unitgrid*s.sz(1);
        y = unitgrid*s.sz(2);
        [X,Y] = ndgrid(x,y);
      end
      coef2 = repmat(s.cz,1,ngrid^2);
      coef1 = repmat(-2,1,ngrid^2);
      coef0 = (s.c(1)*X(:).^2 + s.c(2)*Y(:).^2)';
      Z = solve_quadratic(coef2,coef1,coef0,-1);
      XYZ = [X(:) Y(:) Z(:)]';
      XYZ = s.transform_shift_and_rotate_inverse(XYZ);
      X = reshape(XYZ(1,:),size(X));
      Y = reshape(XYZ(2,:),size(X));
      Z = reshape(XYZ(3,:),size(X));
    end
    
    %% Surface geometry in 2d
    function [y,z] = surf2d(s,perp,ngrid)
      if (nargin < 3)
        ngrid = 101;
      end
      % Find the vector that is perpendicular to both perp and the optic
      % axis
      ax = cross(perp,s.optic_axis);
      % Convert ax into relative coordinates and extract part that is
      % orthogonal to optic axis
      axrel = s.transform_rotate(ax);
      axnorm = sum(axrel.^2);
      if (axnorm == 0)
        error('Cannot project down to 2 dimensions along optic axis');
      end
      axrel = axrel/sqrt(axnorm);
      unitgrid = linspace(-0.5,0.5,ngrid);
      if isscalar(s.sz)
        l = unitgrid*s.sz;
      else
        L = sqrt(sum((axrel(1:2).*s.sz)^2));
        l = unitgrid*L;
      end
      S = axrel*l;
      coef2 = repmat(s.cz,1,ngrid);
      coef1 = repmat(-2,1,ngrid);
      coef0 = s.c(1)*S(1,:).^2 + s.c(2)*S(2,:).^2;
      S(3,:) = solve_quadratic(coef2,coef1,coef0,-1);
      S = s.transform_shift_and_rotate_inverse(S);
      z = sum(bsxfun(@times,S,s.optic_axis));
      y = sum(bsxfun(@times,S,ax));
    end
    
    %% span_along_optic_axis
    function w = span_along_optic_axis(s)
      if isscalar(s.sz)
        w = sort([0 solve_quadratic(s.cz,-2,s.cz*s.sz^2/4,-1)]);
      else
        w = sort([0 solve_quadratic(s.cz,-2,(s.c'*s.sz.^2)/4,-1)]);
      end
    end
    
    %% Flip
    function s = flip(s)
      s.c = -s.c;
      s.cz = -s.cz;
    end
  end
end
