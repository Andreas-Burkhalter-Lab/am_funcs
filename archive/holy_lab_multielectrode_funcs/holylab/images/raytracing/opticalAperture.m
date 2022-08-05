classdef opticalAperture < coordinateTransform3dRigid
% opticalAperture: location, orientation, and aperture information for optical components
%
% This is the base class for all optical components. It handles the
% geometry and transformations as well as tracing of rays to the aperture
% plane. Components that inherit from the class should be defined assuming
% that the optic axis is the z-axis (i.e., [0 0 1]); this class can handle
% all geometric transforms when that is not the case.
%
% Constructor syntax:
%   oa = opticalAperture(diameter)
%   oa = opticalAperture([l1 l2])
% The first syntax is for elements with a circular aperture (diameter must
% be a scalar), and the second is for elements with a rectangular aperture.
%
%   oa = opticalAperture(...,center_coordinates)
% This sets the position of the center of the aperture in 3-dimensional
% space.
%
%   oa = opticalAperture(...,center_coordinates,optic_axis)
% This additionally sets the direction of the optic axis of the component.
% The aperture is in the plane perpendicular to the optic axis.
%
% Methods:
%   rbout = oa.trace(rb,n)
%   [rbout,prel] = oa.trace(rb,n)
%   rbout = oa.trace(rb,n,npost)
% Trace a set of rays to their intersection with the aperture plane,
% through a medium of refractive index n (supplied as a scalar or 1 per
% ray). rb is of type opticalRayBundle. Rays that do not pass through the
% aperture are set to be invalid.
% The second syntax also returns prel, the position of the valid rays
% expressed in relative coordinates.
% The third syntax is present for consistency with tracing across surfaces
% separating media of different refractive indices; npost is not used at
% all.
%
% Methods for coordinate manipulations:
%   oa.optic_axis = v
% Sets the optic axis to be parallel to v, as a rotation from [0 0 1] along
% the axis perpendicular to both
%
%   oa = oa.spin(theta)
% Spin the component on its axis
%
%   oaOut = oa.rotate(q)
% Rotate the component by a quaternion q
%
% opticalAperture inherits from coordinateTransform3dRigid, so any
% methods defined for that class are applicable here as well. In
% particular,
%   prelative = oa.transform_shift_and_rotate(p)
% expresses a set of positions (e.g., ray starting positions) in terms of
% coordinates relative to the aperture and its optic axis, and
%   erelative = oa.transform_rotate(e)
% expresses a set of ray vectors (direction of propagation) in relative
% coordinates. There are also inverse transforms available.
%
%   [X,Y,Z] = oa.surf3d
%   [X,Y,Z] = oa.surf3d(ngrid)
%   [y,z] = s.surf2d(perp)
%   [y,z] = s.surf2d(perp,ngrid)
% For now, these return empty arrays.
%
% See also: coordinateTransform3dRigid.

% Copyright 2010 by Timothy E. Holy

  properties (SetAccess = public, GetAccess = public)
    sz = [];
  end
  properties (Dependent = true, Access = public)
    optic_axis;
  end
  methods
    %% Constructor
    function oa = opticalAperture(szIn,p,optax)
      if (nargin > 0)
        if (length(szIn) < 3 && ~isempty(szIn))
          oa.sz = szIn(:);
        else
          error(['sz must either be a scalar (for round apertures) or a' ...
            ' 2-vector (for rectangular apertures)']);
        end
        if (nargin > 1 && ~isempty(p))
          oa.origin = p(:);
        end
        if (nargin > 2 && ~isempty(optax))
          oa.optic_axis = optax;
        end
      end
    end % constructor
    
    %% Set & get methods
    function oa = set.optic_axis(oa,optax)
      if (length(optax) ~= 3)
        error('optic_axis must be a 3-vector');
      end
      optax_norm = sum(optax.^2);
      if (optax_norm == 0)
        error('optic_axis must not have zero length');
      end
      optax = optax(:) / sqrt(optax_norm);
      % Find the quaternion that will rotate this optic axis to [0 0 1]
      oa = oa.vecs2rotation(optax,[0 0 1]');
    end
    
    function optax = get.optic_axis(oa)
      optax = rotate_quaternion([0 0 1]',oa.q);
    end
    
    %% rotate
    function oa = rotate(oa,qIn)
      oa.q = multiply_quaternion(qIn,oa.q);
    end
    
    %% spin
    function oa = spin(oa,theta)
      optax = oa.optic_axis;
      qspin = [cos(theta/2); sin(theta/2)*optax];
      oa = oa.rotate(qspin);
    end
    
    %% trace
    function [rb,pout] = trace(oa,rb,n,~)
      p = rb.p;
      e = rb.e;
      % Calculate the distance needed to travel
      dp = bsxfun(@minus,p,oa.origin);
      dp_dot_n = sum(bsxfun(@times,dp,oa.optic_axis));
      e_dot_n = sum(bsxfun(@times,e,oa.optic_axis));
      goodFlag = e_dot_n ~= 0;
      rb.valid = goodFlag;
      t = -dp_dot_n(goodFlag)./e_dot_n(goodFlag);
      % Advance the ray
      [rb,p] = rb.propagate(t,n);
      % Determine which rays make it through the aperture
      p = oa.transform_shift_and_rotate(p);
      if (length(oa.sz) == 1)
        % circular aperture
        r2 = sum(p(1:2,:).^2,1);
        validFlag = r2 <= (oa.sz/2)^2;
      else
        % rectangular aperture
        validFlag = all(bsxfun(@le,abs(p(1:2,:)),oa.sz/2),1);
      end
      rb.valid = validFlag;
      if (nargout > 1)
        pout = p(:,validFlag);
      end
    end % trace
    
    %% Functions for compatibility with derived classes (virtual methods)
    function [X,Y,Z] = surf3d(~,~)
      X = [];
      Y = [];
      Z = [];
    end
    
    function [y,z] = surf2d(~,~,~)
      y = [];
      z = [];
    end
    
    function w = span_along_optic_axis(~)
      w = [0 0];
    end
    
    function oa = flip(oa)
    end
  end % methods
end
