classdef coordinateTransform3dRigid
% coordinateTransform3dRigid: translations and rotations in 3d space
% Constructor syntax:
%   ct = coordinateTransform3dRigid(origin,q)
% Define in terms of the origin and associated quaternion. Coordinates x'
% in the transformed space are related to coordinates in the original space
% by
%   x' = q*(x-origin)*qinv
% (See Wikipedia page on "quaternions and spatial rotation" for an
% explanation.)
%
% You can also set these values directly by
%   ct.origin = p;
%   ct.q = quat;
% The value of q will be normalized, even if quat is not.
%
% Methods:
%   ct = ct.vecs2rotation(e1,e2)
% This sets q to be the quaternion that rotates the unit vector e1 to the
% unit vector e2. You must ensure that e1 and e2 are normalized.
%
%   xrel = ct.transform_shift_and_rotate(xabs)
%   vrel = ct.transform_rotate(vabs)
% These perform the operations indicated by their name, to put "absolute"
% coordinates xabs into the coordinate system specified by ct (in the
% formula above, x = xabs and x' = xrel).  For "positions" x you probably
% want to call transform_shift_and_rotate, for "tangent vectors" v or
% derivatives you probably want to call transform_rotate.
%   
%   xabs = ct.transform_shift_and_rotate_inverse(xrel)
%   vabs = ct.transform_rotate_inverse(vrel)
% These perform the inverse transforms of the ones above.
%
% Copyright 2010 by Timothy E. Holy

  properties (Access = public)
    origin = [0 0 0]';
    q = [1 0 0 0]';
  end
  methods
    %% Constructor
    function ct = coordinateTransform3dRigid(p,qIn)
      if (nargin > 0)
        ct.origin = p;
        if (nargin > 1)
          ct.q = qIn;
        end
      end
    end
    %% Set functions
    function ct = set.origin(ct,p)
      if (length(p) ~= 3)
        error('The origin must be a 3-vector');
      end
      ct.origin = p(:);
    end
    function ct = set.q(ct,qIn)
      if (length(qIn) ~= 4)
        error('q must be a 4-vector');
      end
      % Normalize
      ct.q = qIn(:) / sqrt(sum(qIn.^2));
    end
    %% Alternative methods of definition
    function ct = vecs2rotation(ct,e1,e2)
      quat = multiply_quaternion([0; e1(:)],[0; e2(:)]);
      ct.q = [sqrt((1+quat(1))/2);sqrt((1-quat(1))/2)*quat(2:end)];
    end
    %% Coordinate manipulation
    function ct = change_coordinates(ct,ctNew)
      ct.origin = ctNew.transform_shift_and_rotate(ct.origin);
      ct.q = multiply_quaternion(ctNew.q,ct.q);
    end
    %% Coordinate transformation functions
    function xrel = transform_shift_and_rotate(ct,xabs)
      xrel = bsxfun(@minus,xabs,ct.origin);
      if ~isequal(ct.q,[1 0 0 0]')
        xrel = rotate_quaternion(xrel,ct.q);
      end
    end
    function xrel = transform_rotate(ct,xabs)
      xrel = xabs;
      if ~isequal(ct.q,[1 0 0 0]')
        xrel = rotate_quaternion(xrel,ct.q);
      end
    end
    function xabs = transform_shift_and_rotate_inverse(ct,xrel)
      xabs = ct.transform_rotate_inverse(xrel);
      xabs = bsxfun(@plus,xabs,ct.origin);
    end
    function xabs = transform_rotate_inverse(ct,xrel)
      xabs = xrel;
      if ~isequal(ct.q,[1 0 0 0]')
        qconj = [ct.q(1); -ct.q(2:4)];
        xabs = rotate_quaternion(xabs,qconj);
      end
    end
  end
end
