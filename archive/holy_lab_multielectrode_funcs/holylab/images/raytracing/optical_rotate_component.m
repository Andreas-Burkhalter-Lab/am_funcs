function cOut = optical_rotate_component(cIn,q)
% optical_rotate_component: rotate the optic axis of a component
%
% Syntax:
%   cOut = optical_rotate_component(cIn,q)
% where cIn and cOut are cell arrays specifying the component layout, and
% q is a quaternion for performing the rotation. For compound components,
% the component is rotated around its first axial point.
%
% See also: optical_shift_component, optical_flip_component,
% optical_spin_component.

% Copyright 2010 by Timothy E. Holy

  if iscell(cIn)
    cOut = cIn;
    oaFlag = cellfun(@(o) isa(o,'opticalAperture'),cOut);
    oaIndex = find(oaFlag);
    if isempty(oaIndex)
      return
    end
    p0 = cOut{oaIndex(1)}.origin;
    ct = coordinateTransform3dRigid(p0,q);
    for i = oaIndex
      cOut{i} = cOut{i}.change_coordinates(ct);
    end
  else
    cOut = cIn;
    if isa(cIn,'opticalAperture')
      p0 = cOut.origin;
      ct = coordinateTransform3dRigid(p0,q);
      cOut = cOut.change_coordinates(ct);
    end
  end
  