function cOut = optical_flip_component(cIn)
% optical_flip_component: flip a component along its optic axis
%
% This function reverses the order of surfaces and materials and flips
% the orientation of each surface.  For compound components, the
% "reflection" occurs around the midpoint.
%
% Syntax:
%   cOut = optical_flip_component(cIn)
% where cIn and cOut are cell arrays specifying the component layout.
%
% See also: optical_shift_component, optical_rotate_component.

% Copyright 2010 by Timothy E. Holy

  if iscell(cIn)
    cOut = cIn(end:-1:1);
    oaFlag = cellfun(@(o) isa(o,'opticalAperture'),cOut);
    oaIndex = find(oaFlag);
    if isempty(oaIndex)
      return
    end
%     p0 = (cOut{oaIndex(1)}.p + cOut{oaIndex(end)}.p)/2;
    p0 = (cOut{oaIndex(1)}.origin + cOut{oaIndex(end)}.origin)/2;
    for i = oaIndex
      cOut{i}.origin = 2*p0 - cOut{i}.origin;
      cOut{i} = cOut{i}.flip;
    end
  else
    cOut = cIn;
    if isa(cOut,'opticalAperture')
      cOut = cOut.flip;
    end
  end
  