function cOut = optical_shift_component(cIn,p)
% optical_shift_component: change the position of a component
%
% Syntax:
%   cOut = optical_shift_component(cIn,p)
% where cIn and cOut are cell arrays specifying the component layout, and
% p is a 3-vector of the displacement. 
%
% See also: optical_flip_component, optical_rotate_component,
% optical_spin_component.

% Copyright 2010 by Timothy E. Holy

  if iscell(cIn)
    cOut = cIn;
    oaFlag = cellfun(@(o) isa(o,'opticalAperture'),cOut);
    oaIndex = find(oaFlag);
    if isempty(oaIndex)
      return
    end
    for i = oaIndex
      cOut{i}.origin = cOut{i}.origin + p(:);
    end
  else
    cOut = cIn;
    if isa(cIn,'opticalAperture')
      cOut.origin = cOut.origin + p(:);
    end
  end
  