function c = symsph(xc,orient,R,thick,h,matl,mat1,mat2)
% SYMSPH: symmetric spherical lens
% Syntax:
%   c = symsph(xc,h,orient,R,thick,matl,mat1,mat2)
% where
%   xc is the location of the center of the planar face
%   h is the height of the lens (half the diameter)
%   orient is the optic axis
%   R is the radius of curvatures of the spherical face
%   thick is the center thickness of the lens (thickness along optic
%     axis)
%   matl is the lens material (e.g., 'bk7')
%   mat1 is the immersion medium on the planar surface (e.g., 'saline')
%   matc is the immersion medium on the spherical side (e.g., 'air')
  
  if (nargin < 7 || isempty(mat1))
    matp = 'air';
  end
  if (nargin < 8 || isempty(mat2))
    matc = 'air';
  end
  theta = asin(abs(h/R));
  % First curved surface
  p = struct('center',xc + (R-thick/2)*orient,...
             'Rvec',-orient*R,...
             'theta',theta,...
             'mat1',matl,...
             'mat2',mat1);
  c{1} = {@opt2dcircle,p};
  % Second curved surface
  p = struct('center',xc - (R-thick/2)*orient,...
             'Rvec',orient*R,...
             'theta',theta,...
             'mat1',matl,...
             'mat2',mat2);
  c{2} = {@opt2dcircle,p};
