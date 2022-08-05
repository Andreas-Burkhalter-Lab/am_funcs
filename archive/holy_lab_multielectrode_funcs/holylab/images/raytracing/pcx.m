function c = pcx(xp,orient,R,thick,matl,matp,matc)
% PCX: make a plano-convex lens
% Syntax:
%   c = pcx(xp,orient,R,thick,matl,matp,matc)
% where
%   xp is the location of the center of the planar face
%   orient is the optic axis, pointing from the center of the planar face
%     to the center of the spherical face
%   R is the radius of curvatures of the spherical face
%   thick is the center thickness of the lens (thickness along optic
%     axis)
%   matl is the lens material (e.g., 'bk7')
%   matp is the immersion medium on the planar surface (e.g., 'saline')
%   matc is the immersion medium on the spherical side (e.g., 'air')
  
  if (nargin < 6 || isempty(matp))
    matp = 'air';
  end
  if (nargin < 7 || isempty(matc))
    matc = 'air';
  end
  theta = acos(1-thick/R);
  r = R*sin(theta);
  % The planar surface of the PCX lens
  p = struct('center',xp,...
             'normal',orient,...
             'mat1',matp,...
             'mat2',matl,...
             'apr',r);
  c{1} = {@opt2dline,p};
  % Now do the curved surface of the PCX lens
  % The spherical surface of the PCX lens
  p = struct('center',xp - (R-thick)*orient,...
             'Rvec',orient*R,...
             'theta',theta,...
             'mat1',matl,...
             'mat2',matc);
  c{2} = {@opt2dcircle,p};
