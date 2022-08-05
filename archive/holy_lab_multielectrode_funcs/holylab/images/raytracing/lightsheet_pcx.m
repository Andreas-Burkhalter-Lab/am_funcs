function [c,r] = lightsheet_pcx(s)
% lightsheet_pcx: light sheet illumination using plano-convex collimator
% Syntax:
%   [c,r] = lightsheet_pcx(s)
% where s is a structure with the following fields:
%   s.colR = radius of curvature of collimating lens
%   s.colCT = center thickness of collimating lens
%   s.colPosition = 2-vector giving position of planar surface of
%     collimating lens ([x y], i.e., [alongaxis perpaxis]) relative to
%     the source location at [0 0]. If instead this is supplied as a
%     scalar, it specifies y = 0.
%   s.colangle = angle between optic axis of source and optic axis of
%     lens (0 = perfect alignment)
%   s.colmat = string representing collimating lens material
%
%   s.cylR = radius of curvature of cylindrical lens
%   s.cylCT = center thickness of cylindrical lens
%   s.cylPosition = 2-vector giving position of planar surface of
%     cylindrical lens ([x y], i.e., [alongaxis perpaxis]) relative to
%     source location. If a scalar, y = 0.
%   s.cylangle = angle between optic axis of source and optic axis of
%     lens (default: 0, meaning perfect alignment)
%   s.cylmat = string representing cylindrical lens material
%
%   s.immerse = string representing the material ('air' or 'saline') in
%     the back surface of the cylindrical lens is immersed.
%
%   s.fiberna = the numerical aperture of the fiber emission
%   s.nrays = the number of rays to use in tracing (default: 41)
%   s.w = the wavelength (in nanometers) to use for the tracing (default: 488)
%
% On output,
%   c is a cell array of components;
%   r is a ray structure array.
  
% Copyright 2004-2010 by Timothy E. Holy

  s = default(s,'nrays',41,'w',488,'cylangle',0);
  if (length(s.colPosition) == 1)
    s.colPosition(2) = 0;
  end
  if (length(s.cylPosition) == 1)
    s.cylPosition(2) = 0;
  end
  

  % Create the input rays
  r0.x0 = [0, 0];  % source is at 0,0
  r0.w = s.w;
  r0.e = [0 0];
  r0.I = 1;
  theta = linspace(-asin(s.fiberna),asin(s.fiberna),s.nrays);
  r = r0;
  for i = 1:length(theta)
    r(i) = r0;
    r(i).e = [cos(theta(i)) sin(theta(i))];
    r(i).I = 1; %exp(-theta(i)^2/(2*s.fiberna^2));
  end
  
  
  % Collimating lens
  col_rot = [cos(s.colangle) sin(s.colangle);...
             -sin(s.colangle) cos(s.colangle)];   % rotation matrix
  c = pcx(s.colPosition,(col_rot*[1;0])',s.colR,s.colCT,s.colmat,'air','air');
  
  % Cylindrical lens
  cyl_rot = [cos(s.cylangle) sin(s.cylangle);...
             -sin(s.cylangle) cos(s.cylangle)];   % rotation matrix
  ctmp = pcx(s.cylPosition,(cyl_rot*[-1;0])',s.cylR,s.cylCT,s.cylmat,s.immerse,'air');
  c = [c ctmp([2 1])];
  
  % Project onto a screen
  % Now an extra planar surface to act as the "screen"
  screenpos = s.cylPosition(1)+10;
  p = struct('c',screenpos,...
             'normal',[1 0],...
             'mat1',s.immerse,...
             'mat2',s.immerse,...
             'apc',[screenpos 0],...
             'apr',3);
  c{5} = {@opt2dline,p};
  
