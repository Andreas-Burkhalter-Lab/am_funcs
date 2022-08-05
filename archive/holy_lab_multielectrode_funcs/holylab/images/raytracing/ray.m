function r = ray
% RAY: a structure defining a ray
% A ray structure has the following fields:
%   x0: position in space (2-d or 3-d) at the starting point
%   e: a unit vector indicating the direction of propagation
%   valid: will be true if a ray is still "in play," i.e., has
%     successfully passed through all of the optical elements
%   w: wavelength, in nm
%   rgb: color used for rendering the ray on the screen
%   I: intensity (not used currently)
%
% Note x0 and e must be _column vectors_.
%
% Note the ray traces out the points x0+t*e.  Interactions with surfaces
% return the value of t at which contact occurs.
%
% See also: PLANAR, SPHERICAL, CYLINDRICAL.
  
  r = struct('x0',[0; 0],'e',[1; 0],'valid',true,'w',488,'I',1);
