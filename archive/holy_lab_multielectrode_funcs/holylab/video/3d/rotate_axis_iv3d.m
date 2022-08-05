function [zyz,str] = rotate_axis_iv3d(u)
% ROTATE_AXIS_IV3D: compute parameters for rotations around a particular axis in ImageVis3D
%
% In writing ImageVis3D scripts, one often wants to rotate around an
% arbitrary axis. This function computes the rotations in terms of XYZ
% rotations that achieves this.
%
% Syntax:
%   [zyz,str] = rotate_axis_iv3d(u)
% where
%   u is a vector around which you want to rotate
% and
%   zyz is a 3-vector that contains parameters for a 1-degree right-handed
%     rotation around u
%   str is a format string, used this way:
%       fprintf(fid,str,n_deg*zyz)
%     where fid is the file identifier for an open script file. 
%
% See also: T2X, X2T.

% Copyright 2010 by Timothy E. Holy

  u = u(:) / sqrt(sum(u.^2));
  zyz = t2x(x2t([0; 0; 0; 1; u; 1],'van'),'rpm');
  zyz = zyz(5:7);
  str = 'rotateZ %g\nrotateY %g\nrotateZ %g\n';
end
  