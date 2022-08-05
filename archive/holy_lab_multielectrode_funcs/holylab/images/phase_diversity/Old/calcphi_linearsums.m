function s = calcphi_linearsums(A,offsets)
% CALCPHI_LINEARSUMS: manage linear combinations of phases
% 
% This is a utility function for calcphi2d, which returns the conversion
% functions param2phi (converting parameters to phases) and gphi2gparam
% (converting derivatives with respect to phases into derivatives with
% respect to parameters).  The particular model being represented here is
% a linear mixing of phase matrices,
%    phi1 = A11 * p1 + A12 * p2 + ... + A1n * pn + offset1
%    phi2 = A21 * p1 + A22 * p2 + ... + A2n * pn + offset2
%     ....
% where each of the pk and offsetk would be a suitable phase matrix on
% its own.
%
% For example, suppose you have a three images taken in which a
% single actuator is set to three different voltages, say [-0.1, 0,
% 0.1].  You want to model the phase for each image as an unknown
% constant aberration plus an aberration whose amplitude depends upon
% voltage.  In that case you might call this function with
%   A = [1 -0.1; 1 0; 1 0.1];
%   offsets = 0;
% The first column of A represents the constant unknown aberration, the
% second column the voltage-dependent component.
%
% As another example, suppose you have two images which differ by a known
% phase aberration, and you're wanting to solve for a common unknown
% phase aberration (this is the case treated in Paxman et al 1992).  Then
%   A = [1;1];
%   offsets = zeros(h,w,2); offsets(:,:,2) = known_phase_difference;
%
% Syntax:
%   s = calcphi_linearsums(A,offsets)
% where the fields of s are set to contain the two conversion functions
% needed for CALCPHI2D.
% 
% See also: CALCPHI2D, CALCPHI_ZERNIKE.
  
% Copyright 2009 by Timothy E. Holy
  
  if (nargin < 2)
    offsets = 0;
  end
  
  s.param2phi = @(p) remix(p,A',offsets);
  s.gphi2gparam = @(grad) remix(grad,A,0);
end

function Y = remix(X,A,offsets)
  szX = size(X);
  Xrs = reshape(X,[szX(1)*szX(2) szX(3:end)]);
  Yrs = Xrs * A;
  Y = reshape(Yrs,[szX(1) szX(2) size(A,2)]);
  Y = Y + offsets;
end
