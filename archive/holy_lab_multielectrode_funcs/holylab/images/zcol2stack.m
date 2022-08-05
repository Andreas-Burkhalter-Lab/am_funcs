function s = zcol2stack(m,sz)
% ZCOL2STACK: convert an image stack to a z-column matrix
% Syntax:
%   s = zcol2stack(m,sz)
% where
%   m is a zcolumn matrix, of size [nz nx*ny]
%   sz is either a 3-vector, [nx ny nz], or a 2-vector, [nx ny]
% and
%   s is an image stack, of size [nx ny nz]
%
% See also: STACK2ZCOL.
  
% Copyright 2006 by Timothy E. Holy
  
  szm = size(m);
  if (length(szm) ~= 2)
    error('zcol must be two-dimensional');
  end
  if (length(sz) == 3 && sz(3) ~= szm(1))
    error(['Mismatch in size of m and the specified stack ' ...
           'dimensionality']);
  end
  
  s = reshape(m,[szm(1) sz(1) sz(2)]);
  s = shiftdim(s,1);
