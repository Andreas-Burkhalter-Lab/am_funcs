function m = stack2zcol(s)
% STACK2ZCOL: convert an image stack to a z-column matrix
% Syntax:
%   m = stack2zcol(s)
% where
%   s is an image stack, of size [nx ny nz]
% and
%   m is a zcolumn matrix, of size [nz nx*ny]
%
% See also: ZCOL2STACK.
  
% Copyright 2006 by Timothy E. Holy
  
  sz = size(s);
  if (length(sz) ~= 3)
    error('Stack must be three-dimensional');
  end
  
  m = shiftdim(s,2);
  m = reshape(m,[sz(3) sz(1)*sz(2)]);
  
  