function n = vimbufsize
% VIMBUFSIZE: return the current size of the vimage buffer
% Syntax:
%   n = vimbufsize
% where n is the number of images that can be stored in the buffer.
%
% See also: VIMBUFRESIZE.
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_BUFFER
  
  n = length(VIMAGE_BUFFER);
