function vimbufresize(n)
% VIMBUFRESIZE: resize the vimage buffer
% Syntax:
%   vimbufresize(n)
% where n is the number of images you want to store on the buffer. Note
% this function clears the buffer when executed.
%
% See also: VIMBUFSIZE.
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_BUFFER VIMAGE_BUFFERPOS
  
  if ~isscalar(n)
    error('Input must be a scalar');
  end
  VIMAGE_BUFFER = repmat(vimagebufferobj,1,n);
  VIMAGE_BUFFERPOS = 1;
