function unlock(vim)
% VIMAGE/UNLOCK: free locked image data from memory
%
% Syntax:
%   unlock(vim)
% where vim is vimage.
%
% See also: VIMAGE/LOCK, VIMAGE/ISLOCKED.
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST

  VIMAGE_LIST(vim.index).image = [];
  