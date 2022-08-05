function s = islocked(vim)
% VIMAGE/ISLOCKED: determine whether there is locked image data
%
% Syntax:
%   s = islocked(vim)
% where vim is the vimage.
%
% See also: VIMAGE/LOCK, VIMAGE/UNLOCK.
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST

  s = ~isempty(VIMAGE_LIST(vim.index).image);
  