function lock(vim)
% VIMAGE/LOCK: lock a real image in memory
% This is useful if you will need a certain image repeatedly for a
% calculation.  A typical example might be if you want to subtract a
% single background image, and don't want to have to re-load it from disk
% every time.
%
% Syntax:
%   lock(vim)
% This evals the vimage "vim" and stores the result so that it will not
% disappear from the image buffer.
%
% See also: VIMAGE/UNLOCK, VIMAGE/ISLOCKED.
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST
  
  VIMAGE_LIST(vim.index).image = eval(vim);
  