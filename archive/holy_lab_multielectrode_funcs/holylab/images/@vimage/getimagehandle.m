function getimagehandle(vim,himg,clearflag)
% VIMAGE/GETIMAGEHANDLE: retrieve the image handle for a vimage
%
% Syntax:
%   himg = getimagehandle(vim)
% where
%   vim is the vimage
%   himg is a handle to a graphics object of type 'image'
%
% See also: VIMAGE/SETIMAGEHANDLE.
  
  global VIMAGE_LIST

  index = vim.index;
  himg = VIMAGE_LIST(index).imagehandle;
