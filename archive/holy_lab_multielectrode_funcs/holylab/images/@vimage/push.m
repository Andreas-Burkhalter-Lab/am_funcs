function push(vim,func,varargin)
% VIMAGE/PUSH: add a new method for obtaining an image
% Syntax:
%   push(vimage,func,args...)
% where
%   vimage is a virtual image (vimage)
%   func is a string with the name of the function
%   args are the arguments to be passed to the function

% Formerly but no longer works:
% Alternative syntax (internal use):
%   push(index,func,varargin)
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST

  tmp = VIMAGE_LIST(vim.index);
  if ~isempty(tmp.func)
    % There was already a method defined, put it on the history
    histmethod = rmfield(tmp,'image');       % Strip off any image data
    tmp.history = [tmp.history histmethod];  % Append to history
  end
  
  % Now fill in the new method
  tmp.func = func;
  tmp.argin = varargin;
  tmp.count = 0;
  
  % Put back on the list
  VIMAGE_LIST(vim.index) = tmp;
