function clearimagehandle(varargin)
% VIMAGE/CLEARIMAGEHANDLE: forget the image handle for a vimage
%
% This erases the imagehandle data from a vimage. Note that by default,
% SETIMAGEHANDLE causes this function to be called automatically when the
% image object is deleted.
%
% Syntax:
%   clearimagehandle(vim)
% where
%   vim is the vimage
%
% See also: VIMAGE/SETIMAGEHANDLE, VIMAGE/EVAL.
  
  global VIMAGE_LIST

  if (nargin == 0)
    [VIMAGE_LIST.imagehandle] = [];    % Clear them all
    return;
  elseif (nargin == 1)
    vim = varargin{1};
  else
    vim = varargin{3};   % Callback mode
  end
  index = vim.index;
  if (index <= length(VIMAGE_LIST))   % Make sure the vimageobject exists
    if (isempty(gcbo) || ...          % Not callback mode
        VIMAGE_LIST(index).imagehandle == gcbo)    % Callback mode, sane
      VIMAGE_LIST(index).imagehandle = [];
    else
      % It's a callback, and it's not the object we expected; issue
      % warning
      warning('vimage:clearimagehandle',...
              'imagehandle doesn''t match, not deleted.');
    end
  end
  