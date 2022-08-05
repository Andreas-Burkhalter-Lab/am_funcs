function setimagehandle(vim,himg,clearflag)
% VIMAGE/SETIMAGEHANDLE: specify image handle for a vimage
%
% This allows image data to be obtained from the CData property of an
% image, providing another fast route to vimage evaluation.
%
% Syntax:
%   setimagehandle(vim,himg)
% where
%   vim is the vimage
%   himg is a handle to a graphics object of type 'image'
%
% By default, this function sets the 'DeleteFcn' property of the
% image. When the image graphics object is deleted, it will clear the
% imagehandle from the vimage. You can prevent this from happening with
% the syntax
%
%   setimagehandle(vim,himg,'noautoclear')
%
% Note that EVAL will clear the imagehandle field if it discovers
% that the handle is no longer valid.
%
% Note also that VIMAGELOAD clears all the imagehandles upon loading, so
% there is no concern of confusion across Matlab sessions.
%
% See also: VIMAGE/CLEARIMAGEHANDLE, VIMAGE/EVAL.
  
  global VIMAGE_LIST

  index = vim.index;
  if (ishandle(himg) && strcmp(get(himg,'Type'),'image'))
    VIMAGE_LIST(index).imagehandle = himg;
    % Set up automatic clearing
    if (nargin < 3 || ~ischar(clearflag) || ...
        ~strcmp(clearflag,'noautoclear'))
      set(himg,'DeleteFcn',{@clearimagehandle,vim});
    end
  end
  