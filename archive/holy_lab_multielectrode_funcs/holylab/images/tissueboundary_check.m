function tissueboundary_check(stk,zi,options)
% TISSUEBOUNDARY_CHECK: verify that the tissue boundary has been calculated
% correctly
% Syntax:
%   tissueboundary_check(stk,zi,options)
% where
%   stk is the same stack passed to tissueboundary
%   zi is the interpolated tissue height (an output of tissueboundary).
%
%Possible fields of options:
%  pauseeach (default false): if true, will pause after each frame;
%  pausefirst (default false): if true, will pause after the first frame.
%
% This will open a figure window; by pressing a key you will step through
% the stack and see where the boundary was drawn.
%
% See also: TISSUEBOUNDARY.

% Copyright 2006 by Timothy E. Holy

%revision history - image dimensions and tissue orientation are rectified
%and figure plays through as a movie instead of pausing on each
%frame(5/9/06).TFH

  if (nargin <3)
      options = struct;
  end
  if ~isfield(options,'pauseeach')
      options.pauseeach = false;
  end
  if ~isfield(options,'pausefirst')
      options.pausefirst = false;
  end
  maxstk = max(stk(:));%makes the tissue boundary line pixel intensity value equivalent to the brightest pixel in the image.
  sz_stk = size(stk);
  for i = 1:prod(sz_stk(2:end))
    stk(round(zi(i)),i) = maxstk;
  end
  figure
  colormap(gray(256))
    if (ndims(stk) == 3)
    for i = 1:sz_stk(3);
      im = squeeze(stk(:,:,i));
      im = flipud(im);
      imagesc(im,[0 maxstk]);
      axis image
      title([num2str(i) ' of ' num2str(sz_stk(3))])
      if options.pauseeach || (i == 1 && options.pausefirst)
          pause
      else
          drawnow
      end
    end
  else
    imagesc(stk,[0 maxstk]);
  end