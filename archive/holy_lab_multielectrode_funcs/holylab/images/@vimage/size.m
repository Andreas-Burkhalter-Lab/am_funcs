function [szo,varargout] = size(vim,pos)
% VIMAGE/SIZE: obtain the size of a vimage (once evaluated)
% This function calculates size quickly, without having to evaluate the
% vimage.
%
% Syntax:
%   sz = size(vim)
% where vim is the vimage.
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE

  if (length(vim) > 1)
    error('Must be a single vimage');
  end
  index = vim.index;
  tmp = VIMAGE_LIST(index);
  
  % First try direct evaluation
  if ~isempty(tmp.image)
    sz = size(tmp.image);
  else
    % Calculate from the function used to obtain this image
    switch tmp.func
     case 'imload'
      sz = tmp.argin{1}.size;
     case 'imresize'
      sz = floor(tmp.argin{2}*size(tmp.argin{1}));
     case 'imcrop'
      sz = tmp.argin{2}(3:4);
     otherwise
      % Find ancestor = first vimage in arguments
      index = 1;
      while ~isa(tmp.argin{index},'vimage')
        index = index+1;  % This shouldn't overflow, because we're in
                          % vimage/size 
      end
      if (nargin == 1)
        [szo,varargout{1:nargout-1}] = size(tmp.argin{index});
      else
        [szo,varargout{1:nargout-1}] = size(tmp.argin{index},pos);
      end
      return;
    end
  end
  
  % Now massage it into the correct output format
  if (nargin > 1)
    sz = sz(pos);
  end
  nout = nargout;
  if (nout == 0)
    nout = 1;
  end
  if (nout == 1)
    tmpout{1} = sz;
  else
    for i = 1:min(nout,length(sz))
      tmpout{i} = sz(i);
    end
    if (nout < length(sz))
      tmpout{nout} = prod(sz(nout:end));
    elseif (nout > length(sz))
      [tmpout{length(sz)+1:nargout}] = deal(1);
    end
  end
  szo = tmpout{1};
  varargout = tmpout(2:end);
