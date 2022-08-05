function im = eval(vim)
% VIMAGE/EVAL: evaluate a virtual image
% Syntax:
%   im = eval(vim)
% where vim is a virtual image, and im is a structure with im.image
% containing the image data.

% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE

  if (length(vim) > 1)
    error('Must be a single vimage');
  end
  index = vim.index;
  tmp = VIMAGE_LIST(index);
  
  % First try direct evaluation
  im = tmp.image;
  if ~isempty(im)
    return
  end

  % Now look on the buffer
  bufindex = find([VIMAGE_BUFFER.index] == index);
  if ~isempty(bufindex)
    im = VIMAGE_BUFFER(bufindex).image;
    return;
  end

  % Next, try the image handle
  if ~isempty(tmp.imagehandle)
    if ishandle(tmp.imagehandle)
      im = get(tmp.imagehandle,'CData');
      return;
    else
      % Handle expired, clear it
      VIMAGE_LIST(index).imagehandle = [];
    end
  end

  % OK, we really do have to work for this image
  % Calculate all the virtual images that this one depends upon
  argin = tmp.argin;
  for i = 1:length(argin)
    if isa(argin{i},'vimage')
      argin{i} = eval(argin{i});  % Evaluate the dependencies
    end
  end
  % Now we have all the values we need, we can execute the function
  im = feval(tmp.func,argin{:});
  % Are we profiling? If so, increase the counter
  if (VIMAGE_PROFILE)
    VIMAGE_LIST(index).count = VIMAGE_LIST(index).count + 1;
  end
  % Put it on the circular buffer
  VIMAGE_BUFFER(VIMAGE_BUFFERPOS).index = index;
  VIMAGE_BUFFER(VIMAGE_BUFFERPOS).image = im;
  VIMAGE_BUFFERPOS = VIMAGE_BUFFERPOS + 1;
  if (VIMAGE_BUFFERPOS > length(VIMAGE_BUFFER))
    VIMAGE_BUFFERPOS = 1;
  end
