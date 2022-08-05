function cl = imclass(vim)
% VIMAGE/IMCLASS
% Syntax
%   cl = imclass(vim)
% where cl is a string ('uint16', 'single', etc)
  
% Copyright 2005 by Timothy E. Holy
  
  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE

  if (length(vim) > 1)
    error('Must be a single vimage');
  end
  index = vim.index;
  tmp = VIMAGE_LIST(index);
  
  % First try direct evaluation
  if ~isempty(tmp.image)
    cl = class(tmp.image);
    return
  end

%   % Now look on the buffer
%   % Perhaps should skip this one, since it will often fail,
%   % and last method is also fast?
%   bufindex = find([VIMAGE_BUFFER.index] == index);
%   if ~isempty(bufindex)
%     cl = class(VIMAGE_BUFFER(bufindex).image);
%     return
%   end
% 
  % Don't do image handle because data are changed to double

  % Calculate from the function used to obtain this image
  switch tmp.func
   case 'imload'
    cl = tmp.argin{1}.prec;
   case {'immean','imwsum','single'}
    cl = 'single';
   otherwise
    % Find ancestor = first vimage in arguments
    index = 1;
    while ~isa(tmp.argin{index},'vimage')
      index = index+1;  % This shouldn't overflow, because we're in vimage/size
    end
    cl = imclass(tmp.argin{index});
  end
